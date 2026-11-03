"""Vertiport-location optimization models and stochastic search methods.

Implements the exact capacitated facility-location model, the SMA-inspired
binary search (SIBS), variable neighborhood search (VNS), and their hybrid for
the study scenarios, using IBGE 2022 population-grid data. The serialized method
identifiers ``SMA`` and ``SMA+VNS`` denote SIBS and the SIBS-VNS hybrid.
"""

import time
import os
import json
import hashlib
import platform
from pathlib import Path
from datetime import datetime, timezone
from typing import Tuple, Dict, Any, List
import warnings

import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from gurobipy import Model, GRB, quicksum


def collect_environment_metadata() -> Dict[str, Any]:
    '''Capture the software/hardware identity used by a runner.'''
    import sys
    meta = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'python_version': sys.version.split()[0],
        'platform_system': platform.system(),
        'platform_release': platform.release(),
        'platform_version': platform.version(),
        'platform_processor': platform.processor(),
        'platform_machine': platform.machine(),
        'cpu_count': os.cpu_count(),
        'numpy_version': np.__version__,
        'pandas_version': pd.__version__,
        'geopandas_version': gpd.__version__,
        'gurobipy_version': '.'.join(map(str, GRB.VERSION_MAJOR_MINOR_TECHNICAL))
        if hasattr(GRB, 'VERSION_MAJOR_MINOR_TECHNICAL') else None,
    }
    try:
        import scipy
        meta['scipy_version'] = scipy.__version__
    except Exception:
        pass
    try:
        import matplotlib
        meta['matplotlib_version'] = matplotlib.__version__
    except Exception:
        pass
    try:
        import gurobipy
        meta['gurobipy_version'] = '.'.join(map(str, gurobipy.gurobi.version()))
    except Exception:
        pass
    try:
        import psutil
        meta['ram_gb'] = round(psutil.virtual_memory().total / (1024 ** 3), 1)
    except Exception:
        meta['ram_gb'] = None
    return meta

# ============================================================
# Output paths
# ============================================================

# Output directory (same folder as this script)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'optimization_results')

# Create the output directory if needed
os.makedirs(OUTPUT_DIR, exist_ok=True)



# ============================================================
# General flags
# ============================================================

META_ONLY = False
SEEDS_EXPERIMENTAIS = [42, 123, 456, 789, 1011]

# ============================================================
# Local input paths
# ============================================================

SOURCE_DATA_DIR = os.path.join(SCRIPT_DIR, "source_data")
CAMINHO_HELIPONTOS_XLS = os.path.join(SOURCE_DATA_DIR, "anac_helipad_registry_snapshot.xls")
GRID_1KM_URL = os.path.join(SOURCE_DATA_DIR, "ibge_statistical_grid_2022_corrected_1km.zip")
MUN_BR_URL = os.path.join(SOURCE_DATA_DIR, "sao_paulo_municipal_boundaries_2022.zip")
ESTADO_BR_URL = MUN_BR_URL

# IBGE municipality codes
COD_SAO_PAULO = "3550308"
COD_CAMPINAS = "3509502"
COD_GUARULHOS = "3518800"
COD_ESTADO_SP = "35"

# ============================================================
# General utility functions
# ============================================================

def haversine_km(lat1, lon1, lat2, lon2):
    """Haversine distance (km) between (lat, lon) points given in degrees."""
    R = 6371.0
    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    dlat = lat2_rad - lat1_rad
    dlon = np.radians(lon2 - lon1)

    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1_rad) * np.cos(lat2_rad) * np.sin(dlon / 2.0) ** 2
    )
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return R * c


# ============================================================
# Candidate-site input
# ============================================================

def carregar_helipontos(caminho_xls: str) -> pd.DataFrame:
    """Read the candidate-site table."""
    print("Loading candidate sites...")

    try:
        df = pd.read_csv(caminho_xls, sep="\t", encoding="latin1", decimal=",")
    except Exception:
        df = pd.read_excel(caminho_xls)

    col_renames = {}
    for c in df.columns:
        c_clean = c.strip()
        if c_clean in ["Código OACI", "Codigo OACI", "CÓDIGO OACI"]:
            col_renames[c] = "Codigo"
        elif c_clean.lower() == "nome":
            col_renames[c] = "Nome"
        elif "município" in c_clean.lower():
            col_renames[c] = "Municipio"
        elif c_clean.lower() == "tipo":
            col_renames[c] = "Tipo"
        elif "latgeopoint" in c_clean.replace(" ", "").lower():
            col_renames[c] = "LatGeoPoint"
        elif "longeopoint" in c_clean.replace(" ", "").lower():
            col_renames[c] = "LonGeoPoint"
        elif c_clean.lower() == "latitude":
            col_renames[c] = "Latitude_txt"
        elif c_clean.lower() == "longitude":
            col_renames[c] = "Longitude_txt"

    df = df.rename(columns=col_renames)

    if "LatGeoPoint" not in df.columns or "LonGeoPoint" not in df.columns:
        if "Latitude_txt" in df.columns and "Longitude_txt" in df.columns:
            df["LatGeoPoint"] = pd.to_numeric(df["Latitude_txt"], errors="coerce")
            df["LonGeoPoint"] = pd.to_numeric(df["Longitude_txt"], errors="coerce")
        else:
            raise ValueError("Could not find latitude/longitude columns in the candidate-site file.")

    if "Tipo" not in df.columns:
        df["Tipo"] = "No solo"

    df["LatGeoPoint"] = pd.to_numeric(df["LatGeoPoint"], errors="coerce")
    df["LonGeoPoint"] = pd.to_numeric(df["LonGeoPoint"], errors="coerce")

    helipontos_df = (
        df[["LatGeoPoint", "LonGeoPoint", "Tipo"]]
        .dropna(subset=["LatGeoPoint", "LonGeoPoint"])
        .reset_index(drop=True)
    )

    helipontos_df = helipontos_df.rename(
        columns={"LatGeoPoint": "Latitude", "LonGeoPoint": "Longitude"}
    )

    print(f"Loaded {len(helipontos_df)} candidate sites")
    return helipontos_df


# ============================================================
# Demand: IBGE 1 km grid (multiple municipalities)
# ============================================================

def _filtrar_municipios(mun_br: gpd.GeoDataFrame, codigos: List[str]) -> gpd.GeoDataFrame:
    """Filter municipalities by a list of IBGE codes."""
    candidatos_codigo = ["CD_MUN", "CD_MUNICIP", "CD_GEOCMU", "CD_GEOCMUN", "GEOCODIGO"]
    
    for col in candidatos_codigo:
        if col in mun_br.columns:
            valores = mun_br[col].astype(str).str.zfill(7)
            sel = mun_br[valores.isin(codigos)]
            if len(sel) > 0:
                print(f"Municipalities identified by column '{col}'")
                return sel
    
    for col in mun_br.columns:
        try:
            valores = mun_br[col].astype(str).str.zfill(7)
            sel = mun_br[valores.isin(codigos)]
            if len(sel) > 0:
                print(f"Municipalities identified by column '{col}' (fallback)")
                return sel
        except Exception:
            continue
    
    raise ValueError(f"Could not find municipalities: {codigos}")


def _filtrar_estado(mun_br: gpd.GeoDataFrame, cod_uf: str) -> gpd.GeoDataFrame:
    """Filter a state by UF code (e.g., '35' for SP)."""
    candidatos_uf = ["CD_UF", "SIGLA_UF", "CD_GEOCUF"]
    
    for col in candidatos_uf:
        if col in mun_br.columns:
            valores = mun_br[col].astype(str)
            sel = mun_br[valores.str.startswith(cod_uf)]
            if len(sel) > 0:
                print(f"State identified by column '{col}'")
                return sel
    
    # Fallback: filter by municipality code starting with the UF code
    for col in mun_br.columns:
        try:
            if "CD_" in col or "GEOCOD" in col:
                valores = mun_br[col].astype(str).str.zfill(7)
                sel = mun_br[valores.str.startswith(cod_uf)]
                if len(sel) > 0:
                    print(f"State identified by column '{col}' (fallback)")
                    return sel
        except Exception:
            continue
    
    raise ValueError(f"Could not filter the state: {cod_uf}")


def carregar_demanda_grade_ibge_1km(codigos_municipios: List[str] = None, estado: str = None) -> Tuple[pd.DataFrame, gpd.GeoDataFrame]:
    """
    Read the 1 km statistical grid and the 2022 municipal boundaries.
    
    Args:
        codigos_municipios: list of IBGE municipality codes (e.g., ["3550308"])
        estado: UF code (e.g., "35" for SP); uses all municipalities in the state
    """
    print("Reading the 1 km national grid (this may take a while)...")
    grid = gpd.read_file(GRID_1KM_URL)
    print(f"Grid loaded. CRS: {grid.crs}, total cells: {len(grid)}")

    col_pop = None
    if "TOTAL" in grid.columns:
        col_pop = "TOTAL"
    else:
        for c in grid.columns:
            if "POP" in c.upper():
                col_pop = c
                break

    if col_pop is None:
        raise ValueError("Could not find a population column in the IBGE grid.")

    print(f"Using population column: {col_pop}")

    print("Reading the 2022 municipal boundaries...")
    mun_br = gpd.read_file(MUN_BR_URL)
    print(f"Boundaries loaded. CRS: {mun_br.crs}, total: {len(mun_br)}")

    if estado:
        mun_selecionados = _filtrar_estado(mun_br, estado)
        print(f"Municipalities in state {estado}: {len(mun_selecionados)}")
    elif codigos_municipios:
        mun_selecionados = _filtrar_municipios(mun_br, codigos_municipios)
        print(f"Selected municipalities: {len(mun_selecionados)}")
    else:
        raise ValueError("Deve fornecer 'codigos_municipios' ou 'estado'")

    geom_uniao = mun_selecionados.unary_union

    print("Clipping the grid to the selected region...")
    grid_regiao = grid[grid.intersects(geom_uniao)].copy()
    grid_regiao.reset_index(drop=True, inplace=True)
    print(f"1 km cells in the region: {len(grid_regiao)}")

    grid_regiao_proj = grid_regiao.to_crs(epsg=31983)
    grid_regiao["Area_km2"] = grid_regiao_proj.geometry.area / 1_000_000.0

    grid_regiao["Populacao"] = pd.to_numeric(grid_regiao[col_pop], errors="coerce").fillna(0.0)

    grid_regiao_latlon = grid_regiao.to_crs(epsg=4326)
    cent = grid_regiao_latlon.geometry.centroid

    demanda_df = pd.DataFrame({
        "Latitude": cent.y.values,
        "Longitude": cent.x.values,
        "Area_km2": grid_regiao["Area_km2"].values,
        "Populacao": grid_regiao["Populacao"].values,
    })

    mun_latlon = mun_selecionados.to_crs(epsg=4326)

    area_total = demanda_df["Area_km2"].sum()
    pop_total = demanda_df["Populacao"].sum()
    print(
        f"Demand (1 km grid): {len(demanda_df)} cells; "
        f"area ~ {area_total:.1f} km^2; population ~ {pop_total/1e6:.2f} million"
    )

    return demanda_df, mun_latlon


# ============================================================
# DENSIDADE POPULACIONAL
# ============================================================

def calcular_densidade_populacional(dem_df: pd.DataFrame, hel_df: pd.DataFrame, raio_km: float = 5.0) -> np.ndarray:
    """Calcula densidade populacional em raio de cada vertiport."""
    J = len(hel_df)
    lat_d = dem_df["Latitude"].values
    lon_d = dem_df["Longitude"].values
    pop_d = dem_df["Populacao"].values
    lat_h = hel_df["Latitude"].values
    lon_h = hel_df["Longitude"].values
    
    densidade_j = np.zeros(J)
    
    for j in range(J):
        dist_ij = haversine_km(lat_d, lon_d, lat_h[j], lon_h[j])
        dentro_raio = dist_ij <= raio_km
        densidade_j[j] = pop_d[dentro_raio].sum()
    
    return densidade_j


# ============================================================
# MATRIZES DE ALCANCE
# ============================================================

def build_reach_matrix_aircraft(dem_df: pd.DataFrame, hel_df: pd.DataFrame, raio_aeronave_km: float = 10.0) -> np.ndarray:
    """Matriz de alcance territorial (aeronave)."""
    I = len(dem_df)
    J = len(hel_df)
    reach = np.zeros((I, J), dtype=bool)

    lat_d = dem_df["Latitude"].values
    lon_d = dem_df["Longitude"].values
    lat_h = hel_df["Latitude"].values
    lon_h = hel_df["Longitude"].values

    for j in range(J):
        dist_ij = haversine_km(lat_d, lon_d, lat_h[j], lon_h[j])
        reach[:, j] = dist_ij <= raio_aeronave_km

    return reach


def build_reach_matrix_facility(dem_df: pd.DataFrame, hel_df: pd.DataFrame, raio_solo_km: float, raio_elev_km: float) -> np.ndarray:
    """Matriz de alcance populacional (facilidade)."""
    I = len(dem_df)
    J = len(hel_df)
    reach = np.zeros((I, J), dtype=bool)

    lat_d = dem_df["Latitude"].values
    lon_d = dem_df["Longitude"].values
    lat_h = hel_df["Latitude"].values
    lon_h = hel_df["Longitude"].values
    tipos = hel_df["Tipo"].values

    for j in range(J):
        if isinstance(tipos[j], str) and ("solo" in tipos[j].lower()):
            raio = raio_solo_km
        else:
            raio = raio_elev_km
        dist_ij = haversine_km(lat_d, lon_d, lat_h[j], lon_h[j])
        reach[:, j] = dist_ij <= raio

    return reach


def compute_dist_vv(hel_df: pd.DataFrame) -> np.ndarray:
    """Pairwise vertiport distances (km): a J x J matrix."""
    J = len(hel_df)
    lat = hel_df["Latitude"].values
    lon = hel_df["Longitude"].values
    dist_vv = np.zeros((J, J), dtype=float)
    for j in range(J):
        dist_vv[j, :] = haversine_km(lat[j], lon[j], lat, lon)
    return dist_vv


def territorial_coverage_denominators(
    area_km2: np.ndarray,
    reach_aircraft: np.ndarray,
) -> Tuple[float, float, np.ndarray]:
    '''Return total-region area, candidate-reachable area, and its mask.'''
    area = np.asarray(area_km2, dtype=float).reshape(-1)
    reach = np.asarray(reach_aircraft, dtype=bool)
    if reach.shape[0] != area.size:
        raise ValueError('reach_aircraft rows must match area_km2')
    coverable = np.any(reach, axis=1)
    area_region = float(area.sum())
    area_candidate_reachable = float(area[coverable].sum())
    if area_region <= 0.0:
        raise ValueError('Total region area must be positive')
    if area_candidate_reachable <= 0.0:
        raise ValueError('No grid area is reachable from the candidate set')
    return area_region, area_candidate_reachable, coverable


def spacing_diagnostics(
    x_sol: np.ndarray,
    dist_vv: np.ndarray,
    min_sep_km: float,
) -> Dict[str, Any]:
    '''Compute the common smooth spacing term and its diagnostics.'''
    selected = np.flatnonzero(np.asarray(x_sol).astype(bool))
    if selected.size < 2:
        return {
            'spacing_pair_count': 0,
            'minimum_selected_distance_km': None,
            'spacing_penalty_raw': 0.0,
        }
    sub = np.asarray(dist_vv, dtype=float)[np.ix_(selected, selected)]
    upper = sub[np.triu_indices(selected.size, k=1)]
    below = upper < float(min_sep_km)
    raw = np.maximum(0.0, 1.0 - upper[below] / float(min_sep_km)).sum()
    return {
        'spacing_pair_count': int(below.sum()),
        'minimum_selected_distance_km': float(upper.min()),
        'spacing_penalty_raw': float(raw),
    }


def verify_allocation_certificate(
    *,
    x_sol: np.ndarray,
    allocation_i: np.ndarray,
    allocation_j: np.ndarray,
    allocation_q: np.ndarray,
    demand: np.ndarray,
    capacities: np.ndarray,
    area_km2: np.ndarray,
    reach_facility: np.ndarray,
    reach_aircraft: np.ndarray,
    dist_vv: np.ndarray,
    eta_demanda: float,
    eta_territorio: float,
    min_sep_km: float,
    pref_score_j: np.ndarray = None,
    w_open: float = 1.0,
    w_pref: float = 0.0,
    w_sep_penalty: float = 0.0,
    tolerance: float = 1e-7,  # certificate tolerance in native constraint units
) -> Dict[str, Any]:
    '''Verify and summarize a sparse continuous q_ij allocation.'''
    x = np.asarray(x_sol, dtype=int).reshape(-1)
    q_i = np.asarray(allocation_i, dtype=np.int64).reshape(-1)
    q_j = np.asarray(allocation_j, dtype=np.int64).reshape(-1)
    q = np.asarray(allocation_q, dtype=float).reshape(-1)
    D = np.asarray(demand, dtype=float).reshape(-1)
    C = np.asarray(capacities, dtype=float).reshape(-1)
    A = np.asarray(area_km2, dtype=float).reshape(-1)
    reach_fac = np.asarray(reach_facility, dtype=bool)
    reach_air = np.asarray(reach_aircraft, dtype=bool)

    if not (q_i.size == q_j.size == q.size):
        raise ValueError('Allocation arrays must have equal lengths')
    if x.size != C.size or reach_fac.shape != (D.size, C.size):
        raise ValueError('Certificate dimensions do not match the instance')

    finite_q = np.isfinite(q)
    valid_index = (q_i >= 0) & (q_i < D.size) & (q_j >= 0) & (q_j < C.size)
    positive = finite_q & valid_index & (q > tolerance)
    by_cell = np.bincount(q_i[positive], weights=q[positive], minlength=D.size)
    by_facility = np.bincount(q_j[positive], weights=q[positive], minlength=C.size)
    total_allocated = float(q[positive].sum())
    demand_total = float(D.sum())

    negative_q_violation = float(max(0.0, -np.min(q[finite_q]))) if finite_q.any() else 0.0
    nonfinite_q_count = int((~finite_q).sum())
    invalid_index_count = int((~valid_index).sum())
    cell_max_violation = float(max(0.0, np.max(by_cell - D))) if D.size else 0.0
    capacity_limit = C * (x > 0)
    capacity_max_violation = float(max(0.0, np.max(by_facility - capacity_limit))) if C.size else 0.0
    access_violation_count = int((~reach_fac[q_i[positive], q_j[positive]]).sum())
    closed_facility_allocation_count = int((x[q_j[positive]] <= 0).sum())

    area_region, area_candidate_reachable, _ = territorial_coverage_denominators(A, reach_air)
    selected = np.flatnonzero(x > 0)
    covered = np.zeros(A.size, dtype=bool)
    for j in selected:
        covered[np.flatnonzero(reach_air[:, j])] = True
    area_covered = float(A[covered].sum())
    cov_dem = total_allocated / (demand_total + 1e-15)
    cov_ter_svc = area_covered / area_candidate_reachable
    cov_ter_region = area_covered / area_region

    spacing = spacing_diagnostics(x, dist_vv, min_sep_km)
    pref = np.zeros_like(C) if pref_score_j is None else np.asarray(pref_score_j, dtype=float).reshape(-1)
    objective = (
        float(w_open) * float((x > 0).sum())
        + float(w_pref) * float(np.dot(pref, x > 0))
        + float(w_sep_penalty) * spacing['spacing_penalty_raw']
    )
    demand_floor_margin = total_allocated - float(eta_demanda) * demand_total
    territorial_floor_margin_km2 = area_covered - float(eta_territorio) * area_candidate_reachable
    feasible = (
        nonfinite_q_count == 0
        and invalid_index_count == 0
        and negative_q_violation <= tolerance
        and cell_max_violation <= tolerance
        and capacity_max_violation <= tolerance
        and access_violation_count == 0
        and closed_facility_allocation_count == 0
        and demand_floor_margin >= -tolerance
        and territorial_floor_margin_km2 >= -tolerance
    )

    digest = hashlib.sha256()
    for arr in (
        x.astype(np.int8), q_i[positive].astype('<i8'),
        q_j[positive].astype('<i8'), q[positive].astype('<f8'),
    ):
        digest.update(np.ascontiguousarray(arr).tobytes())

    summary = {
        'certificate_sha256': digest.hexdigest(),
        'feasibility_tolerance_native_units': float(tolerance),
        'feasible': bool(feasible),
        'facilities': int((x > 0).sum()),
        'allocation_arc_count': int(positive.sum()),
        'allocated_demand': total_allocated,
        'total_demand': demand_total,
        'demand_coverage': float(cov_dem),
        'demand_floor_margin': float(demand_floor_margin),
        'demand_floor_residual': float(cov_dem - float(eta_demanda)),
        'cell_allocation_max_violation': cell_max_violation,
        'capacity_max_violation': capacity_max_violation,
        'access_violation_count': access_violation_count,
        'closed_facility_allocation_count': closed_facility_allocation_count,
        'negative_q_max_violation': negative_q_violation,
        'nonfinite_q_count': nonfinite_q_count,
        'invalid_index_count': invalid_index_count,
        'area_covered_km2': area_covered,
        'area_candidate_reachable_km2': area_candidate_reachable,
        'area_region_km2': area_region,
        'candidate_reachable_territorial_coverage': float(cov_ter_svc),
        'total_region_territorial_coverage': float(cov_ter_region),
        'territorial_floor_margin_km2': float(territorial_floor_margin_km2),
        'territorial_floor_residual': float(cov_ter_svc - float(eta_territorio)),
        'spacing_pair_count': spacing['spacing_pair_count'],
        'minimum_selected_distance_km': spacing['minimum_selected_distance_km'],
        'spacing_penalty_raw': spacing['spacing_penalty_raw'],
        'spacing_penalty_contribution': float(w_sep_penalty) * spacing['spacing_penalty_raw'],
        'objective': float(objective),
    }
    return {
        'x': x.astype(np.int8),
        'allocation_i': q_i[positive],
        'allocation_j': q_j[positive],
        'allocation_q': q[positive],
        'summary': summary,
    }


def save_feasibility_certificate(
    certificate: Dict[str, Any],
    output_prefix: str,
) -> Tuple[str, str]:
    '''Save raw certificate arrays and a human-readable summary.'''
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    npz_path = str(Path(str(prefix) + '.npz'))
    json_path = str(Path(str(prefix) + '.json'))
    np.savez_compressed(
        npz_path,
        x=np.asarray(certificate['x'], dtype=np.int8),
        allocation_i=np.asarray(certificate['allocation_i'], dtype=np.int64),
        allocation_j=np.asarray(certificate['allocation_j'], dtype=np.int64),
        allocation_q=np.asarray(certificate['allocation_q'], dtype=float),
    )
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(certificate['summary'], f, indent=2, ensure_ascii=False)
    return npz_path, json_path


def update_final_results_manifest(
    label: str,
    paths: List[str],
    status: str,
    metadata: Dict[str, Any] = None,
    manifest_path: str = None,
) -> str:
    '''Atomically update the machine-readable final-results manifest.'''
    if manifest_path is None:
        manifest_path = os.environ.get(
            'IEEE_ACCESS_MANIFEST_PATH',
            os.path.join(SCRIPT_DIR, 'final_outputs', 'final_results_manifest.json'),
        )
    manifest_file = Path(manifest_path)
    manifest_file.parent.mkdir(parents=True, exist_ok=True)
    if manifest_file.exists():
        with open(manifest_file, 'r', encoding='utf-8') as f:
            manifest = json.load(f)
    else:
        manifest = {
            'schema_version': '1.0',
            'status': 'in_progress',
            'outputs': {},
        }

    file_entries = []
    for raw_path in paths:
        if raw_path is None:
            continue
        path = Path(raw_path).resolve()
        entry = {
            'path': os.path.relpath(str(path), SCRIPT_DIR).replace('\\', '/'),
            'exists': path.is_file(),
        }
        if path.is_file():
            digest = hashlib.sha256()
            with open(path, 'rb') as f:
                for chunk in iter(lambda: f.read(1024 * 1024), b''):
                    digest.update(chunk)
            entry['sha256'] = digest.hexdigest()
            entry['size_bytes'] = path.stat().st_size
        file_entries.append(entry)

    now = datetime.now(timezone.utc).isoformat()
    manifest['outputs'][label] = {
        'status': status,
        'updated_utc': now,
        'files': file_entries,
        'metadata': metadata or {},
    }
    manifest['updated_utc'] = now
    required = {'light_medium_meta', 'light_medium_exact_q', 'heavy_component'}
    manifest['status'] = (
        'complete'
        if required.issubset(manifest['outputs'])
        and all(manifest['outputs'][k]['status'] == 'complete' for k in required)
        else 'in_progress'
    )
    temp_path = manifest_file.with_suffix('.json.tmp')
    with open(temp_path, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    os.replace(temp_path, manifest_file)
    return str(manifest_file)


# ============================================================
# Location-preference score
# ============================================================

def compute_preferencia_melhorada(
    hel_df: pd.DataFrame,
    dem_df: pd.DataFrame,
    densidade_j: np.ndarray,
    dist_vv: np.ndarray,
    w_densidade: float = 1.0,
    w_tipo: float = 0.5,
    w_overlap: float = 0.3,
    raio_prox_km: float = 3.0,
) -> np.ndarray:
    """Preference score: lower is more preferable."""
    J = len(hel_df)
    
    if densidade_j.max() > 0:
        dens_norm = densidade_j / densidade_j.max()
    else:
        dens_norm = np.zeros(J)
    pen_densidade = w_densidade * (1.0 - dens_norm)
    
    tipos = hel_df["Tipo"].values
    is_solo = np.array([("solo" in str(t).lower()) for t in tipos], dtype=float)
    pen_tipo = w_tipo * (1.0 - is_solo)
    
    overlap = np.zeros(J)
    for j in range(J):
        overlap[j] = (dist_vv[j, :] < raio_prox_km).sum() - 1
    if overlap.max() > 0:
        overlap_norm = overlap / overlap.max()
    else:
        overlap_norm = np.zeros(J)
    pen_overlap = w_overlap * overlap_norm
    
    score = pen_densidade + pen_tipo + pen_overlap
    
    return score


# ============================================================
# Exact model: capacitated facility location
# ============================================================

class VertiportCFLP:
    """Capacitated location-allocation model."""

    def __init__(
        self,
        hel_df: pd.DataFrame,
        dem_df: pd.DataFrame,
        reach_aircraft: np.ndarray,
        reach_facility: np.ndarray,
        D_frac: float,
        eta_demanda: float,
        eta_territorio: float,
        w_open: float,
        pref_score_j: np.ndarray,
        w_pref: float,
        dist_vv: np.ndarray,
        min_sep_km: float = 0.3,
        w_sep_penalty: float = 0.0,
        binary_assign: bool = False,
    ):
        self.hel = hel_df.reset_index(drop=True)
        self.dem = dem_df.reset_index(drop=True)
        self.reach_air = reach_aircraft.astype(bool)
        self.reach_fac = reach_facility.astype(bool)
        self.J = len(self.hel)
        self.I = len(self.dem)
        self.D_frac = float(D_frac)
        self.eta_dem = float(eta_demanda)
        self.eta_ter = float(eta_territorio)
        self.w_open = float(w_open)
        self.pref_score = pref_score_j.reshape(-1)
        self.w_pref = float(w_pref)
        self.dist_vv = dist_vv
        self.min_sep = float(min_sep_km)
        self.w_sep_penalty = float(w_sep_penalty)
        # Backwards-compatible keyword only.  The authoritative model always
        # uses continuous allocated-demand q_ij variables.
        if binary_assign:
            warnings.warn(
                "binary_assign is deprecated and ignored; q_ij is continuous",
                RuntimeWarning,
            )
        self.binary_assign = False

        if "Populacao" in self.dem.columns:
            pop = self.dem["Populacao"].values.astype(float)
        else:
            pop = np.ones(self.I, dtype=float)
        self.D = self.D_frac * pop
        self.D_tot = self.D.sum()

        if "Area_km2" in self.dem.columns:
            self.A = self.dem["Area_km2"].values.astype(float)
        else:
            self.A = np.ones(self.I, dtype=float)
        self.A_tot = self.A.sum()  # legacy alias: total-region area
        (
            self.A_tot_region,
            self.A_max_coverable,
            self.candidate_coverable_mask,
        ) = territorial_coverage_denominators(self.A, self.reach_air)

        tipos = self.hel["Tipo"].values
        cap_base = np.where(
            np.array([("solo" in str(t).lower()) for t in tipos]),
            1500.0,
            800.0,
        )
        D_local_j = (self.reach_fac * self.D.reshape(-1, 1)).sum(axis=0)
        if np.any(D_local_j > 0):
            media_local = np.mean(D_local_j[D_local_j > 0])
        else:
            media_local = 1.0
        fator = 1.0 + 0.8 * (D_local_j / (media_local if media_local > 0 else 1.0))
        fator = np.clip(fator, 1.0, 3.0)
        self.C = cap_base * fator

        self.pairs_fac = np.argwhere(self.reach_fac)
        self.num_pairs_fac = len(self.pairs_fac)

        self.model = None
        self.x = None
        self.s_ter = None
        self.q = None
        self.y = None  # deprecated alias set to q after build
        self.z_sep = None

        close = np.argwhere(np.triu(self.dist_vv < self.min_sep, k=1))
        self.close_pairs = close.astype(int)
        if len(close):
            self.close_pair_coeff = (
                1.0 - self.dist_vv[close[:, 0], close[:, 1]] / self.min_sep
            )
        else:
            self.close_pair_coeff = np.zeros(0, dtype=float)

    def build(self, verbose: bool = False):
        if self.model is not None:
            return
        print("[CFLP] Building the exact model...")
        m = Model("VertiportCFLP_v3")
        m.Params.OutputFlag = 1 if verbose else 0

        self.model = m
        self.x = m.addVars(self.J, vtype=GRB.BINARY, name="x")
        self.s_ter = m.addVars(self.I, vtype=GRB.BINARY, name="s_ter")
        self.q = m.addVars(
            self.num_pairs_fac, lb=0.0, vtype=GRB.CONTINUOUS, name="q"
        )
        self.y = self.q

        x = self.x
        s_ter = self.s_ter
        q = self.q

        obj_base = (
            self.w_open * quicksum(x[j] for j in range(self.J))
            + self.w_pref * quicksum(self.pref_score[j] * x[j] for j in range(self.J))
        )
        
        if self.w_sep_penalty > 0 and len(self.close_pairs):
            # Linearization of z_jk=x_j*x_k.  Positive coefficients and
            # minimization make the lower linking inequality sufficient.
            self.z_sep = m.addVars(
                len(self.close_pairs), lb=0.0, ub=1.0,
                vtype=GRB.CONTINUOUS, name="z_sep",
            )
            for p, (j, k) in enumerate(self.close_pairs):
                m.addConstr(
                    self.z_sep[p] >= x[int(j)] + x[int(k)] - 1,
                    name=f"link_sep_{p}",
                )
            obj_penalty = self.w_sep_penalty * quicksum(
                float(self.close_pair_coeff[p]) * self.z_sep[p]
                for p in range(len(self.close_pairs))
            )
            obj = obj_base + obj_penalty
        else:
            obj = obj_base
        
        m.setObjective(obj, GRB.MINIMIZE)

        from collections import defaultdict
        idx_by_i_fac = defaultdict(list)
        idx_by_j_fac = defaultdict(list)
        for k, (ii, jj) in enumerate(self.pairs_fac):
            idx_by_i_fac[int(ii)].append(k)
            idx_by_j_fac[int(jj)].append(k)

        for i in range(self.I):
            ks = idx_by_i_fac.get(i, [])
            if ks:
                m.addConstr(
                    quicksum(q[k] for k in ks) <= self.D[i],
                    name=f"alloc_cell_{i}",
                )

        for j in range(self.J):
            ks = idx_by_j_fac.get(j, [])
            if ks:
                m.addConstr(
                    quicksum(q[k] for k in ks) <= self.C[j] * x[j],
                    name=f"cap_{j}",
                )

        m.addConstr(
            quicksum(q[k] for k in range(self.num_pairs_fac))
            >= self.eta_dem * self.D_tot,
            name="cobertura_min_demanda",
        )

        for i in range(self.I):
            j_possiveis = np.where(self.reach_air[i, :])[0]
            if len(j_possiveis) > 0:
                m.addConstr(s_ter[i] <= quicksum(x[j] for j in j_possiveis), name=f"cov_ter_{i}")
            else:
                m.addConstr(s_ter[i] == 0, name=f"no_reach_ter_{i}")

        m.addConstr(
            quicksum(self.A[i] * s_ter[i] for i in range(self.I))
            >= self.eta_ter * self.A_max_coverable,
            name="cobertura_min_territorio",
        )

        print("[CFLP] Model built")

    def solve(self, time_limit: float = 3600.0, verbose: bool = False) -> Dict[str, Any]:
        if self.model is None:
            self.build(verbose=verbose)
        m = self.model
        m.Params.TimeLimit = time_limit

        start = time.time()
        m.optimize()
        tempo = time.time() - start
        status = m.status

        res = {
            "status": status,
            "tempo": tempo,
            "num_vertiports": None,
            "cov_dem_rel": 0.0,
            "cov_ter_rel": 0.0,
            "cov_ter_svc": 0.0,
            "cov_ter_region": 0.0,
            "D_cov": 0.0,
            "A_cov": 0.0,
            "A_max_coverable": self.A_max_coverable,
            "A_tot_region": self.A_tot_region,
            "x": None,
            "s_dem": None,
            "s_ter": None,
            "allocation_i": None,
            "allocation_j": None,
            "allocation_q": None,
            "certificate_summary": None,
            "obj": None,
            "gap": None,
        }

        if status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL] and m.SolCount > 0:
            x_sol = np.array([self.x[j].X for j in range(self.J)])
            s_ter_sol = np.array([self.s_ter[i].X for i in range(self.I)])
            x_sol = (x_sol > 0.5).astype(int)
            s_ter_sol = (s_ter_sol > 0.5).astype(int)

            num_v = x_sol.sum()
            q_all = np.array(
                [self.q[k].X for k in range(self.num_pairs_fac)], dtype=float
            )
            used = q_all > 1e-9
            allocation_pairs = self.pairs_fac[used]
            if used.any():
                allocation_i = allocation_pairs[:, 0].astype(np.int64)
                allocation_j = allocation_pairs[:, 1].astype(np.int64)
            else:
                allocation_i = np.zeros(0, dtype=np.int64)
                allocation_j = np.zeros(0, dtype=np.int64)
            allocation_q = q_all[used]
            D_cov = float(allocation_q.sum())
            selected = np.flatnonzero(x_sol)
            covered = np.zeros(self.I, dtype=bool)
            for j in selected:
                covered[np.flatnonzero(self.reach_air[:, j])] = True
            A_cov = float(self.A[covered].sum())
            cov_dem_rel = D_cov / (self.D_tot + 1e-9)
            cov_ter_svc = A_cov / (self.A_max_coverable + 1e-9)
            cov_ter_region = A_cov / (self.A_tot_region + 1e-9)

            certificate = verify_allocation_certificate(
                x_sol=x_sol,
                allocation_i=allocation_i,
                allocation_j=allocation_j,
                allocation_q=allocation_q,
                demand=self.D,
                capacities=self.C,
                area_km2=self.A,
                reach_facility=self.reach_fac,
                reach_aircraft=self.reach_air,
                dist_vv=self.dist_vv,
                eta_demanda=self.eta_dem,
                eta_territorio=self.eta_ter,
                min_sep_km=self.min_sep,
                pref_score_j=self.pref_score,
                w_open=self.w_open,
                w_pref=self.w_pref,
                w_sep_penalty=self.w_sep_penalty,
            )

            if m.SolCount > 0 and hasattr(m, "MIPGap"):
                try:
                    res["gap"] = float(m.MIPGap)
                except Exception:
                    res["gap"] = None

            res.update({
                "num_vertiports": int(num_v),
                "cov_dem_rel": cov_dem_rel,
                "cov_ter_rel": cov_ter_svc,
                "cov_ter_svc": cov_ter_svc,
                "cov_ter_region": cov_ter_region,
                "D_cov": D_cov,
                "A_cov": A_cov,
                "x": x_sol,
                "s_ter": s_ter_sol,
                "allocation_i": allocation_i,
                "allocation_j": allocation_j,
                "allocation_q": allocation_q,
                "certificate_summary": certificate["summary"],
                "obj": float(m.ObjVal),
            })
            print(
                f"[CFLP] Status={status}, Verts={num_v}, "
                f"cov_dem={cov_dem_rel*100:.2f}%, "
                f"cov_ter_svc={cov_ter_svc*100:.2f}%, "
                f"cov_ter_region={cov_ter_region*100:.2f}%, "
                f"obj={res['obj']:.3f}, tempo={tempo:.1f}s"
            )
        else:
            print(f"[CFLP] Status={status} (no feasible solution).")

        return res


# ============================================================
# Base class for the stochastic methods
# ============================================================

class MetaheuristicaBase:
    """Base class for the stochastic search methods."""
    
    def __init__(
        self,
        helipontos_df: pd.DataFrame,
        demanda_df: pd.DataFrame,
        reach_aircraft: np.ndarray,
        reach_facility: np.ndarray,
        dist_vv: np.ndarray,
        D_frac: float,
        eta_demanda: float,
        eta_territorio: float,
        w_verts: float,
        pref_score_j: np.ndarray,
        w_pref: float,
        lambda_cov_dem: float,
        lambda_cov_ter: float,
        lambda_sep: float,
        min_sep_km: float = 0.3,
        seed: int = 42,
        fitness_call_budget: int = None,
    ):
        np.random.seed(seed)
        self.hel = helipontos_df.reset_index(drop=True)
        self.dem = demanda_df.reset_index(drop=True)
        self.reach_air = reach_aircraft.astype(bool)
        self.reach_fac = reach_facility.astype(bool)
        self.dist_vv = dist_vv
        self.J = len(self.hel)
        self.I = len(self.dem)
        self.seed = seed

        if "Populacao" in self.dem.columns:
            pop = self.dem["Populacao"].values.astype(float)
        else:
            pop = np.ones(self.I, dtype=float)
        self.D = D_frac * pop
        self.D_tot = self.D.sum()

        if "Area_km2" in self.dem.columns:
            self.A = self.dem["Area_km2"].values.astype(float)
        else:
            self.A = np.ones(self.I, dtype=float)
        self.A_tot = self.A.sum()  # legacy alias: total-region area
        (
            self.A_tot_region,
            self.A_max_coverable,
            self.candidate_coverable_mask,
        ) = territorial_coverage_denominators(self.A, self.reach_air)

        tipos = self.hel["Tipo"].values
        cap_base = np.where(
            np.array([("solo" in str(t).lower()) for t in tipos]),
            1500.0,
            800.0,
        )
        D_local_j = (self.reach_fac * self.D.reshape(-1, 1)).sum(axis=0)
        if np.any(D_local_j > 0):
            media_local = np.mean(D_local_j[D_local_j > 0])
        else:
            media_local = 1.0
        fator = 1.0 + 0.8 * (D_local_j / (media_local if media_local > 0 else 1.0))
        fator = np.clip(fator, 1.0, 3.0)
        self.C = cap_base * fator

        self.eta_dem = float(eta_demanda)
        self.eta_ter = float(eta_territorio)
        self.w_verts = float(w_verts)
        self.pref_score = pref_score_j.reshape(-1)
        self.w_pref = float(w_pref)
        self.lambda_cov_dem = float(lambda_cov_dem)
        self.lambda_cov_ter = float(lambda_cov_ter)
        self.lambda_sep = float(lambda_sep)
        self.min_sep = min_sep_km
        self.fitness_calls = 0
        self.fitness_call_budget = (
            None if fitness_call_budget is None else int(fitness_call_budget)
        )
        if self.fitness_call_budget is not None and self.fitness_call_budget <= 0:
            raise ValueError("fitness_call_budget must be positive")

        # These lists eliminate repeated I-row scans and sorts in every
        # fitness call.  The change is output-equivalent and is essential for
        # the 250,827-cell Heavy instance.
        self.reach_fac_sorted_indices = []
        self.reach_air_indices = []
        for j in range(self.J):
            idx_fac = np.flatnonzero(self.reach_fac[:, j])
            if idx_fac.size:
                idx_fac = idx_fac[np.argsort(-self.D[idx_fac])]
            self.reach_fac_sorted_indices.append(idx_fac)
            self.reach_air_indices.append(np.flatnonzero(self.reach_air[:, j]))

        close = np.argwhere(np.triu(self.dist_vv < self.min_sep, k=1))
        self.close_pairs = close.astype(int)
        if len(close):
            self.close_pair_coeff = (
                1.0 - self.dist_vv[close[:, 0], close[:, 1]] / self.min_sep
            )
        else:
            self.close_pair_coeff = np.zeros(0, dtype=float)

    def _can_evaluate(self, count: int = 1) -> bool:
        return (
            self.fitness_call_budget is None
            or self.fitness_calls + int(count) <= self.fitness_call_budget
        )

    def _assign_demand_greedy(self, sol: np.ndarray, return_allocation: bool = False):
        sol = sol.astype(bool)
        ativos = np.where(sol)[0]
        if len(ativos) == 0:
            empty = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64),
                     np.zeros(0, dtype=float))
            if return_allocation:
                return 0.0, np.zeros(self.I, dtype=float), empty
            return 0.0, np.zeros(self.I, dtype=float)

        D_rem = self.D.copy()
        cap_rem = self.C[ativos].copy()
        servida = np.zeros(self.I, dtype=float)
        alloc_i, alloc_j, alloc_q = [], [], []

        for idx_j, j in enumerate(ativos):
            idx_i = self.reach_fac_sorted_indices[j]
            if len(idx_i) == 0:
                continue
            for i in idx_i:
                if cap_rem[idx_j] <= 0:
                    break
                if D_rem[i] <= 0:
                    continue
                aloca = min(D_rem[i], cap_rem[idx_j])
                servida[i] += aloca
                D_rem[i] -= aloca
                cap_rem[idx_j] -= aloca
                if return_allocation and aloca > 0:
                    alloc_i.append(int(i))
                    alloc_j.append(int(j))
                    alloc_q.append(float(aloca))

        D_cov = servida.sum()
        if return_allocation:
            allocation = (
                np.asarray(alloc_i, dtype=np.int64),
                np.asarray(alloc_j, dtype=np.int64),
                np.asarray(alloc_q, dtype=float),
            )
            return D_cov, servida, allocation
        return D_cov, servida

    def _compute_territorial_coverage(self, sol: np.ndarray):
        sol = sol.astype(bool)
        ativos = np.where(sol)[0]
        if len(ativos) == 0:
            return 0.0

        coberta = np.zeros(self.I, dtype=bool)
        for j in ativos:
            coberta[self.reach_air_indices[j]] = True
        A_cov = self.A[coberta].sum()
        return float(A_cov)

    def _penalizacao_proximidade_suave(self, sol: np.ndarray) -> float:
        sol = sol.astype(bool)
        if not len(self.close_pairs):
            return 0.0
        active_pairs = (
            sol[self.close_pairs[:, 0]] & sol[self.close_pairs[:, 1]]
        )
        return float(self.close_pair_coeff[active_pairs].sum())

    def fitness(self, sol: np.ndarray) -> Tuple[float, Dict[str, float]]:
        if not self._can_evaluate():
            raise RuntimeError("fitness-call budget exhausted")
        self.fitness_calls += 1
        sol = sol.astype(int)
        num_v = sol.sum()

        D_cov, _ = self._assign_demand_greedy(sol)
        cov_dem_rel = D_cov / (self.D_tot + 1e-9)
        shortfall_dem = max(0.0, self.eta_dem - cov_dem_rel)
        pen_cov_dem = self.lambda_cov_dem * shortfall_dem if shortfall_dem > 1e-6 else 0.0

        A_cov = self._compute_territorial_coverage(sol)
        cov_ter_svc = A_cov / (self.A_max_coverable + 1e-9)
        cov_ter_region = A_cov / (self.A_tot_region + 1e-9)
        shortfall_ter = max(0.0, self.eta_ter - cov_ter_svc)
        pen_cov_ter = self.lambda_cov_ter * shortfall_ter if shortfall_ter > 1e-6 else 0.0

        pen_prox = self._penalizacao_proximidade_suave(sol)
        pen_sep = self.lambda_sep * pen_prox
        spacing = spacing_diagnostics(sol, self.dist_vv, self.min_sep)

        pref_term = self.w_pref * float((self.pref_score * sol).sum())

        fit = self.w_verts * num_v + pref_term + pen_cov_dem + pen_cov_ter + pen_sep

        info = {
            "num_v": float(num_v),
            "cov_dem_rel": float(cov_dem_rel),
            "cov_ter_rel": float(cov_ter_svc),
            "cov_ter_svc": float(cov_ter_svc),
            "cov_ter_region": float(cov_ter_region),
            "D_cov": float(D_cov),
            "A_cov": float(A_cov),
            "A_max_coverable": float(self.A_max_coverable),
            "A_tot_region": float(self.A_tot_region),
            "shortfall_dem": float(shortfall_dem),
            "shortfall_ter": float(shortfall_ter),
            "pen_cov_dem": float(pen_cov_dem),
            "pen_cov_ter": float(pen_cov_ter),
            "pen_prox": float(pen_prox),
            "pen_sep": float(pen_sep),
            "spacing_pair_count": spacing["spacing_pair_count"],
            "minimum_selected_distance_km": spacing["minimum_selected_distance_km"],
            "fitness_call_at_evaluation": int(self.fitness_calls),
        }
        return fit, info

    def build_feasibility_certificate(self, sol: np.ndarray) -> Dict[str, Any]:
        '''Construct and verify the deterministic sparse q_ij allocation.'''
        _, _, allocation = self._assign_demand_greedy(
            np.asarray(sol, dtype=int), return_allocation=True
        )
        allocation_i, allocation_j, allocation_q = allocation
        return verify_allocation_certificate(
            x_sol=sol,
            allocation_i=allocation_i,
            allocation_j=allocation_j,
            allocation_q=allocation_q,
            demand=self.D,
            capacities=self.C,
            area_km2=self.A,
            reach_facility=self.reach_fac,
            reach_aircraft=self.reach_air,
            dist_vv=self.dist_vv,
            eta_demanda=self.eta_dem,
            eta_territorio=self.eta_ter,
            min_sep_km=self.min_sep,
            pref_score_j=self.pref_score,
            w_open=self.w_verts,
            w_pref=self.w_pref,
            w_sep_penalty=self.lambda_sep,
        )


# ============================================================
# SMA-INSPIRED BINARY SEARCH (legacy serialized ID: SMA)
# ============================================================

class MetaheuristicaSMA(MetaheuristicaBase):
    """Elite-based SMA-inspired binary search (not canonical SMA)."""
    
    def optimize(
        self,
        n_pop: int = 30,
        max_iters: int = 100,
        p_init_open: float = 0.3,
        patience: int = 50,
        verbose: bool = True,
    ):
        J = self.J
        rng = np.random.RandomState(self.seed)
        self.fitness_calls = 0
        stop_reason = "maximum_iterations"
        if not self._can_evaluate(n_pop):
            raise ValueError("fitness_call_budget must allow the initial SMA-inspired binary-search population")

        pop = (rng.rand(n_pop, J) < p_init_open).astype(int)
        fitness_vals = np.zeros(n_pop)
        infos = [None] * n_pop

        best_fit = float("inf")
        best_sol = None
        best_info = None
        last_improvement = 0

        history = []
        t0 = time.time()

        def avalia_pop(iter_num: int):
            nonlocal best_fit, best_sol, best_info, last_improvement
            for idx in range(n_pop):
                f, inf = self.fitness(pop[idx])
                fitness_vals[idx] = f
                infos[idx] = inf
                if f < best_fit:
                    best_fit = f
                    best_sol = pop[idx].copy()
                    best_info = inf
                    last_improvement = iter_num
            history.append({
                "iter": iter_num,
                "tempo": time.time() - t0,
                "best_fit": best_fit,
                "num_v": best_info["num_v"],
                "cov_dem_rel": best_info["cov_dem_rel"],
                "cov_ter_rel": best_info["cov_ter_rel"],
                "cov_ter_svc": best_info["cov_ter_svc"],
                "cov_ter_region": best_info["cov_ter_region"],
                "fitness_calls": int(self.fitness_calls),
            })

        avalia_pop(iter_num=0)
        if verbose:
            print(
                f"[SIBS] initial: fit={best_fit:.4f}, "
                f"verts={best_info['num_v']:.0f}, "
                f"cov_dem={best_info['cov_dem_rel']*100:.2f}%, "
                f"cov_ter={best_info['cov_ter_rel']*100:.2f}%"
            )

        for it in range(1, max_iters + 1):
            if not self._can_evaluate(n_pop):
                stop_reason = "fitness_call_budget"
                break
            ordem = np.argsort(fitness_vals)
            pop = pop[ordem]
            fitness_vals = fitness_vals[ordem]
            infos = [infos[i] for i in ordem]

            elite = pop[0]
            second = pop[1] if n_pop > 1 else elite

            for idx in range(n_pop):
                if idx == 0:
                    continue
                base = elite if rng.rand() < 0.7 else second
                candidato = base.copy()
                n_flip = max(1, int(0.1 * J))
                flip_idx = rng.choice(J, size=n_flip, replace=False)
                candidato[flip_idx] = 1 - candidato[flip_idx]
                pop[idx] = candidato

            avalia_pop(iter_num=it)

            if verbose and it % 10 == 0:
                print(
                    f"[SIBS] it={it}, fit={best_fit:.4f}, "
                    f"verts={best_info['num_v']:.0f}, "
                    f"cov_dem={best_info['cov_dem_rel']*100:.2f}%, "
                    f"cov_ter={best_info['cov_ter_rel']*100:.2f}%"
                )

            if it - last_improvement >= patience:
                stop_reason = "sma_patience"
                if verbose:
                    print(f"[SIBS] stopping: no improvement for {patience} iterations.")
                break

        hist_df = pd.DataFrame(history)
        best_info = dict(best_info)
        best_info["fitness_calls"] = int(self.fitness_calls)
        best_info["stopping_reason"] = stop_reason
        return best_sol, best_fit, best_info, hist_df


# ============================================================
# VNS PURO
# ============================================================

class MetaheuristicaVNS(MetaheuristicaBase):
    """Variable Neighborhood Search puro."""
    
    def _shake(self, sol: np.ndarray, k: int) -> np.ndarray:
        sol_new = sol.copy()
        J = len(sol)
        k = min(k, J)
        idx = np.random.choice(J, size=k, replace=False)
        sol_new[idx] = 1 - sol_new[idx]
        return sol_new

    def _local_search(self, sol: np.ndarray, max_moves: int = 50) -> Tuple[np.ndarray, float, Dict[str, float]]:
        current = sol.copy()
        if not self._can_evaluate():
            raise RuntimeError("fitness-call budget exhausted before local search")
        best_fit, best_info = self.fitness(current)
        J = len(sol)
        for _ in range(max_moves):
            if not self._can_evaluate():
                break
            j = np.random.randint(0, J)
            neighbor = current.copy()
            neighbor[j] = 1 - neighbor[j]
            f, inf = self.fitness(neighbor)
            if f < best_fit:
                current = neighbor
                best_fit = f
                best_info = inf
        return current, best_fit, best_info

    def optimize(
        self,
        max_iters: int = 50,
        k_max: int = 5,
        p_init_open: float = 0.3,
        verbose: bool = True,
    ):
        J = self.J
        rng = np.random.RandomState(self.seed)
        self.fitness_calls = 0
        stop_reason = "maximum_iterations"
        if not self._can_evaluate():
            raise ValueError("fitness_call_budget must allow VNS initialization")

        initial_sol = (rng.rand(J) < p_init_open).astype(int)
        best_sol, best_fit, best_info = self._local_search(initial_sol)

        history = []
        t0 = time.time()

        history.append({
            "iter": 0,
            "tempo": time.time() - t0,
            "best_fit": best_fit,
            "num_v": best_info["num_v"],
            "cov_dem_rel": best_info["cov_dem_rel"],
            "cov_ter_rel": best_info["cov_ter_rel"],
            "cov_ter_svc": best_info["cov_ter_svc"],
            "cov_ter_region": best_info["cov_ter_region"],
            "fitness_calls": int(self.fitness_calls),
        })

        if verbose:
            print(
                f"[VNS] initial: fit={best_fit:.4f}, "
                f"verts={best_info['num_v']:.0f}, "
                f"cov_dem={best_info['cov_dem_rel']*100:.2f}%, "
                f"cov_ter={best_info['cov_ter_rel']*100:.2f}%"
            )

        budget_exhausted = False
        for it in range(max_iters):
            for k in range(1, k_max + 1):
                if not self._can_evaluate():
                    budget_exhausted = True
                    break
                shaken = self._shake(best_sol, k)
                local_sol, local_fit, local_info = self._local_search(shaken)
                if local_fit < best_fit:
                    best_sol = local_sol.copy()
                    best_fit = local_fit
                    best_info = local_info
                    if verbose:
                        print(
                            f"[VNS] it={it}, k={k}, new best: "
                            f"fit={best_fit:.4f}, verts={best_info['num_v']:.0f}, "
                            f"cov_dem={best_info['cov_dem_rel']*100:.2f}%, "
                            f"cov_ter={best_info['cov_ter_rel']*100:.2f}%"
                        )
                    break

            history.append({
                "iter": it + 1,
                "tempo": time.time() - t0,
                "best_fit": best_fit,
                "num_v": best_info["num_v"],
                "cov_dem_rel": best_info["cov_dem_rel"],
                "cov_ter_rel": best_info["cov_ter_rel"],
                "cov_ter_svc": best_info["cov_ter_svc"],
                "cov_ter_region": best_info["cov_ter_region"],
                "fitness_calls": int(self.fitness_calls),
            })
            if budget_exhausted:
                stop_reason = "fitness_call_budget"
                break

        hist_df = pd.DataFrame(history)
        best_info = dict(best_info)
        best_info["fitness_calls"] = int(self.fitness_calls)
        best_info["stopping_reason"] = stop_reason
        return best_sol, best_fit, best_info, hist_df


# ============================================================
# SMA-INSPIRED BINARY SEARCH + VNS RELAY HYBRID
# ============================================================

class MetaheuristicaSMA_VNS(MetaheuristicaBase):
    """Relay hybrid: SMA-inspired binary search followed by VNS."""
    
    def _shake(self, sol: np.ndarray, k: int) -> np.ndarray:
        sol_new = sol.copy()
        J = len(sol)
        k = min(k, J)
        idx = np.random.choice(J, size=k, replace=False)
        sol_new[idx] = 1 - sol_new[idx]
        return sol_new

    def _local_search(self, sol: np.ndarray, max_moves: int = 50) -> Tuple[np.ndarray, float, Dict[str, float]]:
        current = sol.copy()
        if not self._can_evaluate():
            raise RuntimeError("fitness-call budget exhausted before local search")
        best_fit, best_info = self.fitness(current)
        J = len(sol)
        for _ in range(max_moves):
            if not self._can_evaluate():
                break
            j = np.random.randint(0, J)
            neighbor = current.copy()
            neighbor[j] = 1 - neighbor[j]
            f, inf = self.fitness(neighbor)
            if f < best_fit:
                current = neighbor
                best_fit = f
                best_info = inf
        return current, best_fit, best_info

    def optimize(
        self,
        n_pop: int = 30,
        max_iters_sma: int = 100,
        max_iters_vns: int = 50,
        k_max: int = 5,
        p_init_open: float = 0.3,
        patience_sma: int = 50,
        verbose: bool = True,
    ):
        J = self.J
        rng = np.random.RandomState(self.seed)
        self.fitness_calls = 0
        stop_reason = "maximum_iterations"
        if not self._can_evaluate(n_pop):
            raise ValueError("fitness_call_budget must allow the initial SMA-inspired binary-search population")

        pop = (rng.rand(n_pop, J) < p_init_open).astype(int)
        fitness_vals = np.zeros(n_pop)
        infos = [None] * n_pop

        best_fit = float("inf")
        best_sol = None
        best_info = None
        last_improvement_sma = 0

        history = []
        t0 = time.time()

        def avalia_pop(iter_sma: int):
            nonlocal best_fit, best_sol, best_info, last_improvement_sma
            for idx in range(n_pop):
                f, inf = self.fitness(pop[idx])
                fitness_vals[idx] = f
                infos[idx] = inf
                if f < best_fit:
                    best_fit = f
                    best_sol = pop[idx].copy()
                    best_info = inf
                    last_improvement_sma = iter_sma
            history.append({
                "fase": "SMA",
                "iter": iter_sma,
                "tempo": time.time() - t0,
                "best_fit": best_fit,
                "num_v": best_info["num_v"],
                "cov_dem_rel": best_info["cov_dem_rel"],
                "cov_ter_rel": best_info["cov_ter_rel"],
                "cov_ter_svc": best_info["cov_ter_svc"],
                "cov_ter_region": best_info["cov_ter_region"],
                "fitness_calls": int(self.fitness_calls),
            })

        avalia_pop(iter_sma=0)
        if verbose:
            print(
                f"[SIBS] initial: fit={best_fit:.4f}, "
                f"verts={best_info['num_v']:.0f}, "
                f"cov_dem={best_info['cov_dem_rel']*100:.2f}%, "
                f"cov_ter={best_info['cov_ter_rel']*100:.2f}%"
            )

        for it in range(1, max_iters_sma + 1):
            if not self._can_evaluate(n_pop):
                stop_reason = "fitness_call_budget"
                break
            ordem = np.argsort(fitness_vals)
            pop = pop[ordem]
            fitness_vals = fitness_vals[ordem]
            infos = [infos[i] for i in ordem]

            elite = pop[0]
            second = pop[1] if n_pop > 1 else elite

            for idx in range(n_pop):
                if idx == 0:
                    continue
                base = elite if rng.rand() < 0.7 else second
                candidato = base.copy()
                n_flip = max(1, int(0.1 * J))
                flip_idx = rng.choice(J, size=n_flip, replace=False)
                candidato[flip_idx] = 1 - candidato[flip_idx]
                pop[idx] = candidato

            avalia_pop(iter_sma=it)

            if verbose and it % 10 == 0:
                print(
                    f"[SIBS] it={it}, fit={best_fit:.4f}, "
                    f"verts={best_info['num_v']:.0f}, "
                    f"cov_dem={best_info['cov_dem_rel']*100:.2f}%, "
                    f"cov_ter={best_info['cov_ter_rel']*100:.2f}%"
                )

            if it - last_improvement_sma >= patience_sma:
                stop_reason = "sma_patience"
                if verbose:
                    print(f"[SIBS] stopping: no improvement for {patience_sma} iterations.")
                break

        tempo_transicao_sma_vns = time.time() - t0
        sma_stopping_reason = stop_reason
        stop_reason = "maximum_iterations"

        if verbose:
            print("[VNS] starting variable neighborhood search...")

        best_sol_vns = best_sol.copy()
        best_fit_vns = best_fit
        best_info_vns = best_info

        budget_exhausted = False
        for it_vns in range(max_iters_vns):
            for k in range(1, k_max + 1):
                if not self._can_evaluate():
                    budget_exhausted = True
                    break
                shaken = self._shake(best_sol_vns, k)
                local_sol, local_fit, local_info = self._local_search(shaken)
                if local_fit < best_fit_vns:
                    best_sol_vns = local_sol.copy()
                    best_fit_vns = local_fit
                    best_info_vns = local_info
                    if verbose:
                        print(
                            f"[VNS] it={it_vns}, k={k}, new best: "
                            f"fit={best_fit_vns:.4f}, verts={best_info_vns['num_v']:.0f}, "
                            f"cov_dem={best_info_vns['cov_dem_rel']*100:.2f}%, "
                            f"cov_ter={best_info_vns['cov_ter_rel']*100:.2f}%"
                        )
                    break

            history.append({
                "fase": "VNS",
                "iter": it_vns,
                "tempo": time.time() - t0,
                "best_fit": best_fit_vns,
                "num_v": best_info_vns["num_v"],
                "cov_dem_rel": best_info_vns["cov_dem_rel"],
                "cov_ter_rel": best_info_vns["cov_ter_rel"],
                "cov_ter_svc": best_info_vns["cov_ter_svc"],
                "cov_ter_region": best_info_vns["cov_ter_region"],
                "fitness_calls": int(self.fitness_calls),
            })
            if budget_exhausted:
                stop_reason = "fitness_call_budget"
                break

        best_sol = best_sol_vns
        best_fit = best_fit_vns
        best_info = best_info_vns

        hist_df = pd.DataFrame(history)
        hist_df["tempo_transicao"] = tempo_transicao_sma_vns

        best_info = dict(best_info)
        best_info["fitness_calls"] = int(self.fitness_calls)
        best_info["stopping_reason"] = stop_reason
        best_info["sma_stopping_reason"] = sma_stopping_reason

        return best_sol, best_fit, best_info, hist_df


# ============================================================
# Scenarios
# ============================================================

CENARIOS = {
    "light": {
        "TIPO": "municipios",
        "CODIGOS": [COD_SAO_PAULO],
        "NOME": "São Paulo",
        "KAPPA": 0.005,
        "PHI": 0.005,  # deprecated compatibility alias
        "ETA_DEMANDA": 0.20,
        "ETA_TERRITORIO": 0.70,
        "RAIO_AERONAVE_KM": 10.0,
        "RAIO_SOLO_KM": 7.0,
        "RAIO_ELEV_KM": 6.0,
        "EXACT_TIME_LIMIT_SEC": 1800.0,
    },
    "medium": {
        "TIPO": "municipios",
        "CODIGOS": [COD_SAO_PAULO, COD_CAMPINAS, COD_GUARULHOS],
        "NOME": "SP + Campinas + Guarulhos",
        "KAPPA": 0.010,
        "PHI": 0.010,  # deprecated compatibility alias
        "ETA_DEMANDA": 0.25,
        "ETA_TERRITORIO": 0.75,
        "RAIO_AERONAVE_KM": 10.0,
        "RAIO_SOLO_KM": 6.0,
        "RAIO_ELEV_KM": 5.0,
        "EXACT_TIME_LIMIT_SEC": 1800.0,
    },
}

SMA_POP = 30
SMA_ITERS = 100
SMA_PATIENCE = 50
VNS_ITERS = 50
VNS_KMAX = 5

META_W_VERTS = 1.0
META_W_PREF = 0.5
META_LAMBDA_COV_DEM = 1e4
META_LAMBDA_COV_TER = 5e3
META_LAMBDA_SEP = 500.0

MIN_SEP_KM = 0.3


# ============================================================
# MAIN
# ============================================================

def main():
    print("=" * 70)
    print("Vertiport optimization: method comparison")
    print("Exact model vs SIBS vs VNS vs SIBS-VNS hybrid")
    print("=" * 70)
    
    print("Loading candidate sites...")
    helipontos_df = carregar_helipontos(CAMINHO_HELIPONTOS_XLS)

    resultados_exato = {}
    resultados_meta = {}

    for nome_cen, params in CENARIOS.items():
        print("\n" + "=" * 70)
        print(f"Scenario: {nome_cen.upper()} - {params['NOME']}")
        print("=" * 70)

        if params["TIPO"] == "municipios":
            demanda_df, mun_gdf = carregar_demanda_grade_ibge_1km(
                codigos_municipios=params["CODIGOS"]
            )
        elif params["TIPO"] == "estado":
            demanda_df, mun_gdf = carregar_demanda_grade_ibge_1km(
                estado=params["ESTADO"]
            )
        else:
            raise ValueError(f"Tipo desconhecido: {params['TIPO']}")

        geom_uniao = mun_gdf.unary_union
        hel_gdf = gpd.GeoDataFrame(
            helipontos_df.copy(),
            geometry=gpd.points_from_xy(
                helipontos_df["Longitude"], helipontos_df["Latitude"]
            ),
            crs="EPSG:4326",
        )
        hel_regiao = hel_gdf[hel_gdf.within(geom_uniao)].copy()
        hel_sub = pd.DataFrame(hel_regiao.drop(columns="geometry")).reset_index(drop=True)

        print(f"Candidate sites in the region: {len(hel_sub)}")
        print(f"Demand cells: {len(demanda_df)}")

        densidade_sub = calcular_densidade_populacional(demanda_df, hel_sub, raio_km=5.0)
        dist_vv_sub = compute_dist_vv(hel_sub)
        pref_score_sub = compute_preferencia_melhorada(
            hel_sub, demanda_df, densidade_sub, dist_vv_sub,
            w_densidade=1.0, w_tipo=0.5, w_overlap=0.3, raio_prox_km=3.0
        )

        reach_air = build_reach_matrix_aircraft(demanda_df, hel_sub, params["RAIO_AERONAVE_KM"])
        reach_fac = build_reach_matrix_facility(demanda_df, hel_sub, params["RAIO_SOLO_KM"], params["RAIO_ELEV_KM"])

        # EXATO
        res_exato = None
        if not META_ONLY:
            print("[EXACT] Building and solving the CFLP model...")
            exato = VertiportCFLP(
                hel_df=hel_sub,
                dem_df=demanda_df,
                reach_aircraft=reach_air,
                reach_facility=reach_fac,
                D_frac=params["KAPPA"],
                eta_demanda=params["ETA_DEMANDA"],
                eta_territorio=params["ETA_TERRITORIO"],
                w_open=1.0,
                pref_score_j=pref_score_sub,
                w_pref=META_W_PREF,
                dist_vv=dist_vv_sub,
                min_sep_km=MIN_SEP_KM,
                w_sep_penalty=META_LAMBDA_SEP,
                binary_assign=False,
            )
            exato.build(verbose=False)
            res_exato = exato.solve(time_limit=params["EXACT_TIME_LIMIT_SEC"], verbose=False)
            resultados_exato[nome_cen] = res_exato

        # Stochastic methods
        resultados_meta[nome_cen] = {
            'SMA': {},
            'VNS': {},
            'SMA+VNS': {}
        }

        for seed in SEEDS_EXPERIMENTAIS:
            print(f"[SEARCH] Running stochastic methods (seed={seed})...")
            
            # SMA-inspired search-only (legacy serialized ID: SMA)
            print("  SIBS...")
            sma = MetaheuristicaSMA(
                helipontos_df=hel_sub,
                demanda_df=demanda_df,
                reach_aircraft=reach_air,
                reach_facility=reach_fac,
                dist_vv=dist_vv_sub,
                D_frac=params["KAPPA"],
                eta_demanda=params["ETA_DEMANDA"],
                eta_territorio=params["ETA_TERRITORIO"],
                w_verts=META_W_VERTS,
                pref_score_j=pref_score_sub,
                w_pref=META_W_PREF,
                lambda_cov_dem=META_LAMBDA_COV_DEM,
                lambda_cov_ter=META_LAMBDA_COV_TER,
                lambda_sep=META_LAMBDA_SEP,
                min_sep_km=MIN_SEP_KM,
                seed=seed,
            )
            t0 = time.time()
            best_sol_sma, best_fit_sma, info_sma, hist_sma = sma.optimize(
                n_pop=SMA_POP,
                max_iters=SMA_ITERS,
                p_init_open=0.3,
                patience=SMA_PATIENCE,
                verbose=False,
            )
            t_sma = time.time() - t0
            
            resultados_meta[nome_cen]['SMA'][seed] = {
                "resultado": info_sma,
                "fitness": best_fit_sma,
                "tempo": t_sma,
                "x_sol": best_sol_sma,
                "hist": hist_sma,
            }
            
            # VNS puro
            print("  VNS...")
            vns = MetaheuristicaVNS(
                helipontos_df=hel_sub,
                demanda_df=demanda_df,
                reach_aircraft=reach_air,
                reach_facility=reach_fac,
                dist_vv=dist_vv_sub,
                D_frac=params["KAPPA"],
                eta_demanda=params["ETA_DEMANDA"],
                eta_territorio=params["ETA_TERRITORIO"],
                w_verts=META_W_VERTS,
                pref_score_j=pref_score_sub,
                w_pref=META_W_PREF,
                lambda_cov_dem=META_LAMBDA_COV_DEM,
                lambda_cov_ter=META_LAMBDA_COV_TER,
                lambda_sep=META_LAMBDA_SEP,
                min_sep_km=MIN_SEP_KM,
                seed=seed,
            )
            t0 = time.time()
            best_sol_vns, best_fit_vns, info_vns, hist_vns = vns.optimize(
                max_iters=VNS_ITERS,
                k_max=VNS_KMAX,
                p_init_open=0.3,
                verbose=False,
            )
            t_vns = time.time() - t0
            
            resultados_meta[nome_cen]['VNS'][seed] = {
                "resultado": info_vns,
                "fitness": best_fit_vns,
                "tempo": t_vns,
                "x_sol": best_sol_vns,
                "hist": hist_vns,
            }
            
            # SIBS-VNS hybrid (serialized ID: SMA+VNS)
            print("  SIBS-VNS hybrid...")
            sma_vns = MetaheuristicaSMA_VNS(
                helipontos_df=hel_sub,
                demanda_df=demanda_df,
                reach_aircraft=reach_air,
                reach_facility=reach_fac,
                dist_vv=dist_vv_sub,
                D_frac=params["KAPPA"],
                eta_demanda=params["ETA_DEMANDA"],
                eta_territorio=params["ETA_TERRITORIO"],
                w_verts=META_W_VERTS,
                pref_score_j=pref_score_sub,
                w_pref=META_W_PREF,
                lambda_cov_dem=META_LAMBDA_COV_DEM,
                lambda_cov_ter=META_LAMBDA_COV_TER,
                lambda_sep=META_LAMBDA_SEP,
                min_sep_km=MIN_SEP_KM,
                seed=seed,
            )
            t0 = time.time()
            best_sol_comb, best_fit_comb, info_comb, hist_comb = sma_vns.optimize(
                n_pop=SMA_POP,
                max_iters_sma=SMA_ITERS,
                max_iters_vns=VNS_ITERS,
                k_max=VNS_KMAX,
                p_init_open=0.3,
                patience_sma=SMA_PATIENCE,
                verbose=False,
            )
            t_comb = time.time() - t0
            
            resultados_meta[nome_cen]['SMA+VNS'][seed] = {
                "resultado": info_comb,
                "fitness": best_fit_comb,
                "tempo": t_comb,
                "x_sol": best_sol_comb,
                "hist": hist_comb,
            }
            
            print(
                f"    SIBS: fit={best_fit_sma:.2f}, verts={int(info_sma['num_v'])}, t={t_sma:.1f}s | "
                f"VNS: fit={best_fit_vns:.2f}, verts={int(info_vns['num_v'])}, t={t_vns:.1f}s | "
                f"SIBS-VNS hybrid: fit={best_fit_comb:.2f}, verts={int(info_comb['num_v'])}, t={t_comb:.1f}s"
            )

    # Save results
    import pickle
    arquivo_resultados = os.path.join(OUTPUT_DIR, 'resultados_completos.pkl')
    with open(arquivo_resultados, 'wb') as f:
        pickle.dump({
            'exato': resultados_exato,
            'meta': resultados_meta,
            'cenarios': CENARIOS,
        }, f)
    
    print("\n" + "="*70)
    print("Run complete.")
    print(f"Results saved to: {arquivo_resultados}")
    print("="*70)


if __name__ == "__main__":
    main()
