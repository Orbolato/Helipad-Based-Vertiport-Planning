"""Independently check saved allocations, feasibility classifications, and summaries."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import geopandas as gpd
from analysis_common import ROOT,load_instance,engine,opt
from run_spatial import refined

def main():
    count=0;classified={};tol=1e-7
    for scenario in ['light','medium','heavy']:
        inst=load_instance(scenario);base=engine(inst);C=base.C.copy();D=base.D.copy();A=base.A.copy();del base
        variants={1000:inst}
        if scenario!='heavy':
            variants[500]=refined(inst,2);variants[250]=refined(inst,4)
            grid=gpd.read_file(ROOT/'data'/f'{scenario}_grid.gpkg').to_crs(31983)
            boundary=gpd.read_file(ROOT/'data'/f'{scenario}_boundaries.gpkg').to_crs(31983).union_all()
            fraction=(grid.geometry.intersection(boundary).area/grid.area).clip(0,1).to_numpy()
        else:
            grid=gpd.read_file(ROOT/'data/heavy_grid.gpkg').to_crs(4326);h=np.zeros(len(grid))
            for k in range(4):
                corners=np.array([list(g.exterior.coords)[k] for g in grid.geometry])
                h=np.maximum(h,opt.haversine_km(inst['dem'].Latitude.to_numpy(),inst['dem'].Longitude.to_numpy(),corners[:,1],corners[:,0]))
            h+=.001
        for path in sorted((ROOT/'results').glob(scenario+'_*.json')):
            d=json.loads(path.read_text(encoding='utf-8'))
            if 'certificate_sha256' not in d or 'experiment' not in d:continue
            z=np.load(path.with_suffix('.npz'));x=z['x'];qi=z['allocation_i'];qj=z['allocation_j'];q=z['allocation_q']
            experiment=d['experiment'];current=inst;dem=D.copy();area=A.copy()
            if experiment=='intensity':dem*=d['factor']
            elif experiment.startswith('composition'):
                dem=pd.read_csv(ROOT/'results'/f"{scenario}_composition_{d['pattern']}_demand.csv").planning_demand.to_numpy()
            elif experiment.startswith('resolution'):
                current=variants[d['resolution_m']];dem=inst['params']['KAPPA']*current['dem'].Populacao.to_numpy();area=current['dem'].Area_km2.to_numpy()
            elif experiment=='boundary_fixed':dem*=fraction;area*=fraction
            assert np.isin(x,[0,1]).all() and np.isfinite(q).all() and (q>=0).all(),path.name
            assert ((qi>=0)&(qi<len(dem))).all() and ((qj>=0)&(qj<len(C))).all(),path.name
            cell=np.bincount(qi,weights=q,minlength=len(dem));load=np.bincount(qj,weights=q,minlength=len(C))
            access=current['fac'][qi,qj]
            if experiment=='corner_envelope':
                distances=opt.haversine_km(inst['dem'].Latitude.to_numpy()[qi],inst['dem'].Longitude.to_numpy()[qi],inst['hel'].Latitude.to_numpy()[qj],inst['hel'].Longitude.to_numpy()[qj])
                radius=np.array([5. if 'solo' in str(t).lower() else 4. for t in inst['hel'].Tipo])[qj]
                access=(distances+h[qi]<=radius) if d['kind']=='contracted' else (distances-h[qi]<=radius)
                assert d['reach_ratio_lower']<=d['candidate_reachable_territorial_coverage']<=d['reach_ratio_upper']
            covered=current['air'][:,x.astype(bool)].any(axis=1);reachable=current['air'].any(axis=1)
            demand_margin=q.sum()-inst['params']['ETA_DEMANDA']*dem.sum()
            area_margin=area[covered].sum()-inst['params']['ETA_TERRITORIO']*area[reachable].sum()
            feasible=bool((cell-dem).max()<=tol and (load-C*x).max()<=tol and access.all() and (x[qj]>0).all() and demand_margin>=-tol and area_margin>=-tol)
            assert feasible==d['feasible'],(path.name,feasible,demand_margin,area_margin)
            np.testing.assert_allclose(q.sum(),d['allocated_demand'],atol=1e-6,rtol=1e-12)
            np.testing.assert_allclose(q.sum()/dem.sum(),d['demand_coverage'],atol=1e-12)
            np.testing.assert_allclose(area[covered].sum()/area[reachable].sum(),d['candidate_reachable_territorial_coverage'],atol=1e-12)
            objective=x.sum()+.5*np.dot(x,inst['pref'])+500*opt.spacing_diagnostics(x,inst['vv'],.3)['spacing_penalty_raw']
            np.testing.assert_allclose(objective,d['objective'],rtol=1e-12,atol=1e-8)
            digest=hashlib.sha256()
            for a in [x.astype(np.int8),qi.astype('<i8'),qj.astype('<i8'),q.astype('<f8')]:digest.update(np.ascontiguousarray(a).tobytes())
            assert digest.hexdigest()==d['certificate_sha256'],path.name
            count+=1;classified[experiment]=classified.get(experiment,0)+1
        print(scenario,'allocation checks passed',flush=True)
    assert count==147,(count,classified)
    summary=pd.read_csv(ROOT/'results/intensity_summary.csv')
    runs=pd.DataFrame([json.loads(p.read_text()) for p in (ROOT/'results').glob('*_intensity_*_seed_*.json')])
    assert len(runs)==45 and runs.feasible.sum()==44
    for (scenario,factor),g in runs.groupby(['scenario','factor']):
        row=summary[(summary.scenario==scenario)&(abs(summary.factor-factor)<1e-9)].iloc[0]
        np.testing.assert_allclose([row.sites_all_mean,row.sites_all_sd],[g.facilities.mean(),g.facilities.std(ddof=1)],rtol=1e-12)
        assert row.feasible==g.feasible.sum()
    report={'allocations_checked':count,'by_experiment':classified,'intensity_runs':45,'intensity_feasible':44,'native_tolerance':tol,'status':'passed'}
    (ROOT/'validation_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
