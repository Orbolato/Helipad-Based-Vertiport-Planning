"""Controlled spatial-demand, grid-support, and selected-site analyses."""
from analysis_common import *
import geopandas as gpd
from pyproj import Transformer
from copy import copy
import argparse,gc
ALBERS='+proj=aea +lat_1=-2 +lat_2=-22 +lat_0=-12 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=GRS80 +units=m +no_defs'

def originals(scenario):
    result=[]
    if scenario!='heavy':
        x=np.load(ROOT/'data/original_solutions'/scenario/'exact_optimum.npz')['x']
        result.append(('exact',0,x))
    for seed in SEEDS:
        x=np.load(ROOT/'results'/f'{scenario}_intensity_1.0_seed_{seed}.npz')['x']
        result.append(('relay',seed,x))
    return result

def exact(inst,C,stem,start=None):
    p=inst['params']
    m=opt.VertiportCFLP(hel_df=inst['hel'],dem_df=inst['dem'],reach_aircraft=inst['air'],reach_facility=inst['fac'],D_frac=p['KAPPA'],eta_demanda=p['ETA_DEMANDA'],eta_territorio=p['ETA_TERRITORIO'],w_open=1,pref_score_j=inst['pref'],w_pref=.5,dist_vv=inst['vv'],min_sep_km=.3,w_sep_penalty=500)
    m.C=C.copy();m.build(verbose=False)
    m.model.Params.Threads=2
    m.model.Params.FeasibilityTol=1e-9
    m.model.Params.IntFeasTol=1e-9
    m.model.Params.MIPGap=0
    if start is not None:
        for j,v in enumerate(start):m.x[j].Start=int(v)
    result=m.solve(time_limit=300,verbose=False)
    if result['x'] is None:raise RuntimeError('No incumbent: '+stem)
    cert=verify(inst,result['x'],p['KAPPA']*inst['dem'].Populacao.to_numpy(),C,stem)
    meta=dict(solver_status=int(m.model.Status),mip_gap=float(m.model.MIPGap),objective_bound=float(m.model.ObjBound),solver_seconds=float(m.model.Runtime))
    m.model.dispose()
    return cert,meta

def composition(scenario):
    inst=load_instance(scenario);base=engine(inst);C=base.C.copy();del base
    pop=inst['dem'].Populacao.to_numpy();density=pop/inst['dem'].Area_km2.to_numpy()
    order=np.argsort(density);threshold=density[order[np.searchsorted(np.cumsum(pop[order]),pop.sum()/2)]]
    dense=density>=threshold
    for pattern in ['concentrated','dispersed']:
        weight=np.where(dense,1.5,.5)
        if pattern=='dispersed':weight=2-weight
        adjusted=pop*weight*pop.sum()/np.dot(pop,weight)
        assert np.isclose(adjusted.sum(),pop.sum(),rtol=1e-14)
        alt=copy(inst);alt['dem']=inst['dem'].copy();alt['dem']['Populacao']=adjusted
        pd.DataFrame({'cell_index':np.arange(len(pop)),'weight':weight,'planning_demand':adjusted*inst['params']['KAPPA']}).to_csv(ROOT/'results'/f'{scenario}_composition_{pattern}_demand.csv',index=False)
        for method,seed,x in originals(scenario):
            stem=f'{scenario}_composition_{pattern}_fixed_{method}_{seed}'
            cert=verify(alt,x,adjusted*inst['params']['KAPPA'],C,stem)
            write_record(stem,cert['summary'],scenario=scenario,experiment='composition_fixed',pattern=pattern,method=method,seed=seed,population_weighted_density_median=float(threshold),redistributed_demand_fraction=float(.5*np.abs(adjusted-pop).sum()/pop.sum()))
        stem=f'{scenario}_composition_{pattern}_reoptimized'
        if scenario=='heavy':
            cert,meta=run_relay(alt,42,1.,stem,capacities=C);baseline=originals(scenario)[0][2]
            meta['method']='relay';meta['seed']=42
        else:
            baseline=originals(scenario)[0][2];cert,meta=exact(alt,C,stem,baseline);meta['method']='exact'
        x=cert['x'];union=np.count_nonzero(x|baseline);intersection=np.count_nonzero(x&baseline)
        write_record(stem,cert['summary'],scenario=scenario,experiment='composition_reoptimized',pattern=pattern,site_jaccard=intersection/union,retained_sites=int(intersection),baseline_sites=int(baseline.sum()),**meta)
        del alt;gc.collect()

def refined(inst,subdivision):
    grid=gpd.read_file(ROOT/'data'/f"{inst['scenario']}_grid.gpkg").to_crs(ALBERS)
    b=grid.bounds
    xy=[]
    for a in range(subdivision):
        for c in range(subdivision):
            xy.append((b.minx.to_numpy()+(a+.5)*(b.maxx-b.minx).to_numpy()/subdivision,b.miny.to_numpy()+(c+.5)*(b.maxy-b.miny).to_numpy()/subdivision))
    x=np.stack([v[0] for v in xy],axis=1).ravel();y=np.stack([v[1] for v in xy],axis=1).ravel()
    lon,lat=Transformer.from_crs(ALBERS,4326,always_xy=True).transform(x,y)
    count=subdivision**2
    dem=pd.DataFrame({'Longitude':lon,'Latitude':lat,'Populacao':np.repeat(inst['dem'].Populacao.to_numpy()/count,count),'Area_km2':np.repeat(inst['dem'].Area_km2.to_numpy()/count,count)})
    alt=copy(inst);alt['dem']=dem
    p=inst['params'];alt['air']=opt.build_reach_matrix_aircraft(dem,inst['hel'],10)
    alt['fac']=opt.build_reach_matrix_facility(dem,inst['hel'],p['RAIO_SOLO_KM'],p['RAIO_ELEV_KM'])
    assert np.isclose(dem.Populacao.sum(),inst['dem'].Populacao.sum(),rtol=1e-13)
    assert np.isclose(dem.Area_km2.sum(),inst['dem'].Area_km2.sum(),rtol=1e-13)
    return alt

def resolution(scenario):
    inst=load_instance(scenario);base=engine(inst);C=base.C.copy();del base
    # Native centroid baseline and deterministic subcell quadrature.
    for n in [1,2,4]:
        alt=inst if n==1 else refined(inst,n)
        for method,seed,x in originals(scenario):
            stem=f'{scenario}_resolution_{1000//n}m_fixed_{method}_{seed}'
            cert=verify(alt,x,alt['dem'].Populacao.to_numpy()*inst['params']['KAPPA'],C,stem)
            write_record(stem,cert['summary'],scenario=scenario,experiment='resolution_fixed',resolution_m=1000//n,method=method,seed=seed,cells=len(alt['dem']))
        if n>1:
            stem=f'{scenario}_resolution_{1000//n}m_reoptimized'
            baseline=originals(scenario)[0][2]
            cert,meta=exact(alt,C,stem,baseline)
            x=cert['x']
            write_record(stem,cert['summary'],scenario=scenario,experiment='resolution_reoptimized',resolution_m=1000//n,site_jaccard=float(np.count_nonzero(x&baseline)/np.count_nonzero(x|baseline)),**meta)
        if n>1:del alt;gc.collect()
    # Clip boundary-cell mass and area under uniform within-cell density.
    grid=gpd.read_file(ROOT/'data'/f'{scenario}_grid.gpkg').to_crs(31983)
    boundary=gpd.read_file(ROOT/'data'/f'{scenario}_boundaries.gpkg').to_crs(31983).union_all()
    clipped=grid.geometry.intersection(boundary)
    fraction=(clipped.area/grid.area).clip(0,1).to_numpy()
    alt=copy(inst);alt['dem']=inst['dem'].copy()
    alt['dem']['Populacao']*=fraction;alt['dem']['Area_km2']*=fraction
    for method,seed,x in originals(scenario):
        stem=f'{scenario}_boundary_fixed_{method}_{seed}'
        cert=verify(alt,x,alt['dem'].Populacao.to_numpy()*inst['params']['KAPPA'],C,stem)
        write_record(stem,cert['summary'],scenario=scenario,experiment='boundary_fixed',method=method,seed=seed,retained_population=float(alt['dem'].Populacao.sum()),boundary_cells=int(np.count_nonzero(fraction<1-1e-8)))

def statewide_envelope():
    scenario='heavy';inst=load_instance(scenario);base=engine(inst);C=base.C.copy();D=base.D.copy();del base
    grid=gpd.read_file(ROOT/'data/heavy_grid.gpkg').to_crs(4326)
    h=np.zeros(len(grid))
    # Largest Haversine centroid-to-vertex distance plus 1 m guard.
    for k in range(4):
        coords=np.array([list(g.exterior.coords)[k] for g in grid.geometry])
        h=np.maximum(h,opt.haversine_km(inst['dem'].Latitude.to_numpy(),inst['dem'].Longitude.to_numpy(),coords[:,1],coords[:,0]))
    h+=.001
    nearest=np.full(len(grid),np.inf);selected_nearest={seed:np.full(len(grid),np.inf) for seed in SEEDS}
    xs={seed:x for method,seed,x in originals(scenario)}
    fac_minus=np.zeros_like(inst['fac']);fac_plus=np.zeros_like(inst['fac'])
    for j,row in inst['hel'].iterrows():
        dist=opt.haversine_km(inst['dem'].Latitude.to_numpy(),inst['dem'].Longitude.to_numpy(),row.Latitude,row.Longitude)
        nearest=np.minimum(nearest,dist)
        radius=5 if 'solo' in str(row.Tipo).lower() else 4
        fac_minus[:,j]=dist+h<=radius;fac_plus[:,j]=dist-h<=radius
        for seed in SEEDS:
            if xs[seed][j]:selected_nearest[seed]=np.minimum(selected_nearest[seed],dist)
    A=inst['dem'].Area_km2.to_numpy()
    for kind,fac in [('contracted',fac_minus),('expanded',fac_plus)]:
        alt=copy(inst);alt['fac']=fac
        for seed in SEEDS:
            stem=f'heavy_corner_envelope_{kind}_seed_{seed}'
            cert=verify(alt,xs[seed],D,C,stem)
            # Ratio envelope uses opposite numerator and denominator extremes.
            d=selected_nearest[seed]
            low=A[d+h<=10].sum()/A[nearest-h<=10].sum()
            high=min(1.,A[d-h<=10].sum()/A[nearest+h<=10].sum())
            write_record(stem,cert['summary'],scenario=scenario,experiment='corner_envelope',kind=kind,seed=seed,reach_ratio_lower=float(low),reach_ratio_upper=float(high),max_corner_distance_km=float(h.max()),candidate_boundary_area_fraction=float(A[np.abs(nearest-10)<=h].sum()/A[nearest<=10].sum()))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('analysis',choices=['composition','resolution','envelope']);parser.add_argument('--scenario',default='light');args=parser.parse_args()
    if args.analysis=='composition':composition(args.scenario)
    elif args.analysis=='resolution':resolution(args.scenario)
    else:statewide_envelope()
