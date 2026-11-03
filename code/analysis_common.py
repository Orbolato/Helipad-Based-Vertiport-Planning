"""Shared construction and independent continuous allocation verification."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ.setdefault(name,'1')
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
import vertiport_location_optimizer as opt
ROOT=Path(__file__).resolve().parents[1]
SEEDS=[42,123,456,789,1011]

def load_instance(scenario):
    dem=pd.read_csv(ROOT/'data'/f'{scenario}_cells.csv',float_precision='round_trip')
    hel=pd.read_csv(ROOT/'data'/f'{scenario}_candidates.csv',float_precision='round_trip')
    params=json.loads((ROOT/'data'/f'{scenario}_parameters.json').read_text(encoding='utf-8'))
    fixed=np.load(ROOT/'data'/f'{scenario}_fixed_inputs.npz')
    air=opt.build_reach_matrix_aircraft(dem,hel,10.)
    fac=opt.build_reach_matrix_facility(dem,hel,params['RAIO_SOLO_KM'],params['RAIO_ELEV_KM'])
    return dict(scenario=scenario,dem=dem,hel=hel,params=params,air=air,fac=fac,vv=fixed['dist_vv'],pref=fixed['pref'])

def engine(inst,seed=42,factor=1.):
    p=inst['params']
    return opt.MetaheuristicaSMA_VNS(helipontos_df=inst['hel'],demanda_df=inst['dem'],reach_aircraft=inst['air'],reach_facility=inst['fac'],dist_vv=inst['vv'],D_frac=p['KAPPA']*factor,eta_demanda=p['ETA_DEMANDA'],eta_territorio=p['ETA_TERRITORIO'],w_verts=1.,pref_score_j=inst['pref'],w_pref=.5,lambda_cov_dem=1e4,lambda_cov_ter=5e3,lambda_sep=500.,min_sep_km=.3,seed=seed,fitness_call_budget=15780 if inst['scenario']=='heavy' else None)

def verify(inst,x,D,C,stem=None):
    """Maximize served demand with fixed sites, then verify native residuals."""
    selected=np.flatnonzero(x)
    ii,jj=np.where(inst['fac'][:,selected] & (D[:,None]>0))
    jj=selected[jj]
    cells,inv=np.unique(ii,return_inverse=True)
    n=len(ii); nc=len(cells); nj=len(x)
    if n:
        mat=coo_matrix((np.ones(2*n),(np.r_[inv,nc+jj],np.r_[np.arange(n),np.arange(n)])),shape=(nc+nj,n)).tocsr()
        lp=linprog(-np.ones(n),A_ub=mat,b_ub=np.r_[D[cells],C*x],bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9})
        if not lp.success:raise RuntimeError(lp.message)
        q=lp.x
    else:q=np.zeros(0)
    cert=opt.verify_allocation_certificate(x_sol=x,allocation_i=ii,allocation_j=jj,allocation_q=q,demand=D,capacities=C,area_km2=inst['dem'].Area_km2.to_numpy(),reach_facility=inst['fac'],reach_aircraft=inst['air'],dist_vv=inst['vv'],eta_demanda=inst['params']['ETA_DEMANDA'],eta_territorio=inst['params']['ETA_TERRITORIO'],min_sep_km=.3,pref_score_j=inst['pref'],w_open=1,w_pref=.5,w_sep_penalty=500.)
    if stem:opt.save_feasibility_certificate(cert,str(ROOT/'results'/stem))
    return cert

def run_relay(inst,seed,factor,stem,capacities=None):
    e=engine(inst,seed,factor)
    if capacities is not None:e.C=capacities.copy()
    start=time.perf_counter()
    x,fit,info,hist=e.optimize(n_pop=30,max_iters_sma=100,max_iters_vns=15780 if inst['scenario']=='heavy' else 50,k_max=5,p_init_open=.3,patience_sma=50,verbose=False)
    seconds=time.perf_counter()-start
    cert=verify(inst,x,e.D,e.C,stem)
    return cert,dict(fitness=float(fit),fitness_calls=int(info['fitness_calls']),optimization_seconds=seconds)

def write_record(stem,summary,**metadata):
    record=dict(summary,**metadata)
    (ROOT/'results'/f'{stem}.json').write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')
    print(stem,record['facilities'],record['feasible'],flush=True)
