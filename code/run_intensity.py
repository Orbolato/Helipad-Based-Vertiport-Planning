"""Complete a paired five-seed, three-level demand-intensity experiment."""
from analysis_common import *
import platform,sys

def main():
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--scenario',choices=['light','medium','heavy'])
    parser.add_argument('--seeds',default=','.join(map(str,SEEDS)))
    parser.add_argument('--rerun',action='store_true')
    args=parser.parse_args()
    selected_seeds=[int(v) for v in args.seeds.split(',')]
    for scenario in ['heavy','light','medium']:
        if args.scenario and args.scenario!=scenario:continue
        print('Loading',scenario,flush=True)
        inst=load_instance(scenario)
        for factor in [.8,1.,1.2]:
            for seed in selected_seeds:
                stem=f'{scenario}_intensity_{factor:.1f}_seed_{seed}'
                path=ROOT/'results'/f'{stem}.json'
                if path.exists() and path.with_suffix('.npz').exists() and not args.rerun:
                    # Reverify preserved site vectors against the extracted inputs.
                    old=json.loads(path.read_text(encoding='utf-8'))
                    x=np.load(path.with_suffix('.npz'))['x']
                    e=engine(inst,seed,factor)
                    cert=verify(inst,x,e.D,e.C,stem)
                    assert abs(cert['summary']['objective']-old['objective'])<1e-6
                    assert cert['summary']['feasible']==old['feasible'],(stem,'original feasibility mismatch')
                    # Maximum allocation can exceed a previously recorded feasible allocation.
                    write_record(stem,cert['summary'],scenario=scenario,experiment='intensity',factor=factor,seed=seed,source=old['source'],fitness_calls=old['fitness_calls'])
                    del e
                    continue
                cert,meta=run_relay(inst,seed,factor,stem)
                write_record(stem,cert['summary'],scenario=scenario,experiment='intensity',factor=factor,seed=seed,source='five-seed extension',**meta)
        del inst
    records=[json.loads(p.read_text(encoding='utf-8')) for p in (ROOT/'results').glob('*_intensity_*_seed_*.json')]
    assert len(records)==45
    frame=pd.DataFrame(records)
    frame.to_csv(ROOT/'results/intensity_seed_results.csv',index=False)
    out=[]
    for (s,f),group in frame.groupby(['scenario','factor']):
        feasible=group[group.feasible]
        out.append(dict(scenario=s,factor=f,runs=len(group),feasible=len(feasible),sites_all_mean=group.facilities.mean(),sites_all_sd=group.facilities.std(ddof=1),sites_feasible_mean=feasible.facilities.mean(),sites_feasible_sd=feasible.facilities.std(ddof=1),minimum=int(group.facilities.min()),maximum=int(group.facilities.max())))
    pd.DataFrame(out).to_csv(ROOT/'results/intensity_summary.csv',index=False)
    (ROOT/'results/intensity_environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'timing_use':'Additional run times are recorded for reproducibility only, not compared with the original timing experiment.'},indent=2))
    print(pd.DataFrame(out).to_string(index=False),flush=True)

if __name__=='__main__':main()
