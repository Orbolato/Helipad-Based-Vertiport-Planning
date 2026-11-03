from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'
TABLES=ROOT/'tables'
SCENARIOS=['light','medium','heavy']
LABELS={'light':'Municipal','medium':'Regional','heavy':'Statewide'}
def read(stem):return json.loads((RESULTS/(stem+'.json')).read_text(encoding='utf-8'))
def write(name,value):(TABLES/name).write_text(value,encoding='utf-8')
def table(caption,label,cols,head,rows,notes=''):
    return '\n'.join([r'\begin{table}[htbp]',r'\centering',r'\caption{'+caption+'}',r'\label{'+label+'}',r'\small',r'\begin{tabular}{@{}'+cols+'@{}}',r'\toprule',head+r' \\',r'\midrule',*[' & '.join(row)+r' \\' for row in rows],r'\bottomrule',r'\end{tabular}',r'\par\footnotesize '+notes,r'\end{table}',''])
def main():
    import pandas as pd,numpy as np
    TABLES.mkdir(exist_ok=True)
    if (RESULTS/'intensity_summary.csv').exists():
        frame=pd.read_csv(RESULTS/'intensity_summary.csv');rows=[]
        for s in SCENARIOS:
            row=[LABELS[s]]
            for factor in [.8,1.,1.2]:
                d=frame[(frame.scenario==s)&np.isclose(frame.factor,factor)].iloc[0]
                row.append(f"{d.sites_all_mean:.2f} $\\pm$ {d.sites_all_sd:.2f} ({int(d.feasible)}/5)")
            rows.append(row)
        intro=r'Table~\ref{tab:sensitivity} summarizes all five seeds at each demand intensity. '
        pairs=[]
        for s in ['medium','heavy']:
            g=frame[frame.scenario==s].sort_values('factor')
            pairs.append(f"{LABELS[s].lower()} mean counts increased from {g.iloc[0].sites_all_mean:.2f} to {g.iloc[-1].sites_all_mean:.2f}")
        intro+=('; '.join(pairs)).capitalize()+'.\n'
        total=int(frame.feasible.sum())
        intro+=f'Of the 45 configurations, {total} satisfied both coverage floors after continuous allocation. '
        if total<45:intro+='Counts for infeasible runs are retained in the all-run summaries and are not interpreted as feasible infrastructure requirements. '
        intro+='The high-intensity statewide run with seed 1011 selected 180 sites and attained 79.9928\\% candidate-reachable coverage; its territorial shortfall was 3.268 km$^2$. Among the four feasible high-intensity statewide runs, the count was 182.50 $\\pm$ 3.70. '
        intro+='These differences describe sensitivity to synthetic demand magnitude rather than forecast adoption.\n\n'
        write('intensity_results.tex',intro+table('Demand-intensity sensitivity: selected sites, mean $\\pm$ sample SD (feasible-run counts)','tab:sensitivity','lccc',r'Scenario & $0.8\kappa_s$ & $\kappa_s$ & $1.2\kappa_s$',rows,r'All five seeds contribute to each mean and standard deviation. The infeasible statewide run is retained in the all-run summary.'))
    rows=[]
    for s in SCENARIOS:
        for pattern in ['concentrated','dispersed']:
            fixed=[read(f'{s}_composition_{pattern}_fixed_relay_{seed}') for seed in [42,123,456,789,1011]]
            rerun=read(f'{s}_composition_{pattern}_reoptimized')
            exact='--' if s=='heavy' else ('Yes' if read(f'{s}_composition_{pattern}_fixed_exact_0')['feasible'] else 'No')
            rows.append([LABELS[s],pattern.capitalize(),exact,f"{sum(d['feasible'] for d in fixed)}/5",str(rerun['facilities']),f"{rerun['site_jaccard']:.3f}"])
    gaps=[read(f'medium_composition_{p}_reoptimized')['mip_gap']*100 for p in ['concentrated','dispersed']]
    body=table('Spatial-demand sensitivity at constant total demand and nominal capacity','tab:composition',r'llccc@{\hspace{1em}}c',r'Scenario & Demand & Exact set & Relay sets & Sites & Jaccard',rows,r'Exact set and relay sets report feasibility of unchanged nominal configurations. Sites and Jaccard refer to reoptimized configurations; Jaccard compares with the nominal exact set in the two smaller cases and nominal relay seed 42 statewide. The regional dispersed-demand incumbent has a '+f'{gaps[1]:.2f}'+r'\% optimality gap; the other municipal/regional composition solves reached zero gap. Statewide reoptimization uses one relay run per pattern and has no optimality guarantee.')
    heavy=[read(f'heavy_composition_{p}_reoptimized') for p in ['concentrated','dispersed']]
    body+=r'''
The four-site municipal and 13-site regional exact configurations remained
feasible under both spatial redistributions. In contrast, every nominal
regional and statewide relay configuration failed the demand floor under
concentrated demand. The statewide maximum allocated-demand fractions were
28.01--28.41\%, below the 30\% requirement. Under dispersed demand, all five
regional configurations and three of five statewide configurations remained
feasible. Feasibility therefore depends on the selected sites and spatial
demand composition, even when aggregate demand and nominal capacities are
unchanged.

'''
    body+=f"Statewide reoptimization selected {heavy[0]['facilities']} and {heavy[1]['facilities']} sites for concentrated and dispersed demand, respectively, with Jaccard indices {heavy[0]['site_jaccard']:.3f} and {heavy[1]['site_jaccard']:.3f} relative to the nominal seed-42 selection. "
    body+=r'''These comparisons are conditional on one stochastic run per redistribution.
For context, the ten pairwise Jaccard indices among the five nominal relay
runs range from 0.321 to 0.382 statewide, from 0 to 0.222 municipally, and
from 0.057 to 0.269 regionally. Site variation therefore also arises from
stochastic search and alternative configurations. Similar site counts do not
establish location stability, and the experiment does not attribute all
changes in reoptimized sites exclusively to demand composition.
'''
    write('composition_results.tex',body)
    rows=[]
    for s in ['light','medium']:
        for resolution in [1000,500,250]:
            d=read(f'{s}_resolution_{resolution}m_fixed_exact_0')
            fixed=[read(f'{s}_resolution_{resolution}m_fixed_relay_{seed}') for seed in [42,123,456,789,1011]]
            rows.append([LABELS[s],str(resolution),str(d['facilities']),f"{100*d['demand_coverage']:.2f}",f"{100*d['candidate_reachable_territorial_coverage']:.2f}",f"{sum(v['feasible'] for v in fixed)}/5"])
    gaps=[read(f'medium_resolution_{n}m_reoptimized')['mip_gap']*100 for n in [500,250]]
    body=table('Spatial-resolution sensitivity of unchanged nominal site configurations','tab:resolution','lrrrrc',r'Scenario & Grid (m) & Exact sites & Demand (\%) & Reach (\%) & Relay feasible',rows,r'Demand is the maximum continuous allocation for the nominal exact set; reach is candidate-relative. Both exact sets remain feasible. Reoptimization retains the same municipal and regional site sets at both finer resolutions. Municipal solves reached zero gap; regional 500 m and 250 m incumbents have gaps of '+f'{gaps[0]:.2f}'+r'\% and '+f'{gaps[1]:.2f}'+r'\%, respectively, at the 300 s limit.')
    body+=r'''
The nominal exact configurations remained feasible at 500 m and 250 m.
Their candidate-reachable coverage changed from 70.14\% to 70.36\%
municipally and from 75.52\% to 75.27\% regionally between the 1 km and
250 m representations. However, relay seed 1011 fell below the territorial
floor at both finer resolutions in both cases. At 250 m its coverage was
69.95\% municipally and 74.89\% regionally. Small changes in coverage can
therefore change feasibility when a configuration lies close to a threshold.
The retained exact site sets show conditional stability in these experiments;
the nonzero regional reoptimization gaps preclude asserting their optimality
on the finer grids.
'''
    write('resolution_results.tex',body)
    rows=[]
    for s in ['light','medium']:
        d=read(f'{s}_boundary_fixed_exact_0')
        allruns=[d]+[read(f'{s}_boundary_fixed_relay_{seed}') for seed in [42,123,456,789,1011]]
        rows.append([LABELS[s]+' boundary',f"{100*d['demand_coverage']:.2f}",f"{100*d['candidate_reachable_territorial_coverage']:.2f}",f"{sum(x['feasible'] for x in allruns)}/6"])
    for kind in ['contracted','expanded']:
        ds=[read(f'heavy_corner_envelope_{kind}_seed_{seed}') for seed in [42,123,456,789,1011]]
        key='reach_ratio_lower' if kind=='contracted' else 'reach_ratio_upper'
        rows.append(['Statewide '+kind,f"{100*min(d['demand_coverage'] for d in ds):.2f}--{100*max(d['demand_coverage'] for d in ds):.2f}",f"{100*min(d[key] for d in ds):.2f}--{100*max(d[key] for d in ds):.2f}",f"{sum(d['feasible'] for d in ds)}/5"])
    body=table('Boundary-support and statewide distance-envelope sensitivity','tab:spatial_support','lccc',r'Test & Demand (\%) & Reach (\%) & Feasible',rows,r'Boundary rows report exact-set coverage after weighting cell mass and area by the boundary-intersection fraction; feasibility counts include that set and five relay sets. Statewide demand ranges use contracted or expanded access links. Reach ranges use lower or upper territorial ratio envelopes, respectively. Statewide feasibility counts concern access-link perturbations with nominal territorial coverage; the lower reach envelopes also fall below the 80\% floor.')
    h=max(read(f'heavy_corner_envelope_contracted_seed_{seed}')['max_corner_distance_km'] for seed in [42,123,456,789,1011])
    body+=r'''
Boundary weighting retained feasibility for the nominal exact and relay
configurations in both smaller cases, but changed their coverage values.
These values depend on the assumed within-cell population distribution and
should not be read as corrected municipal demand estimates. Statewide,
contracted access links reduced maximum demand allocation below the 30\%
floor for all five nominal relay configurations, while expanded links retained
feasibility. The territorial envelopes span the 80\% threshold, so the
centroid-model outcomes do not establish feasibility for every within-cell
placement of demand or area. The maximum centroid-to-corner distance including
the 1 m allowance was '''+f'{h:.3f}'+r''' km. The envelopes intentionally combine
unfavourable numerator and denominator cases and are not estimates of expected
coverage error.
'''
    write('spatial_support_results.tex',body)
    allrecords=[json.loads(p.read_text(encoding='utf-8')) for p in RESULTS.glob('*.json') if '_seed_' in p.stem or '_fixed_' in p.stem or '_reoptimized' in p.stem]
    pd.DataFrame(allrecords).to_csv(RESULTS/'all_analysis_records.csv',index=False)
    print('Manuscript tables generated')
if __name__=='__main__':main()
