from pathlib import Path
import pandas as pd,json
ROOT=Path(__file__).resolve().parents[1]
R=ROOT/'results'
def esc(t):return str(t).replace('&',r'\&').replace('%',r'\%').replace('_',r'\_')
def longtable(cols,head,rows):
    return '\n'.join([r'\begin{longtable}{'+cols+'}',r'\toprule',head+r'\\\midrule\endhead',*[' & '.join(map(esc,row))+r' \\' for row in rows],r'\bottomrule',r'\end{longtable}'])
def main():
    parts=[r'''\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=0.85in]{geometry}
\usepackage{booktabs,longtable,array}
\usepackage[hidelinks]{hyperref}
\usepackage{microtype}
\setlength{\emergencystretch}{3em}
\hypersetup{pdftitle={Supplementary information: CSTP-D-26-01291},pdfauthor={}}
\begin{document}
\begin{center}\Large\textbf{Supplementary information}\\[4pt]
\normalsize CSTP-D-26-01291\end{center}

\section{Demand-intensity sensitivity}

Each scenario was evaluated at three demand-intensity levels with seeds 42,
123, 456, 789, and 1011. Candidate capacities and all other planning parameters
were fixed within each scenario. Coverage was assessed after maximizing
continuous allocation for the selected sites. Territorial coverage is relative
to the area reachable from all candidates. Mean and sample SD are calculated
from all five counts; infeasible configurations remain in the record.

\small
''']
    records=[json.loads(p.read_text()) for p in R.glob('*_intensity_*_seed_*.json')]
    f=pd.DataFrame(records)
    rows=[]
    names={'light':'Municipal','medium':'Regional','heavy':'Statewide'}
    for s in ['light','medium','heavy']:
        for _,d in f[f.scenario==s].sort_values(['factor','seed']).iterrows():
            rows.append([names[s],f'{d.factor:.1f}',str(int(d.seed)),str(int(d.facilities)),f'{d.objective:.3f}',f'{100*d.demand_coverage:.3f}',f'{100*d.candidate_reachable_territorial_coverage:.3f}','Yes' if d.feasible else 'No'])
    parts.append(longtable('llrrrrrc',r'Scenario & Factor & Seed & Sites & Objective & Demand (\%) & Reach (\%) & Feasible',rows))
    parts.append(r'''\normalsize
The high-intensity statewide configuration for seed 1011 has 180 selected
sites and is infeasible because candidate-reachable area coverage is
79.9928\%, below the 80\% floor. Its territorial shortfall is 3.268 km$^2$.
The four feasible configurations at that setting select a mean of 182.50
sites (sample SD 3.70). The all-run mean is 182.00 (SD 3.39).

\section{Spatial distribution of the exact configurations}

The following summaries use the original feasible allocations associated with
the nominal four-site municipal and 13-site regional optima. The model's
location objective does not uniquely determine continuous allocations.
Consequently, these mapped allocations can differ from the maximum-demand
allocations used for fixed-site sensitivity assessment. Cells are assigned
to municipalities by centroid; the boundary fringe is reported separately.
''')
    f=pd.read_csv(R/'local_spatial_summary.csv');rows=[]
    for _,d in f.iterrows():
        city='Boundary fringe' if 'Boundary fringe' in d.municipality else d.municipality
        rows.append([names[d.scenario],city,str(d.cells),str(d.selected_sites),f'{d.allocated_demand:.2f}',f'{d.demand_coverage_pct:.2f}',f'{d.covered_grid_area_pct:.2f}'])
    parts.append(r'\small'+longtable('llrrrrr',r'Scenario & Municipality & Cells & Sites & Allocated & Demand (\%) & Area (\%)',rows))
    parts.append(r'''\normalsize
Area percentages in this table use each municipality's centroid-assigned
grid area, not the candidate-reachable denominator. Allocated quantities and
capacities below are synthetic planning units. Utilization is allocated
quantity divided by nominal capacity. Candidate indices are local to each
scenario and identify rows in the supplied candidate tables.
''')
    f=pd.read_csv(R/'exact_selected_site_details.csv');rows=[]
    for _,d in f.iterrows():
        rows.append([names[d.scenario],str(d.candidate_index),d.municipality,f'{d.longitude:.5f}',f'{d.latitude:.5f}',f'{d.capacity:.1f}',f'{d.allocated_demand:.1f}',f'{100*d.utilization:.1f}'])
    parts.append(r'\scriptsize'+longtable('lrlrrrrr',r'Scenario & Index & Municipality & Longitude & Latitude & Capacity & Allocated & Use (\%)',rows))
    parts.append(r'''\normalsize
\section{Variation among nominal relay selections}

The Jaccard index is the number of shared selected candidates divided by the
number selected in either configuration. These comparisons describe variation
among stochastic runs under identical nominal planning inputs. They do not
measure spatial displacement or establish an effect of demand composition.
''')
    f=pd.read_csv(R/'nominal_site_stability.csv');rows=[]
    for s in ['light','medium','heavy']:
        g=f[f.scenario==s].jaccard
        rows.append([names[s],str(len(g)),f'{g.mean():.3f}',f'{g.median():.3f}',f'{g.min():.3f}',f'{g.max():.3f}'])
    parts.append(longtable('lrrrrr',r'Scenario & Pairs & Mean & Median & Minimum & Maximum',rows))
    parts.append(r'''\section{Computational records and interpretation}

The computational supplement contains aligned cell polygons, populations,
areas, candidate coordinates and types, nominal preferences, and scenario
parameters. Continuous allocations and constraint residuals are retained for
all 147 sensitivity and spatial-assessment records, including infeasible
outcomes. The accompanying validation routine independently checks allocation
bounds, capacity, access, both coverage floors, objectives, and content digests.

Composition tests keep total demand, capacities, and preference scores fixed.
Subcell calculations divide each parent cell's population and area equally;
they do not incorporate independent 200 m population observations. Boundary
weighting preserves centroid locations and geometric access while scaling mass
and area by the intersection fraction. Statewide distance envelopes use
centroid-to-corner distances with a 1 m allowance. Their feasibility flags
concern access changes with nominal territorial coverage; the separate reach
ratio envelopes describe territorial uncertainty. Exact reoptimization uses
a 300 s limit and two threads, with remaining gaps retained. Statewide
composition reoptimization uses seed 42.

All graphical outputs are derived from the supplied data and numerical
algorithms. The manuscript describes assistance from OpenAI Codex (GPT-6;
OpenAI) in text preparation and implementation of analysis and visualization
routines. Scientific responsibility remains with the authors.
\end{document}''')
    (ROOT/'Supplementary_Information.tex').write_text('\n\n'.join(parts),encoding='utf-8')
if __name__=='__main__':main()
