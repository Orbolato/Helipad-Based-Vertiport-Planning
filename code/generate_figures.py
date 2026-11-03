"""Generate vector figures and spatial summaries from model inputs and allocations."""
from analysis_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
import geopandas as gpd
import pickle,itertools
FIG=ROOT/'figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,'axes.labelsize':12,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
COLORS={'VNS':'#0072B2','SMA+VNS':'#D55E00'}
NAMES={'VNS':'VNS','SMA+VNS':'Relay'}
CITY={'3550308':'Sao Paulo','3509502':'Campinas','3518800':'Guarulhos'}

def save(fig,name):
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',metadata={'Author':'','Creator':'Matplotlib'})
    fig.savefig(ROOT/'figures'/(name+'.png'),dpi=150,bbox_inches='tight')
    plt.close(fig)

def aggregate(runs,col):
    ticks=np.unique(np.concatenate([r['hist'].fitness_calls.to_numpy() for r in runs.values()]))
    vals=[]
    for r in runs.values():
        h=r['hist'].drop_duplicates('fitness_calls',keep='last').sort_values('fitness_calls')
        idx=np.searchsorted(h.fitness_calls,ticks,side='right')-1
        v=h[col].to_numpy()[np.clip(idx,0,len(h)-1)].astype(float)
        v[ticks<h.fitness_calls.iloc[0]]=np.nan
        vals.append(v)
    arr=np.array(vals)
    return ticks,np.nanmedian(arr,axis=0),np.nanmin(arr,axis=0),np.nanmax(arr,axis=0)

def statewide():
    data=pickle.loads((ROOT/'data/statewide_traces.pkl').read_bytes())
    hel=pd.read_csv(ROOT/'data/heavy_candidates.csv');state=gpd.read_file(ROOT/'data/heavy_boundaries.gpkg').dissolve()
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.8))
    fig.subplots_adjust(left=.08,right=.98,bottom=.15,top=.82,wspace=.23)
    for ax,method in zip(axes,['VNS','SMA+VNS']):
        feasible=[(seed,r) for seed,r in data['meta'][method].items() if r['certificate_summary']['feasible']]
        seed,r=min(feasible,key=lambda t:t[1]['certificate_summary']['objective'])
        state.plot(ax=ax,color='#f6f6f6',edgecolor='#444444',linewidth=.7)
        ax.scatter(hel.Longitude,hel.Latitude,s=8,color='#aaaaaa',alpha=.55)
        sel=hel.loc[np.asarray(r['x_sol'],dtype=bool)]
        ax.scatter(sel.Longitude,sel.Latitude,s=18,c=COLORS[method],edgecolor='white',linewidth=.25)
        ax.set_xlabel('Longitude (degrees)')
        # Horizontal labeling avoids orientation-dependent rendering.
        ax.set_ylabel('Latitude\n(degrees)' if ax is axes[0] else '',rotation=0,labelpad=30,va='center')
        ax.set_title(f'{NAMES[method]}: {len(sel)} sites (seed {seed})')
        ax.tick_params(labelsize=9)
    fig.legend(handles=[Line2D([],[],marker='o',linestyle='',color='#aaa',label='All candidates'),Line2D([],[],marker='o',linestyle='',color=COLORS['VNS'],label='Selected VNS'),Line2D([],[],marker='o',linestyle='',color=COLORS['SMA+VNS'],label='Selected relay')],loc='upper center',ncol=3,frameon=False,bbox_to_anchor=(.53,1.))
    save(fig,'statewide_selected_sites')
    fig,axes=plt.subplots(1,3,figsize=(7.2,3.7))
    fig.subplots_adjust(left=.07,right=.99,bottom=.18,top=.76,wspace=.32)
    for ax,(col,title,floor) in zip(axes,[('cov_dem_rel','(a) Allocated demand',30),('cov_ter_svc','(b) Candidate-reachable\narea',80),('cov_ter_region','(c) Total state\ngrid area',None)]):
        for method in ['VNS','SMA+VNS']:
            t,med,low,high=aggregate(data['meta'][method],col)
            ax.fill_between(t,100*low,100*high,step='post',color=COLORS[method],alpha=.16)
            ax.step(t,100*med,where='post',color=COLORS[method],label=NAMES[method],linewidth=1.7)
        if floor is not None:ax.axhline(floor,color='#444',linestyle=':',linewidth=1)
        ax.set_title(title);ax.set_xlabel('Evaluations');ax.set_ylabel('Coverage (%)' if ax is axes[0] else '')
        ax.ticklabel_format(axis='x',style='sci',scilimits=(3,3));ax.grid(alpha=.15)
    fig.legend(handles=[Line2D([],[],color=COLORS[m],label=NAMES[m]) for m in COLORS]+[Patch(facecolor='#888',alpha=.2,label='Five-seed range'),Line2D([],[],color='#444',linestyle=':',label='Coverage floor')],loc='upper center',bbox_to_anchor=(.52,1.03),ncol=2,frameon=False)
    save(fig,'statewide_coverage_traces')
    fig,ax=plt.subplots(figsize=(7.2,3.7));fig.subplots_adjust(left=.1,right=.98,bottom=.17,top=.79)
    for method in COLORS:
        t,med,low,high=aggregate(data['meta'][method],'best_fit')
        ax.fill_between(t,low,high,step='post',color=COLORS[method],alpha=.16)
        ax.step(t,med,where='post',color=COLORS[method],label=NAMES[method])
    ax.set_yscale('log');ax.set_ylabel('Penalized fitness');ax.set_xlabel('Objective-function evaluations');ax.grid(alpha=.15)
    fig.legend(handles=[Line2D([],[],color=COLORS[m],label=NAMES[m]) for m in COLORS]+[Patch(facecolor='#888',alpha=.2,label='Five-seed range')],loc='upper center',ncol=3,frameon=False)
    save(fig,'statewide_objective_convergence')

def local_maps():
    fig,axes=plt.subplots(2,2,figsize=(7.2,9.1));fig.subplots_adjust(left=.08,right=.96,bottom=.28,top=.94,hspace=.35,wspace=.38)
    summary=[];site_rows=[];allocation_rows=[]
    panels=[('light','3550308','(a) Municipal\nSao Paulo'),('medium','3550308','(b) Regional\nSao Paulo'),('medium','3509502','(c) Regional\nCampinas'),('medium','3518800','(d) Regional\nGuarulhos')]
    instances={}
    for scenario in ['light','medium']:
        inst=load_instance(scenario);e=engine(inst)
        z=np.load(ROOT/'data/original_solutions'/scenario/'exact_optimum.npz')
        x=z['x'];qi=z['allocation_i'];qj=z['allocation_j'];q=z['allocation_q']
        cert=opt.verify_allocation_certificate(x_sol=x,allocation_i=qi,allocation_j=qj,allocation_q=q,demand=e.D,capacities=e.C,area_km2=e.A,reach_facility=inst['fac'],reach_aircraft=inst['air'],dist_vv=inst['vv'],eta_demanda=e.eta_dem,eta_territorio=e.eta_ter,min_sep_km=.3,pref_score_j=e.pref_score,w_open=1,w_pref=.5,w_sep_penalty=500)
        assert cert['summary']['feasible'],scenario
        served=np.bincount(qi,weights=q,minlength=len(e.D));load=np.bincount(qj,weights=q,minlength=len(e.C))
        coverage=inst['air'][:,np.flatnonzero(x)].any(axis=1)
        grid=gpd.read_file(ROOT/'data'/f'{scenario}_grid.gpkg').to_crs(4326)
        mun=gpd.read_file(ROOT/'data'/f'{scenario}_boundaries.gpkg').to_crs(4326)
        cent=gpd.GeoSeries(gpd.points_from_xy(inst['dem'].Longitude,inst['dem'].Latitude),crs=4326)
        sites=gpd.GeoSeries(gpd.points_from_xy(inst['hel'].Longitude,inst['hel'].Latitude),crs=4326)
        for j in np.flatnonzero(x):
            code=next((str(row.CD_MUN) for _,row in mun.iterrows() if sites.iloc[j].within(row.geometry)),'outside')
            site_rows.append(dict(scenario=scenario,candidate_index=int(j),municipality=CITY.get(code,code),longitude=float(inst['hel'].Longitude.iloc[j]),latitude=float(inst['hel'].Latitude.iloc[j]),facility_type=inst['hel'].Tipo.iloc[j],capacity=float(e.C[j]),allocated_demand=float(load[j]),utilization=float(load[j]/e.C[j])))
        for i,j,v in zip(qi,qj,q):allocation_rows.append(dict(scenario=scenario,cell_index=int(i),candidate_index=int(j),allocated_demand=float(v)))
        mask_all=np.zeros(len(e.D),bool)
        for _,row in mun.iterrows():
            mask=cent.within(row.geometry).to_numpy();mask_all|=mask
            select=sites.within(row.geometry).to_numpy() & x.astype(bool)
            summary.append(dict(scenario=scenario,municipality=CITY.get(str(row.CD_MUN),str(row.CD_MUN)),cells=int(mask.sum()),population=float(inst['dem'].Populacao[mask].sum()),selected_sites=int(select.sum()),allocated_demand=float(served[mask].sum()),demand_coverage_pct=float(100*served[mask].sum()/e.D[mask].sum()),covered_grid_area_pct=float(100*e.A[mask & coverage].sum()/e.A[mask].sum())))
        mask=~mask_all
        summary.append(dict(scenario=scenario,municipality='Boundary fringe (centroid outside)',cells=int(mask.sum()),population=float(inst['dem'].Populacao[mask].sum()),selected_sites=0,allocated_demand=float(served[mask].sum()),demand_coverage_pct=float(100*served[mask].sum()/e.D[mask].sum()),covered_grid_area_pct=float(100*e.A[mask & coverage].sum()/e.A[mask].sum())))
        instances[scenario]=(inst,grid,mun,served,e.D,coverage,x)
    for ax,(scenario,code,title) in zip(axes.ravel(),panels):
        inst,grid,mun,served,D,cov,x=instances[scenario]
        geom=mun.loc[mun.CD_MUN.astype(str)==code]
        box=geom.total_bounds;padx=(box[2]-box[0])*.06;pady=(box[3]-box[1])*.06
        visible=grid.intersects(geom.union_all()).to_numpy()
        grid.loc[visible].plot(ax=ax,facecolor='#ededed',edgecolor='none')
        grid.loc[visible & cov].plot(ax=ax,facecolor='#f3d78b',edgecolor='none')
        frac=np.divide(served,D,out=np.zeros_like(D),where=D>0)
        served_cells=visible & (served>1e-7)
        grid.loc[served_cells].plot(ax=ax,column=frac[served_cells],cmap='Blues',vmin=0,vmax=1,edgecolor='none')
        geom.boundary.plot(ax=ax,color='#222',linewidth=.8)
        ax.scatter(inst['hel'].Longitude,inst['hel'].Latitude,s=6,color='#666',alpha=.45)
        selected=np.array([j for j in np.flatnonzero(x) if geom.union_all().contains(__import__('shapely').geometry.Point(inst['hel'].Longitude.iloc[j],inst['hel'].Latitude.iloc[j]))],dtype=int)
        ax.scatter(inst['hel'].Longitude.iloc[selected],inst['hel'].Latitude.iloc[selected],s=62,marker='*',color='#c23b22',edgecolor='white',linewidth=.5,zorder=5)
        for j in selected:
            px=inst['hel'].Longitude.iloc[j];py=inst['hel'].Latitude.iloc[j]
            # Site indices and coordinates are tabulated separately.
        ax.set_xlim(box[0]-padx,box[2]+padx);ax.set_ylim(box[1]-pady,box[3]+pady)
        ax.set_title(title);ax.set_xlabel('Longitude (degrees)');ax.set_ylabel('Latitude (degrees)');ax.tick_params(labelsize=10); ax.xaxis.set_major_locator(__import__('matplotlib').ticker.MaxNLocator(4)); ax.yaxis.set_major_locator(__import__('matplotlib').ticker.MaxNLocator(5))
    fig.legend(handles=[Patch(facecolor='#ededed',label='Outside selected 10 km reach'),Patch(facecolor='#f3d78b',label='Within reach, no allocation'),Patch(facecolor='#2879b9',label='Allocated cells'),Line2D([],[],marker='*',markersize=10,color='#c23b22',linestyle='',label='Selected site')],loc='lower center',bbox_to_anchor=(.52,.0),ncol=1,frameon=False)
    cax=fig.add_axes([.23,.175,.54,.013]);bar=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='Blues'),cax=cax,orientation='horizontal');bar.set_label('Fraction of cell planning demand allocated',fontsize=10)
    save(fig,'municipal_regional_spatial_allocation')
    pd.DataFrame(summary).to_csv(ROOT/'results/local_spatial_summary.csv',index=False)
    pd.DataFrame(site_rows).to_csv(ROOT/'results/exact_selected_site_details.csv',index=False)
    pd.DataFrame(allocation_rows).to_csv(ROOT/'results/exact_allocation_arcs.csv',index=False)
    print(pd.DataFrame(summary).to_string(index=False))

def stability():
    rows=[]
    for scenario in ['light','medium','heavy']:
        xs={seed:np.load(ROOT/'results'/f'{scenario}_intensity_1.0_seed_{seed}.npz')['x'].astype(bool) for seed in SEEDS}
        for a,b in itertools.combinations(SEEDS,2):rows.append(dict(scenario=scenario,seed_a=a,seed_b=b,jaccard=float(np.count_nonzero(xs[a]&xs[b])/np.count_nonzero(xs[a]|xs[b]))))
    pd.DataFrame(rows).to_csv(ROOT/'results/nominal_site_stability.csv',index=False)

if __name__=='__main__':
    statewide();local_maps();stability()
