"""Three-way comparison, numerical sensitivity, and reproducible figures."""
from pathlib import Path
import numpy as np,json,csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True,help='Experiment output directory')
OUT=parser.parse_args().output.resolve()
m=json.loads((OUT/'manifest.json').read_text());t=np.loadtxt(OUT/'fine_mean_pdf.csv',delimiter=',',skiprows=1)[:,0]
def load(name,kind):
 a=np.loadtxt(OUT/'runs'/f'{name}_{kind}.csv',delimiter=',',skiprows=1)
 assert a[0,::2].max()<=t[0] and a[-1,::2].min()>=t[-1]
 return np.stack([np.interp(t,a[:,2*j],a[:,2*j+1]) for j in range(3)],axis=1)
def l2(a,b):return 100*np.sqrt(np.trapz((a-b)**2,t,axis=0)/np.trapz(b*b,t,axis=0))
results={};rows=[];refinements=[]
for kind in ['pdf','cdf']:
 individual=np.stack([load(f'r{i:04d}',kind) for i in range(1,101)])
 data={'Fine-scale mean':np.loadtxt(OUT/f'fine_mean_{kind}.csv',delimiter=',',skiprows=1)[:,1:],
       'Upscale each, then average':individual.mean(0),
       'Average inputs, then upscale':load('mean',kind),
       'Saved mean inputs, rerun':load('archived_mean',kind)}
 results[kind]=data
 np.savetxt(OUT/f'comparison_{kind}.csv',np.column_stack([t,*data.values()]),delimiter=',',header='time,'+','.join(f'{k}_x{x}' for k in data for x in [.5,1.5,2.5]),comments='')
 for name,a in data.items():
  errs=l2(a,data['Fine-scale mean'])
  for j in range(3):
   row={'curve':name,'kind':kind,'x':.5+j,'relative_L2_percent':float(errs[j]),'value_at_end':float(a[-1,j])}
   if kind=='pdf':
    mass=np.trapz(a[:,j],t);mu=np.trapz(t*a[:,j],t)/mass;var=np.trapz((t-mu)**2*a[:,j],t)/mass
    row.update(peak_time=float(t[np.argmax(a[:,j])]),peak_height=float(a[:,j].max()),area=float(mass),window_conditional_mean=float(mu),window_conditional_std=float(np.sqrt(max(var,0))),tail_area_after_5=float(np.trapz(a[t>=5,j],t[t>=5])))
   else:
    for q in [.1,.5,.9]:row[f't{int(q*100)}']=float(np.interp(q,a[:,j],t)) if a[0,j]<=q<=a[-1,j] else None
   rows.append(row)
 for name in ['mean','archived_mean','r0001','r0050','r0100']:
  refined=load(name+'_refined',kind);coarse=load(name,kind)
  for j,v in enumerate(l2(coarse,refined)):refinements.append({'case':name,'kind':kind,'x':.5+j,'coarse_vs_refined_L2_percent':float(v)})
 if kind=='pdf':
  for tail in [False,True]:
   fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
   colors=['black','#2474b5','#d85c20','#9d77a9']
   for j,ax in enumerate(axs):
    for (name,a),col in zip(data.items(),colors):
     if name=='Saved mean inputs, rerun':continue
     ax.plot(t,a[:,j],label=name,color=col,lw=1.8)
    if tail:ax.set_yscale('log');ax.set_ylim(1e-5,10);ax.set_xlim(.01,20)
    else:ax.set_xlim(0,3+3*j)
    ax.set(xlabel='Time (simulation units)',title=f'x = {j+.5}');ax.grid(alpha=.2)
   axs[0].set_ylabel('Derivative of concentration BTC');axs[0].legend(fontsize=8)
   fig.suptitle('100 paired realizations · D = 0')
   fig.savefig(OUT/('comparison_tail.png' if tail else 'comparison_peaks.png'),dpi=180);plt.close(fig)
 else:
  fig,axs=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
  for j,ax in enumerate(axs):
   for name,a in data.items():
    if name=='Saved mean inputs, rerun':continue
    ax.plot(t,a[:,j],label=name)
   ax.set(xlabel='Time (simulation units)',title=f'x = {j+.5}',xlim=(0,10),ylim=(0,1.05));ax.grid(alpha=.2)
  axs[0].set_ylabel('Concentration / inlet concentration');axs[0].legend(fontsize=8)
  fig.savefig(OUT/'comparison_cdf.png',dpi=180);plt.close(fig)
resolution=[]
for kind in ['pdf','cdf']:
 for name in ['mean_space','mean_rank']:
  for j,v in enumerate(l2(load('mean',kind),load(name,kind))):resolution.append({'case':name,'kind':kind,'x':.5+j,'baseline_vs_finer_L2_percent':float(v)})
for name,rs in [('metrics.csv',rows),('time_step_check.csv',refinements),('grid_check.csv',resolution)]:
 fields=list(dict.fromkeys(k for r in rs for k in r))
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rs)
# Reproduction comparison against original archived upscaled PDF.
a=np.loadtxt(Path(m['archive'])/'upscaled_mean/upscaled_BTC_Upscaled_derivative.csv',delimiter=',',skiprows=1)
b=np.stack([np.interp(t,a[:,2*j],a[:,2*j+1]) for j in range(3)],axis=1)
repro=l2(results['pdf']['Saved mean inputs, rerun'],b)
fine=results['pdf']['Fine-scale mean'];e=l2(results['pdf']['Upscale each, then average'],fine);a=l2(results['pdf']['Average inputs, then upscale'],fine)
summary={'time_window':[float(t[0]),float(t[-1])],'individual_then_mean_L2_percent':e.tolist(),'mean_then_solve_L2_percent':a.tolist(),'relative_error_reduction_percent':(100*(1-e/a)).tolist(),'rerun_saved_mean_vs_archived_upscaled_L2_percent':repro.tolist(),'refinement':refinements,'grid_sensitivity':resolution,'fine_mean_audit':json.loads((OUT/'fine_mean_audit.json').read_text())['archived_BTC_mean_relative_L2_percent_vs_rebuilt']}
(OUT/'results.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='refinement'},indent=2))

# Human-readable results are calculated from this run, never embedded in source.
report=['# Averaging-order comparison', '', '| x | Average inputs first: L2 error | Upscale first: L2 error |', '|---|---:|---:|']
for j in range(3):
 report.append(f'| {j+.5} | {a[j]:.3f}% | {e[j]:.3f}% |')
report += ['', 'Errors compare concentration derivatives with the reconstructed fine mean on t=0.005–19.995, without normalization or extrapolation.', '', '![Peaks](comparison_peaks.png)', '', '![Tail](comparison_tail.png)', '', 'See metrics.csv, time_step_check.csv, grid_check.csv, and validation.json for diagnostics. These are sensitivity checks, not a full convergence proof.']
(OUT/'results.md').write_text('\n'.join(report)+'\n')
