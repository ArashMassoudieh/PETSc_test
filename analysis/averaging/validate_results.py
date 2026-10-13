from pathlib import Path
import numpy as np,json,hashlib
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True,help='Experiment output directory')
p=parser.parse_args().output.resolve()
m=json.loads((p/'manifest.json').read_text());names=[x['name'] for x in m['cases']]+[n+'_refined' for n in ['mean','archived_mean','r0001','r0050','r0100']]+['mean_space','mean_rank']
if (p/'runs/r0003_tight_done.json').exists():names.append('r0003_tight')
checks=[]
for n in names:
 done=json.loads((p/'runs'/f'{n}_done.json').read_text());assert done['returncode']==0
 c=np.loadtxt(p/'runs'/f'{n}_cdf.csv',delimiter=',',skiprows=1);d=np.loadtxt(p/'runs'/f'{n}_pdf.csv',delimiter=',',skiprows=1)
 assert np.isfinite(c).all() and np.isfinite(d).all()
 assert c[:,1::2].min()>-1e-5
 assert np.diff(c[:,1::2],axis=0).min()>-1e-5 and d[:,1::2].min()>-1e-5
 assert c[-1,0]>=19.995 and d[-1,0]>=19.995
 checks.append({'case':n,'cdf_min':float(c[:,1::2].min()),'cdf_max':float(c[:,1::2].max()),'pdf_min':float(d[:,1::2].min()),'cdf_end':c[-1,1::2].tolist(),'strict_upper_bound_pass':bool(c[:,1::2].max()<1.00001)})
warnings=[c['case'] for c in checks if not c['strict_upper_bound_pass']]
(p/'validation.json').write_text(json.dumps({'n_runs':len(checks),'status':'completed_with_upper_bound_warning' if warnings else 'passed','strict_upper_bound_failures':warnings,'max_upper_overshoot_percent':100*max(0,max(c['cdf_max'] for c in checks)-1),'result':'All completed and finite with nondecreasing CDFs, nonnegative derivatives within 1e-5, and full interval coverage. The strict unit-inlet upper-bound check is recorded separately; overshoots are retained without clipping.','checks':checks},indent=2)+'\n')
print('Checked',len(checks),'completed runs; strict upper-bound warnings:',len(warnings))
