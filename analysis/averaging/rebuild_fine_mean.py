"""Rebuild equal-weight fine-scale means from all 100 archived realizations."""
from pathlib import Path
import numpy as np,json,hashlib
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True,help='Experiment output directory')
OUT=parser.parse_args().output.resolve()
m=json.loads((OUT/'manifest.json').read_text());arch=Path(m['archive'])
t=np.linspace(.005,19.995,20001)
means={};sources=[]
for typ,ext in [('pdf','_derivative'),('cdf','')]:
 total=np.zeros((len(t),3));left=0.;right=float('inf')
 for i in range(1,101):
  p=arch/f'fine_r{i:04d}'/f'r{i:04d}_BTC_FineScaled{ext}.csv';a=np.loadtxt(p,delimiter=',',skiprows=1)
  assert a.shape[1]==6 and np.isfinite(a).all()
  left=max(left,float(a[0,::2].max()));right=min(right,float(a[-1,::2].min()))
  assert left<=t[0] and right>=t[-1]
  for j in range(3):
   assert np.all(np.diff(a[:,2*j])>0)
   total[:,j]+=np.interp(t,a[:,2*j],a[:,2*j+1])/100
  sources.append({'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'rows':len(a)})
 means[typ]=total
 np.savetxt(OUT/f'fine_mean_{typ}.csv',np.c_[t,total],delimiter=',',header='time,x0.5,x1.5,x2.5',comments='')
 print(typ,'complete',left,right,flush=True)
b=np.loadtxt(arch/'BTC_mean.csv',delimiter=',',skiprows=1)
err=[]
for j in range(3):
 y=np.interp(t,b[:,2*j],b[:,2*j+1]);z=means['pdf'][:,j]
 err.append(float(100*np.sqrt(np.trapz((y-z)**2,t)/np.trapz(z*z,t))))
result={'sources':sources,'shared_time':[float(t[0]),float(t[-1])],'mean':'equal weight arithmetic mean of interpolated per-realization curves','archived_BTC_mean_relative_L2_percent_vs_rebuilt':err}
(OUT/'fine_mean_audit.json').write_text(json.dumps(result,indent=2)+'\n');print('Archived mean vs rebuilt L2 %',err)
