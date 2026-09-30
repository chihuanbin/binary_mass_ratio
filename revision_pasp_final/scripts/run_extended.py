from common import *
from run_convergence_response import sparse_fits
import csv

def main():
 dest=ROOT/'output/reconstruction_extended';dest.mkdir(exist_ok=True,parents=True);q=np.arange(.15,1,.1)
 rising=q.copy();excess=1+2*np.exp(-.5*((q-.9)/.1)**2);md=np.array([1.1-q if i<3 else q for i in range(6)]);md/=md.sum(1)[:,None];design={'R1':np.tile(rising/rising.sum(),(6,1))/6,'R2':np.tile(excess/excess.sum(),(6,1))/6,'R3':md/6}
 savejson(ROOT/'config/extended_truth_design.json',{'frozen_before_sampling':True,'R1':'within-block weights proportional to q-bin centers','R2':'within-block weights 1+2 exp[-0.5((q_center-0.9)/0.1)^2]','R3':'first three mass blocks weights 1.1-q_center; last three weights q_center; each mass marginal 1/6','n_realizations':500,'N':20000,'seed':2026093008,'caps':[200,500,1000,5000]})
 rng=np.random.default_rng(2026093008);ii,jj=np.nonzero(C[:,:-1]);summary={};rows=[];diagnostics=[]
 for name,mat in design.items():
  t=mat.ravel();assert abs(t.sum()-1)<1e-14;np.save(dest/f'{name}_truth.npy',t);n=rng.multinomial(20000,project(t),500);v,ni,de,gu=sparse_fits(n,np.tile(R[ii,jj],(500,1)),ii,jj,np.array([200,500,1000,5000]));np.savez_compressed(dest/f'{name}.npz',truth=t,counts=n,estimates=v,n_iter=ni,final_delta=de,guards=gu)
  def metrics(h):
   h=np.atleast_2d(h);d=h-t;m=(h+t)/2;js=.5*(np.sum(np.where(h>0,h*np.log(np.maximum(h,1e-300)/m),0),1)+np.sum(t*np.log(t/m),1));cond=np.sum(h.reshape(-1,6,9)*q,2)/np.sum(h.reshape(-1,6,9),2);truecond=np.sum(mat*q,1)/mat.sum(1)
   return {'E':np.max(abs(d),1),'L1':np.sum(abs(d),1),'RMSE':np.sqrt(np.mean(d*d,1)),'JS':js,'mean_q_error':d@Q,**{f'conditional_q_error_{i}':cond[:,i]-truecond[i] for i in range(6)}}
  summary[name]={}
  for label,h in [('EM200',v[0]),('EM5000',v[-1]),('uniform',np.full((500,54),1/54))]:
   met=metrics(h);summary[name][label]={k:float(x.mean()) for k,x in met.items()}
   for b in range(500):rows.append({'truth':name,'estimator':label,'realization':b,**{k:float(x[b]) for k,x in met.items()}})
  summary[name]['converged200']=int((de[0]<1e-10).sum());summary[name]['converged5000']=int((de[-1]<1e-10).sum())
  for c,cap in enumerate([200,500,1000,5000]):
   for b,vec in enumerate(v[c]):diagnostics.append(dict(truth=name,realization=b,cap=cap,converged=bool(de[c,b]<1e-10 and not gu[b]),n_iter=int(ni[c,b]),final_max_delta=float(de[c,b]),hit_iteration_cap=bool(ni[c,b]==cap and de[c,b]>=1e-10),numerical_guard_triggered=bool(gu[b]),finite_output=bool(np.isfinite(vec).all()),sum_output=float(vec.sum()),min_output=float(vec.min()),max_output=float(vec.max())))
 with (dest/'per_realization_metrics.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 with (dest/'per_fit_diagnostics.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=diagnostics[0].keys());w.writeheader();w.writerows(diagnostics)
 savejson(dest/'summary.json',summary)
 (ROOT/'reports/EXTENDED_RECONSTRUCTION_BENCHMARK.md').write_text('# EXTENDED RECONSTRUCTION BENCHMARK\n\nIndependent extension, analytically declared before sampling, 500 realizations/truth. Same response and initialization. No interpretation as observed twin population or physical mass dependence. All errors including conditional q means retained.\n\n```json\n'+json.dumps(summary,indent=2)+'\n```\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
