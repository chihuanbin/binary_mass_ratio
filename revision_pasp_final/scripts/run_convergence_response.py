from common import *
import csv,time
from numba import njit,prange,set_num_threads
set_num_threads(4)
@njit(parallel=True)
def sparse_fits(counts,vals,ii,jj,caps):
 B=len(counts);out=np.empty((len(caps),B,54));nit=np.zeros((len(caps),B),np.int64);deltas=np.empty((len(caps),B));guards=np.zeros(B,np.int64)
 for b in prange(B):
  pi=np.full(54,1/54);eps=np.zeros(54)
  for k in range(len(ii)):eps[ii[k]]+=vals[b,k]
  capindex=0;converged=False;delta=1.;it=0
  for it in range(1,caps[-1]+1):
   den=np.zeros(54)
   for k in range(len(ii)):den[jj[k]]+=pi[ii[k]]*vals[b,k]
   weights=np.zeros(54)
   for k in range(len(ii)):
    d=den[jj[k]]
    if d<=0:guards[b]=1;d=1e-300
    weights[ii[k]]+=vals[b,k]*counts[b,jj[k]]/d
   new=pi*weights/eps;s=new.sum()
   if s<=0 or not np.isfinite(s):guards[b]=1;break
   new/=s
   delta=0.0
   for i in range(54):
    diff=abs(new[i]-pi[i])
    if diff>delta:delta=diff
   for i in range(54):pi[i]=new[i]
   if delta<1e-10:converged=True
   if it==caps[capindex] or converged:
    out[capindex,b]=pi;nit[capindex,b]=it;deltas[capindex,b]=delta;capindex+=1
    if converged:
     while capindex<len(caps):out[capindex,b]=pi;nit[capindex,b]=it;deltas[capindex,b]=delta;capindex+=1
     break
   if capindex==len(caps):break
 return out,nit,deltas,guards

def main():
 dest=ROOT/'output/response_prior';dest.mkdir(parents=True,exist_ok=True);conv=ROOT/'output/convergence';conv.mkdir(parents=True,exist_ok=True)
 ii,jj=np.nonzero(C[:,:-1]);alpha=C[ii,jj]+1.;caps=np.array([200,500,1000,5000]);u=truths()[0];s1=u.copy();s1[30:32]*=.75;s1/=s1.sum();rng=np.random.default_rng(2026093006);allrows=[];summary={};long_summary={};started=time.time()
 # Verify optimized update against archived implementation on independent counts.
 import sys
 sys.path.insert(0,str(OLD/'src/legacy_scripts'));from unfolding import dagostini_observable_only
 trial=rng.multinomial(20000,project(u),size=10);vals=np.tile(R[ii,jj],(10,1));o,ni,de,gu=sparse_fits(trial,vals,ii,jj,np.array([200]))
 residual=max(np.max(abs(o[0,b]-dagostini_observable_only(R,trial[b],n_iter=200))) for b in range(10));assert residual<1e-10
 for pop,t in [('S0',u),('S1',s1)]:
  indicators=[];widths=[];longind=[];longwidth=[];point=[];pointlong=[]
  for s in range(20):
   n=rng.multinomial(20000,project(t));outer,its,ds,gs=sparse_fits(n[None,:],R[ii,jj][None,:],ii,jj,caps);fit=outer[0,0];point.append(fit);pointlong.append(outer[-1,0]);pfit=project(fit);pfit=np.maximum(pfit,1e-12);pfit/=pfit.sum();nb=rng.multinomial(20000,pfit,size=1000)
   values=np.empty((1000,len(ii)))
   for i in range(54):
    mask=ii==i;values[:,mask]=rng.dirichlet(alpha[mask],size=1000)
   boot,bit,bdel,bguard=sparse_fits(nb,values,ii,jj,caps)
   for label,outputs,niters,deltas,guard in [('outer',outer,its,ds,gs),('bootstrap',boot,bit,bdel,bguard)]:
    for b in range(outputs.shape[1]):
     for c,cap in enumerate(caps):
      v=outputs[c,b];allrows.append(dict(population=pop,realization=s,kind=label,fit=b,cap=int(cap),converged=bool(deltas[c,b]<1e-10 and guard[b]==0),n_iter=int(niters[c,b]),final_max_delta=float(deltas[c,b]),hit_iteration_cap=bool(niters[c,b]==cap and deltas[c,b]>=1e-10),numerical_guard_triggered=bool(guard[b]),finite_output=bool(np.isfinite(v).all()),sum_output=float(v.sum()),min_output=float(v.min()),max_output=float(v.max()),delta_from_200=float(np.max(abs(v-outputs[0,b])))))
   lo,hi=np.quantile(boot[0],[.16,.84],axis=0);ll,hl=np.quantile(boot[-1],[.16,.84],axis=0);indicators.append((lo<=t)&(t<=hi));widths.append(hi-lo);longind.append((ll<=t)&(t<=hl));longwidth.append(hl-ll)
   np.savez_compressed(dest/f'{pop}_{s:02d}.npz',truth=t,outer_counts=n,outer_estimates=outer[:,0],bootstrap_counts=nb,response_positive_values=values,response_i=ii,response_j=jj,bootstrap_estimates=boot,lo=lo,hi=hi,lo_long=ll,hi_long=hl)
   print(f'{pop} {s+1}/20: elapsed {time.time()-started:.1f}s',flush=True)
  def summarize(ind,w):
   ind=np.array(ind);per=ind.mean(1);crng=np.random.default_rng(2026093007);ci=np.quantile(crng.choice(per,(20000,20)).mean(1),[.025,.975]);return {'k':int(ind.sum()),'n':1080,'coverage':float(ind.mean()),'cluster95':ci.tolist(),'all54':int(ind.all(1).sum()),'n_realizations':20,'mean_width':float(np.mean(w)),'per_realization_coverage':per.tolist()}
  summary[pop]=summarize(indicators,widths);long_summary[pop]=summarize(longind,longwidth)
  point=np.array(point);pl=np.array(pointlong);summary[pop]['mean_E200']=float(np.max(abs(point-t),1).mean());summary[pop]['mean_Elong']=float(np.max(abs(pl-t),1).mean());summary[pop]['mean_q_error200']=float(np.mean((point-t)@Q));summary[pop]['mean_q_errorlong']=float(np.mean((pl-t)@Q))
 with (conv/'per_fit_diagnostics.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]));w.writeheader();w.writerows(allrows)
 savejson(dest/'results.json',{'structural_zero_200':summary,'structural_zero_5000':long_summary,'historical_comparisons':json.loads((L/'exp5_operator_uncertainty.json').read_text()),'independent_new_experiment':True,'seed':2026093006,'cluster_seed':2026093007,'bootstrap_replicates':1000,'caps':caps.tolist(),'optimized_vs_archived_implementation_max_delta':float(residual),'historical_fit_authentication':'BLOCKED: original count vectors and per-state estimates absent'})
 csum={}
 for cap in caps:
  rows=[r for r in allrows if r['cap']==cap];csum[str(cap)]={'n_fits':len(rows),'converged':sum(r['converged'] for r in rows),'hit_cap':sum(r['hit_iteration_cap'] for r in rows),'guard':sum(r['numerical_guard_triggered'] for r in rows),'max_delta_from_200':max(r['delta_from_200'] for r in rows),'mean_delta_from_200':float(np.mean([r['delta_from_200'] for r in rows]))}
 savejson(conv/'summary.json',csum)
 (ROOT/'reports/EM_CONVERGENCE_AUDIT.md').write_text('# EM CONVERGENCE AUDIT\n\nHistorical exact-fit authentication BLOCKED: original counts, estimates and traces absent. Independent new outer and structural-zero bootstrap fits use frozen response, uniform initialization, original guards and tolerance 1e-10. Algebraic optimized update verified against archived implementation on 10 independent counts (not historical fits). Every new fit is audited at 200/500/1000/5000 caps; cap hitting is not convergence.\n\n```json\n'+json.dumps(csum,indent=2)+'\n```\n')
 (ROOT/'reports/RESPONSE_PRIOR_SENSITIVITY.md').write_text('# RESPONSE PRIOR SENSITIVITY\n\nIndependent new experiment: 20 realizations/truth, 1000 bootstraps each, seed 2026093006. Positive-count cells use Dirichlet(C+1), empirical zeros (including NOT_OBSERVED) remain exactly zero. Fitted nominal distribution generates bootstrap counts; perturbed response used only in inverse. Comparison to fixed and all-cell procedures is historical, not paired. 200-iteration results are primary to match historical implementation; long-cap sensitivity is separately reported.\n\n```json\n'+json.dumps({'200':summary,'5000':long_summary},indent=2)+'\n```\n')
 print(json.dumps(csum,indent=2))
if __name__=='__main__':main()
