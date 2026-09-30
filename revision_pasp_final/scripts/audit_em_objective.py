"""Replay retained independent counts/responses; no simulation or RNG draws."""
from common import *
from numba import njit, prange, set_num_threads
import csv, time
set_num_threads(4)
CAPS=np.array([200,500,1000,5000])
@njit(parallel=True)
def audit(counts,vals,ii,jj):
 B=len(counts); out=np.empty((4,B,54)); changes=np.empty((4,B)); mind=np.empty(B); decreases=np.zeros(B,np.int64); guards=np.zeros(B,np.int64)
 for b in prange(B):
  pi=np.full(54,1/54);eps=np.zeros(54)
  for k in range(len(ii)):eps[ii[k]]+=vals[b,k]
  prev=0.;minimum=1e300;ci=0;done=False
  for it in range(5001):
   den=np.zeros(54)
   for k in range(len(ii)):den[jj[k]]+=pi[ii[k]]*vals[b,k]
   eta=den.sum();ll=0.
   for j in range(54):
    if counts[b,j]>0:
     if den[j]<=0:guards[b]=1
     ll+=counts[b,j]*np.log(max(den[j]/eta,1e-300))
   if it>0:
    diff=ll-prev
    if diff<minimum:minimum=diff
    tolerance=1e-10+1e-12*max(abs(prev),abs(ll))
    if diff < -tolerance:decreases[b]+=1
   prev=ll
   if it==5000 or done:break
   weights=np.zeros(54)
   for k in range(len(ii)):
    d=den[jj[k]]
    if d<=0:guards[b]=1;d=1e-300
    weights[ii[k]]+=vals[b,k]*counts[b,jj[k]]/d
   new=pi*weights/eps;s=new.sum()
   if s<=0 or not np.isfinite(s):guards[b]=1;break
   new/=s;delta=0.
   for i in range(54):
    delta=max(delta,abs(new[i]-pi[i]))
   for i in range(54):pi[i]=new[i]
   converged=delta<1e-10
   if it+1==CAPS[ci] or converged:
    out[ci,b]=pi;changes[ci,b]=delta;ci+=1
    if converged:
     while ci<4:out[ci,b]=pi;changes[ci,b]=delta;ci+=1
     done=True
  mind[b]=minimum
 return out,changes,mind,decreases,guards

def main():
 dest=ROOT/'output/em_objective';dest.mkdir(parents=True,exist_ok=True);rows=[];hashes=[];start=time.time();maxres=0.
 ii,jj=np.nonzero(C[:,:-1])
 files=sorted((ROOT/'output/response_prior').glob('S?_??.npz'))+sorted((ROOT/'output/reconstruction_extended').glob('R?.npz'))
 for path in files:
  z=np.load(path);hashes.append({'path':str(path.relative_to(ROOT.parent)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
  batches=[('outer',z['outer_counts'][None,:],R[ii,jj][None,:],z['outer_estimates'][:,None,:]),('bootstrap',z['bootstrap_counts'],z['response_positive_values'],z['bootstrap_estimates'])] if 'outer_counts' in z else [('extended',z['counts'],np.tile(R[ii,jj],(len(z['counts']),1)),z['estimates'])]
  for kind,n,v,archived in batches:
   estimates,deltas,mins,dec,guards=audit(n,v,ii,jj);res=float(np.max(abs(estimates-archived)));maxres=max(maxres,res);assert res<1e-10,(path,kind,res)
   for b in range(len(n)):
    rows.append({'source':path.stem,'kind':kind,'fit':b,'min_delta_log_likelihood':float(mins[b]),'decreases_beyond_tolerance':int(dec[b]),'guard_triggered':int(guards[b]),**{f'final_max_delta_{cap}':float(deltas[c,b]) for c,cap in enumerate(CAPS)}})
  print(f'{path.name}: {len(rows)} fits; {time.time()-start:.1f}s',flush=True)
 with (dest/'per_fit_objective.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
 final=np.array([r['final_max_delta_5000'] for r in rows]);summary={'n_fits':len(rows),'n_decreases_beyond_tolerance':sum(r['decreases_beyond_tolerance'] for r in rows),'fits_with_decreases':sum(r['decreases_beyond_tolerance']>0 for r in rows),'minimum_delta_log_likelihood':min(r['min_delta_log_likelihood'] for r in rows),'max_replay_estimate_residual':maxres,'guard_fits':sum(r['guard_triggered'] for r in rows),'final_max_delta_5000_quantiles':{str(q):float(np.quantile(final,q)) for q in [0,.05,.25,.5,.75,.9,.95,.99,1]},'tolerance':'1e-10 + 1e-12 * max(abs(ll_previous),abs(ll_current)) in total log-likelihood units','objective':'sum_j n_j log[(sum_i t_i R_ij)/(sum_i t_i epsilon_i)]; multinomial count constant omitted','rng_draws':0,'inputs':hashes}
 summary['status']='PASS' if summary['n_decreases_beyond_tolerance']==0 and summary['guard_fits']==0 else 'FAIL'
 savejson(dest/'summary.json',summary)
 (ROOT/'reports/EM_OBJECTIVE_AUDIT.md').write_text('# EM OBJECTIVE AUDIT\n\nReplayed retained independent count vectors and response draws with identical uniform initialization and update. No astrophysical simulations or random draws. Every adjacent likelihood value from initialization through stopping/cap was checked; positive-count terms use normalized detected probabilities. Final arrays match the prior diagnostic arrays within 1e-10.\n\n```json\n'+json.dumps(summary,indent=2)+'\n```\n')
 print(json.dumps({k:v for k,v in summary.items() if k!='inputs'},indent=2))
 if summary['status']!='PASS':raise RuntimeError('STOP: likelihood monotonicity failed')
if __name__=='__main__':main()
