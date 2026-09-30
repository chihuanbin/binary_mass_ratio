from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT.parent/'pasp_cmd_recoverability_revision';L=OLD/'data/derived/legacy_results'
R=np.load(L/'calibration_200k.npz',allow_pickle=True)['R'];C=np.rint(R*3704).astype(int);Q=np.tile(np.arange(.15,1,.1),6)
def savejson(p,z): p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(z,indent=2)+'\n')
def project(t,r=R):
 rd=r[:,:-1];p=rd.T@t/(rd.sum(1)@t);return p/p.sum()
def truths():
 u=np.full(54,1/54);joint=u.copy();joint[30]*=.75;joint/=joint.sum();fixed=u.copy().reshape(6,9);removed=fixed[3,3]*.25;fixed[3,:]+=removed/8;fixed[3,3]=u[30]*.75;return u,joint,fixed.ravel()
def wilson(k,n):
 z=1.959963984540054;p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return [float(c-h),float(c+h)]
def fraction(x):
 k=int(np.sum(x));n=len(x);return dict(k=k,n=n,rate=k/n,wilson95=wilson(k,n))
def em(n,r=R,cap=200,tol=1e-10):
 # Algebraically identical to archived w=(Rd*pi)/denom update, supporting batches.
 n=np.atleast_2d(n).astype(float);b=len(n);rd=r[...,:-1];eps=rd.sum(-1);guard_eps=eps<=0;eps=np.maximum(eps,1e-12);pi=np.full((b,54),1/54);active=np.ones(b,bool);iters=np.zeros(b,int);delta=np.full(b,np.inf);guards=np.zeros(b,bool)
 for it in range(1,cap+1):
  if rd.ndim==2:
   den=pi@rd;guard_den=den<=0;ratio=n/np.maximum(den,1e-300);new=pi*(ratio@rd.T)/eps
  else:
   den=np.einsum('bi,bij->bj',pi,rd);guard_den=den<=0;ratio=n/np.maximum(den,1e-300);new=pi*np.einsum('bj,bij->bi',ratio,rd)/eps
  sums=new.sum(1);bad=(sums<=0)|~np.isfinite(sums);guards|=active&(bad|np.any(guard_den,axis=1));new/=np.maximum(sums[:,None],1e-300)
  d=np.max(abs(new-pi),1);delta[active]=d[active];iters[active]=it;pi[active&~bad]=new[active&~bad];active&=(d>=tol)&~bad
  if not active.any():break
 diag=[dict(converged=bool(delta[i]<tol and not guards[i]),n_iter=int(iters[i]),final_max_delta=float(delta[i]),hit_iteration_cap=bool(iters[i]==cap and delta[i]>=tol),numerical_guard_triggered=bool(guards[i] or np.any(guard_eps)),finite_output=bool(np.isfinite(pi[i]).all()),sum_output=float(pi[i].sum()),min_output=float(pi[i].min()),max_output=float(pi[i].max())) for i in range(b)]
 return pi,diag
