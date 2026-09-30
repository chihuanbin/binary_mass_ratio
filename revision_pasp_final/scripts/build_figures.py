from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['DejaVu Sans','Arial','Helvetica'],'font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
def save(fig,name,arrays):
 fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
 fig.savefig(OUT/f'{name}.svg',bbox_inches='tight')
 fig.savefig(OUT/f'{name}.png',dpi=600,bbox_inches='tight')
 np.savez_compressed(OUT/f'{name}_source.npz',**arrays);plt.close(fig)
def main():
 mass=['1–1.5','1.5–2','2–2.5','2.5–3','3–4','4–5'];fig,ax=plt.subplots(figsize=(7,5),layout='constrained');im=ax.imshow(R,origin='lower',aspect='auto',interpolation='nearest',cmap='viridis',vmin=0,vmax=R.max())
 for k in range(1,6):ax.axhline(k*9-.5,color='white',lw=.5);ax.axvline(k*9-.5,color='white',lw=.5)
 ax.axvline(53.5,color='#d0d0d0',lw=1.3);ax.set_xticks([4,13,22,31,40,49,54],mass+['NOT\nOBS']);ax.set_yticks([4,13,22,31,40,49],mass);ax.tick_params(labelsize=8);ax.set_xlabel('Recovered primary-mass block (solar masses)');ax.set_ylabel('True primary-mass block (solar masses)');ax.set_title('Conditional response: 54 true × 55 recovered categories',fontsize=10);fig.colorbar(im,ax=ax,pad=.02,label='Row probability');save(fig,'response',{'R':R})
 z=np.load(L/'exp34_power_surface.npz',allow_pickle=True);lr=np.load(L/'exp_lr_power_surface.npz',allow_pickle=True);P=z['P_det'][5,:,0];LP=lr['power'][5,:,0];fig,axs=plt.subplots(1,2,figsize=(7,3.6),layout='constrained',sharey=True)
 for ax,a,label in zip(axs,[P,LP],['Pearson distance','Tailored simple LR']):
  im=ax.imshow(a,origin='lower',aspect='auto',cmap='magma',vmin=0,vmax=1,interpolation='nearest');ax.set_xticks(range(8),mass+['2–3','All'],rotation=40,ha='right',rotation_mode='anchor',fontsize=8);ax.set_yticks(range(7),[f'{v:.2f}' for v in z['amps']],fontsize=8);ax.set_title(label,fontsize=10);ax.set_xlabel('Primary-mass support (solar masses)');ax.text(3,4,f'{a[4,3]:.4f}',color='white',ha='center',va='center',fontsize=8)
 axs[0].set_ylabel('Relative weight (sampled levels)');fig.colorbar(im,ax=axs,label='Historical rejection fraction',shrink=.85);save(fig,'historical_power',{'Pearson':P,'LR':LP,'amps':z['amps']})
 a=json.loads((L/'exp5_operator_uncertainty.json').read_text());fig,axs=plt.subplots(1,2,figsize=(7,3),layout='constrained',sharey=True)
 for ax,pop in zip(axs,['S0','S1']):
  e=np.array(a[pop]['max_bin_bias']['all']);base=0 if pop=='S0' else .0044998269297334675;ax.scatter(np.arange(1,21),e,s=16,color='#376b8c');ax.axhline(base,color='.25',ls='--',lw=1,label=f'Uniform baseline: {base:.5f}');ax.set_title(f'{pop}: mean error {e.mean():.5f}',fontsize=10);ax.set_xlabel('Independent realization');ax.set_xticks([1,5,10,15,20]);ax.set_ylim(-.001,.021);ax.legend(fontsize=8,loc='upper right')
 axs[0].set_ylabel('Maximum absolute latent-fraction error');save(fig,'reconstruction',{'S0':np.array(a['S0']['max_bin_bias']['all']),'S1':np.array(a['S1']['max_bin_bias']['all'])})
 u,joint,fixed=truths();s1=u.copy();s1[30:32]*=.75;s1/=s1.sum();dt=np.array([joint-u,s1-u]);dp=np.array([project(joint)-project(u),project(s1)-project(u)]);bound=max(abs(dt).max(),abs(dp).max());fig,axs=plt.subplots(2,2,figsize=(7,5.4),layout='constrained')
 for row,title in enumerate(['One-bin statistical-detection contrast','Two-bin reconstruction contrast']):
  for col,(arr,metric) in enumerate([(dt,'Latent change'),(dp,'Projected change')]):
   ax=axs[row,col];im=ax.imshow(arr[row].reshape(6,9),origin='lower',aspect='auto',cmap='RdBu_r',vmin=-bound,vmax=bound,interpolation='nearest');ax.set_xticks([0,3,6,8],['0.15','0.45','0.75','0.95']);ax.set_yticks(range(6),mass,fontsize=8);ax.set_xlabel('q bin center');ax.set_title(metric,fontsize=10)
  axs[row,0].set_ylabel('Primary mass (solar masses)');axs[row,0].text(0,1.16,title,transform=axs[row,0].transAxes,fontsize=10)
 fig.colorbar(im,ax=axs,label='Signed probability change',shrink=.85);save(fig,'projection',{'latent':dt,'projected':dp})
 prop=json.loads((OLD/'data/derived/proposal_audit/proposal_audit.json').read_text());fig,axs=plt.subplots(1,2,figsize=(7,3),layout='constrained')
 axs[0].barh(['Support rejected','In support'],[50386,245614],color=['#777777','#376b8c']);axs[0].set_xlim(0,315000);axs[0].set_title('All 296,000 proposals',fontsize=10)
 for y,v in enumerate([50386,245614]):axs[0].text(v+3000,y,f'{v:,}',va='center',fontsize=8)
 labels=['Detected','Domain','Inverse','G cut','Uncertainty','Nonfinite'];v=[245614,0,0,0,0,0];axs[1].barh(labels,v,color='#376b8c');axs[1].set_xlim(0,325000);axs[1].set_title('245,614 in-support proposals',fontsize=10)
 for y,x in enumerate(v):axs[1].text(x+3000,y,f'{x:,}',va='center',fontsize=8)
 for ax in axs:ax.set_xlabel('Count');ax.tick_params(labelsize=8);ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
 save(fig,'accounting',{'total':np.array([296000]),'support':np.array([50386,245614]),'in_support_outcomes':np.array(v)})
 fig,ax=plt.subplots(figsize=(6.4,3.2),layout='constrained');arrays={}
 for i,name in enumerate(['R1','R2','R3']):
  z=np.load(ROOT/'output/reconstruction_extended'/f'{name}.npz');t=z['truth'];h=z['estimates'];arrays[name+'_truth']=t;arrays[name+'_estimates']=h
  for j,(c,label,color) in enumerate([(0,'200 iterations','#376b8c'),(3,'5000 iterations','#b56d37')]):
   e=np.max(abs(h[c]-t),1);mean=e.mean();low,high=np.quantile(e,[.16,.84]);ax.errorbar(i+(j-.5)*.17,mean,yerr=[[mean-low],[high-mean]],fmt='o',capsize=3,color=color,label=label if i==0 else None)
  base=np.max(abs(t-1/54));ax.plot([i-.3,i+.3],[base,base],ls='--',color='.3',label='Fixed uniform' if i==0 else None)
 ax.set_xticks(range(3),['R1: rising q','R2: high-q excess','R3: mass-dependent']);ax.set_ylim(0,.08);ax.set_ylabel('Maximum local latent-fraction error');ax.legend(fontsize=9);save(fig,'extended',arrays)
 savejson(OUT/'manifest.json',{'source_response_sha256':hashlib.sha256((L/'calibration_200k.npz').read_bytes()).hexdigest(),'sources':[str(p.relative_to(ROOT.parent)) for p in [L/'calibration_200k.npz',L/'exp34_power_surface.npz',L/'exp_lr_power_surface.npz',L/'exp5_operator_uncertainty.json',OLD/'data/derived/proposal_audit/proposal_audit.json']],'plotting':'frozen numerical values; no interpolation; all realizations','outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.suffix in ['.pdf','.npz']}})
if __name__=='__main__':main()
