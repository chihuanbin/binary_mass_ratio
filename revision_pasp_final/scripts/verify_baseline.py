from pathlib import Path
import hashlib,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; PROJECT=ROOT.parent; OLD=PROJECT/'pasp_cmd_recoverability_revision'; D=OLD/'data/derived'; L=D/'legacy_results'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assets=sorted(p for p in OLD.rglob('*') if p.is_file() and p.name!='.DS_Store')
 assets+=sorted(p for p in (PROJECT/'arch_code_results').rglob('*') if p.is_file() and p.suffix in ('.py','.csv','.json','.npz','.npy') and '__pycache__' not in str(p))
 rows=[]
 for p in assets:
  rel=p.relative_to(PROJECT).as_posix(); purpose='historical source/result'
  if 'overleaf_upload' in rel: purpose='user-designated manuscript baseline'
  if p.name=='calibration_200k.npz': purpose='frozen response: Fig.1, projection Fig.4; detection, reconstruction, intervals'
  if 'power_surface' in rel: purpose='archived power grid: Fig.2 and focused benchmark'
  if p.name=='exp5_operator_uncertainty.json': purpose='archived reconstruction/coverage summaries: Fig.3, Table 1; lacks counts, per-state estimates and endpoints'
  if 'proposal_audit' in rel or 'loss_budget' in rel: purpose='proposal/support accounting: Fig.5'
  rows.append(dict(path=rel,purpose=purpose,size=p.stat().st_size,sha256=sha(p),historical_frozen=True,may_modify=False))
 manifest=ROOT/'provenance/original_asset_hashes.json'
 if manifest.exists():
  frozen=json.loads(manifest.read_text())
  excluded_file=ROOT/'provenance/portable_excluded_assets.json'
  excluded=json.loads(excluded_file.read_text()) if excluded_file.exists() else []
  for entry in frozen:
   if not (PROJECT/entry['path']).exists():
    assert entry['path'] in excluded, 'Required original asset missing: '+entry['path']
    continue
   assert sha(PROJECT/entry['path'])==entry['sha256'], 'Original asset changed: '+entry['path']
 else: manifest.write_text(json.dumps(rows,indent=2)+'\n')
 text='# PROJECT INVENTORY\n\nAll listed original assets are frozen for this revision and must not be modified. SHA256 and sizes recorded before new analysis.\n\n'
 text+='## Reproduction sufficiency\n\nFig.1: response available. Fig.2: archived power arrays available. Fig.3 and Table 1: realization-level errors and coverage summaries available, but original count vectors, fitted 54-state vectors, and interval endpoints are absent. Fig.4: truth vectors and response available for deterministic projection. Fig.5: replay event summaries available; quota linkage absent. Reconstruction/interval historical aggregates can be checked; exact historical fits cannot be authenticated from absent count vectors. Newly seeded experiments will be explicitly independent.\n\n'
 text+='|Path|Purpose|Bytes|SHA256|Frozen|May modify|\n|---|---|---:|---|---|---|\n'
 text+='\n'.join(f"|{r['path']}|{r['purpose']}|{r['size']}|{r['sha256']}|yes|no|" for r in rows)
 if not (ROOT/'reports/PROJECT_INVENTORY.md').exists(): (ROOT/'reports/PROJECT_INVENTORY.md').write_text(text+'\n')
 cal=np.load(L/'calibration_200k.npz',allow_pickle=True);R=cal['R'];n=cal['n_per_state'];assert R.shape==(54,55);assert np.all(n==3704);assert n.sum()==200016;assert np.max(abs(R.sum(1)-1))<1e-14
 C=R*n[:,None];assert np.max(abs(C-np.rint(C)))<1e-8
 pear=np.load(L/'exp34_power_surface.npz',allow_pickle=True);lr=np.load(L/'exp_lr_power_surface.npz',allow_pickle=True)
 assert pear['P_det'][5,4,0,3]==205/2000;assert lr['k_alt'][5,4,0,3]==1363
 null=json.loads((D/'lr_null_validation.json').read_text());assert null['pearson']['rejections']==482;assert null['simple_vs_simple_LR']['rejections']==563
 z=json.loads((L/'exp5_operator_uncertainty.json').read_text());truth=list(csv.DictReader((D/'reconstruction_truth_54state.csv').open()))
 audit={'calibration':{'states':54,'categories':55,'rows_per_state':3704,'Ncal':200016,'counts_integer_residual':float(np.max(abs(C-np.rint(C))))},'historical_focus':{'Pearson':205,'LR':1363,'denominator':2000},'independent_null':{'Pearson':482,'LR':563,'denominator':10000},'populations':{}}
 for pop,e,b,q,c in [('S0',.01149,0,-.001317,913),('S1',.01152,.004500,-.001610,907)]:
  a=z[pop]; ev=np.mean(a['max_bin_bias']['all']);assert round(ev,5)==e;assert round(a['mean_q_bias'],6)==q
  t=np.array([float(r['t_'+pop]) for r in truth]);assert round(max(abs(t-1/54)),6)==b
  cf=a['coverage_fixed_R'];cj=a['coverage_joint_R_n'];assert sum(np.rint(np.array(cf['per_seed'])*54))==c;assert cf['successes']==c;assert cj['successes']==927
  assert cf['simultaneous_all_bins']['successes']==(1 if pop=='S0' else 0);assert cj['simultaneous_all_bins']['successes']==(1 if pop=='S0' else 0)
  audit['populations'][pop]={'mean_E':ev,'mean_q_error_archived':a['mean_q_bias'],'uniform_E':float(max(abs(t-1/54))),'fixed_coverage':c,'perturbed_coverage':927,'denominator':1080,'all54':1 if pop=='S0' else 0,'n_realizations':20}
 prop=json.loads((D/'proposal_audit/proposal_audit.json').read_text());assert prop['N_proposal']==296000;assert prop['N_model_support_rejected']==50386;assert prop['N_in_support']==245614;assert 296000-50386==245614;assert 245614-200016==45598
 audit['proposal']={'total':296000,'support_reject':50386,'in_support':245614,'calibration':200016,'unlinked_difference':45598};audit['status']='PASS';audit['missing_historical_assets']=['outer count vectors','bootstrap count vectors','per-state fitted vectors','interval endpoints','per-fit convergence traces'];audit['focus_definition']={'M1':[2.5,3.0],'q':[.4,.5],'weight':.75,'normalization':'global'}
 (ROOT/'output/baseline.json').write_text(json.dumps(audit,indent=2)+'\n')
 (ROOT/'reports/BASELINE_NUMERICAL_AUDIT.md').write_text('# BASELINE NUMERICAL AUDIT\n\nPASS: all requested baseline numbers match frozen machine-readable arrays/summaries to reported precision. No numerical mismatch. Pearson count 205 is exactly inferred from stored 0.1025 and documented denominator 2000; LR integer count is archived. Mean-q signed errors are verified against archived aggregate, not independently reconstructed from absent per-fit vectors. Coverage counts verified from realization-level counts; stored intervals retained as historical.\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n\nBLOCKED (historical exact-fit convergence authentication): original count vectors and fitted state vectors were not retained. New reruns will not be represented as historical fits.\n')
 print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
