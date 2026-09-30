import sys,json,csv,hashlib,re
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));from common import R,C,Q,L,OLD,truths,project
FINAL=ROOT.parent/'manuscript/pasp_submission_final';B=json.loads((ROOT/'output/baseline.json').read_text());SM=json.loads((ROOT/'output/size_matched_tests/results.json').read_text());PR=json.loads((ROOT/'output/response_prior/results.json').read_text());H=json.loads((L/'exp5_operator_uncertainty.json').read_text())
def test_response_shape():assert R.shape==(54,55)
def test_response_row_normalization():assert np.max(abs(R.sum(1)-1))<1e-14
def test_state_count():assert 6*9==54==len(Q)
def test_calibration_size():assert C.sum()==200016;assert np.all(C.sum(1)==3704)
def test_fixed_m1_marginal_preserved():u,j,f=truths();assert np.max(abs(u.reshape(6,9).sum(1)-f.reshape(6,9).sum(1)))<1e-14;assert f[30]==u[30]*.75
@pytest.mark.parametrize('t',truths())
def test_truth_vectors_normalized(t):assert np.isfinite(t).all();assert (t>=0).all();assert abs(t.sum()-1)<1e-14
@pytest.mark.parametrize('t',truths())
def test_detected_vectors_normalized(t):p=project(t);assert np.isfinite(p).all();assert (p>0).all();assert abs(p.sum()-1)<1e-14
def test_baseline_numbers_match():assert B['status']=='PASS';assert B['historical_focus']=={'Pearson':205,'LR':1363,'denominator':2000}
def test_null_rates_match_archive():z=json.loads((OLD/'data/derived/lr_null_validation.json').read_text());assert z['pearson']['rejections']==482;assert z['simple_vs_simple_LR']['rejections']==563
@pytest.mark.parametrize('p,e,q',[('S0',.01149,-.001317),('S1',.01152,-.001610)])
def test_reconstruction_numbers_match_archive(p,e,q):assert round(np.mean(H[p]['max_bin_bias']['all']),5)==e;assert round(H[p]['mean_q_bias'],6)==q
@pytest.mark.parametrize('p,c,s',[('S0',913,1),('S1',907,0)])
def test_coverage_counts_match_archive(p,c,s):assert H[p]['coverage_fixed_R']['successes']==c;assert H[p]['coverage_joint_R_n']['successes']==927;assert H[p]['coverage_fixed_R']['simultaneous_all_bins']['successes']==s

def test_proposal_accounting():p=B['proposal'];assert p['total']-p['support_reject']==p['in_support'];assert p['in_support']-p['calibration']==45598

def test_figure_source_values():
 assert np.array_equal(np.load(ROOT/'figures/response_source.npz')['R'],R)
 z=np.load(ROOT/'figures/historical_power_source.npz');assert np.array_equal(z['Pearson'],np.load(L/'exp34_power_surface.npz',allow_pickle=True)['P_det'][5,:,0]);assert np.array_equal(z['LR'],np.load(L/'exp_lr_power_surface.npz',allow_pickle=True)['power'][5,:,0]);assert z['Pearson'][4,3]==.1025;assert z['LR'][4,3]==.6815
 e=np.load(ROOT/'figures/reconstruction_source.npz');assert np.array_equal(e['S0'],H['S0']['max_bin_bias']['all']);assert np.array_equal(e['S1'],H['S1']['max_bin_bias']['all'])
 d=np.load(ROOT/'figures/projection_source.npz');assert np.max(abs(d['latent'].sum(1)))<1e-14;assert np.max(abs(d['projected'].sum(1)))<1e-14
 a=np.load(ROOT/'figures/accounting_source.npz');assert a['support'].sum()==296000;assert a['in_support_outcomes'].sum()==245614

def test_table_source_values():
 z=json.loads((ROOT/'tables/detection.json').read_text())
 for row,name in zip(z['rows'],['joint','fixed_m1']):assert row[1:]==[f"{SM[name][t][q]['rate']:.5f}" for t in ['Pearson','LR'] for q in ['null','power']]
 z=json.loads((ROOT/'tables/response_prior.json').read_text());assert z['rows'][2][1:3]==[f"{PR['structural_zero_200'][p]['coverage']:.3f}" for p in ['S0','S1']]
 z=json.loads((ROOT/'tables/reconstruction.json').read_text());assert z['rows'][1][1:]==[f"{B['populations'][p]['mean_E']:.5f}" for p in ['S0','S1']]

def test_no_nan_in_final_outputs():
 for p in (ROOT/'output').rglob('*.npz'):
  with np.load(p) as z:
   for k in z.files:assert np.isfinite(z[k]).all(),f'{p}:{k}'

def test_no_negative_probabilities():
 for p in (ROOT/'output').rglob('*.npy'):assert np.min(np.load(p))>=0
 for p in (ROOT/'output/response_prior').glob('*.npz'):
  with np.load(p) as z:
   for k in ['truth','outer_estimates','response_positive_values','bootstrap_estimates','lo','hi']:assert np.min(z[k])>=0
   assert np.max(abs(z['bootstrap_estimates'].sum(-1)-1))<1e-14

def test_manuscript_key_numbers_match_generated_tables():
 text=(FINAL/'main.tex').read_text()
 for n in ['0.1025','0.6815','0.0482','0.0563','0.01149','0.01152','0.004500','-0.001317','-0.001610','972/1,080','949/1,080','40,040','0.07292','0.03312']:assert n in text
 for a in SM.values():
  for t in ['Pearson','LR']:
   for stat in ['null','power']:assert f"{a[t][stat]['rate']:.5f}" in text
 assert '87.9\\%--90.0\\%' in text

def test_original_asset_hashes_unchanged():
 for e in json.loads((ROOT/'provenance/original_asset_hashes.json').read_text()):
  missing_file=ROOT/'provenance/portable_excluded_assets.json'
  if not (ROOT.parent/e['path']).exists():
   assert missing_file.exists() and e['path'] in json.loads(missing_file.read_text()), e['path']
   continue
  h=hashlib.sha256()
  with (ROOT.parent/e['path']).open('rb') as f:
   for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
  assert h.hexdigest()==e['sha256'],e['path']

def test_size_matched_samples_and_statistics():
 with np.load(ROOT/'output/size_matched_tests/samples_and_statistics.npz') as z:
  assert len(z['null_calibration_counts'])==200000;assert len(z['null_validation_counts'])==100000
  for a in ['joint','fixed_m1']:
   assert len(z[a+'_counts'])==20000;assert np.all(z[a+'_counts'].sum(1)==20000)
   for t,key in [('Pearson','pearson'),('LR','lr')]:
    alt=z[a+'_'+key+'_alt_stat'];nv=z['pearson_validation_stat'] if key=='pearson' else z[a+'_lr_validation_stat'];cut=SM[a][t]['threshold'];assert (alt>cut).sum()==SM[a][t]['power']['k'];assert (nv>cut).sum()==SM[a][t]['null']['k']

def test_interval_endpoints_reproduce_coverage():
 for p in ['S0','S1']:
  successes=0;width=[];alls=0
  for f in (ROOT/'output/response_prior').glob(p+'_*.npz'):
   with np.load(f) as z:
    lo,hi=np.quantile(z['bootstrap_estimates'][0],[.16,.84],axis=0);assert np.max(abs(lo-z['lo']))<1e-14;assert np.max(abs(hi-z['hi']))<1e-14;ind=(lo<=z['truth'])&(z['truth']<=hi);successes+=ind.sum();width.extend(hi-lo);alls+=ind.all()
  assert successes==PR['structural_zero_200'][p]['k'];assert alls==PR['structural_zero_200'][p]['all54'];assert abs(np.mean(width)-PR['structural_zero_200'][p]['mean_width'])<1e-14

def test_structural_zeros_preserved():
 for f in (ROOT/'output/response_prior').glob('*.npz'):
  with np.load(f) as z:
   rr=np.zeros((1000,54,55));rr[:,z['response_i'],z['response_j']]=z['response_positive_values'];assert np.all(rr[:,C==0]==0);assert np.max(abs(rr.sum(-1)-1))<1e-14

def test_convergence_is_not_falsely_claimed():
 rows=list(csv.DictReader((ROOT/'output/convergence/per_fit_diagnostics.csv').open()));assert len(rows)==40040*4
 for r in rows:
  assert (r['converged']=='True')==(float(r['final_max_delta'])<1e-10 and r['numerical_guard_triggered']=='False')
  assert r['finite_output']=='True';assert abs(float(r['sum_output'])-1)<1e-14

def test_optimized_em_matches_archived():
 sys.path.insert(0,str(OLD/'src/legacy_scripts'));from unfolding import dagostini_observable_only
 from run_convergence_response import sparse_fits
 ii,jj=np.nonzero(C[:,:-1]);rng=np.random.default_rng(73294);n=rng.multinomial(20000,project(truths()[0]),size=4);v,_,_,_=sparse_fits(n,np.tile(R[ii,jj],(4,1)),ii,jj,np.array([200]))
 for i in range(4):assert np.max(abs(v[0,i]-dagostini_observable_only(R,n[i])))<1e-10

def test_reader_facing_internal_ids_absent():
 text=(FINAL/'main.tex').read_text();assert not any(s in text for s in ['S3-LR-01','Experiment 5','closure','confirmatory gate','internal ID'])

def test_reference_keys_resolve():
 text=(FINAL/'main.tex').read_text();keys=set(re.findall(r'@\w+\{([^,]+)',(FINAL/'references.bib').read_text()));cites=set(k for c in re.findall(r'\\cite\w*\{([^}]+)\}',text) for k in c.split(','));assert cites<=keys;assert len(keys)==(26 if 'BenchmarkRelease' in cites else 25)

def test_extended_metrics():
 summary=json.loads((ROOT/'output/reconstruction_extended/summary.json').read_text())
 for name in ['R1','R2','R3']:
  with np.load(ROOT/'output/reconstruction_extended'/f'{name}.npz') as z:
   t=z['truth'];assert abs(t.sum()-1)<1e-14;assert np.max(abs(t.reshape(6,9).sum(1)-1/6))<1e-14;assert len(z['counts'])==500;assert abs(np.max(abs(z['estimates'][0]-t),1).mean()-summary[name]['EM200']['E'])<1e-14

def test_latex_compile_consistency():
 log_path=FINAL/'main.log';log=(log_path if log_path.exists() else ROOT/'reports/latex_compile.log').read_text();assert not re.search(r'undefined|Overfull|Emergency stop|Fatal error',log,re.I);assert (FINAL/'main.pdf').stat().st_size>100000
