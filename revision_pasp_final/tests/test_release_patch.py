import csv,json,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def test_historical_width_aggregate_provenance():
 e=json.loads((ROOT/'provenance/historical_width_evidence.json').read_text());h=json.loads((ROOT.parent/e['sources'][0]['path']).read_text())
 for r in e['values']:
  assert r['classification'].startswith('A:')
  assert h[r['population']][r['procedure']]['mean_interval_width']==r['value']
  assert f"{r['value']:.5f}"==r['table_rounding']
 assert 'archived aggregate means' in (ROOT.parent/'manuscript/pasp_submission_final/tables/response_prior.tex').read_text()
def test_em_objective_replay_monotonicity():
 s=json.loads((ROOT/'output/em_objective/summary.json').read_text())
 assert s['status']=='PASS' and s['n_fits']==41540 and s['n_decreases_beyond_tolerance']==0
 assert s['max_replay_estimate_residual']<1e-10 and s['rng_draws']==0
 rows=list(csv.DictReader((ROOT/'output/em_objective/per_fit_objective.csv').open()));assert len(rows)==41540
 for r in rows:
  assert int(r['decreases_beyond_tolerance'])==0
  assert np.isfinite([float(r[f'final_max_delta_{cap}']) for cap in [200,500,1000,5000]]).all()
 assert abs(np.median([float(r['final_max_delta_5000']) for r in rows])-s['final_max_delta_5000_quantiles']['0.5'])<1e-15
 for r in s['inputs']:assert hashlib.sha256((ROOT.parent/r['path']).read_bytes()).hexdigest()==r['sha256']
def test_target_size_wording():
 s=(ROOT.parent/'manuscript/pasp_submission_final/main.tex').read_text()
 assert r'common 5\% target size' in s and 'discrete multinomial count and test-statistic distributions' in s
 assert 'exactly size matched' not in s
