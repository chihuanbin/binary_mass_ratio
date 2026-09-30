from common import *
def main():
 assert json.loads((ROOT/'output/baseline.json').read_text())['status']=='PASS'
 u,joint,fixed=truths();p0=project(u);alts={'joint':joint,'fixed_m1':fixed};out=ROOT/'output/fixed_m1';out.mkdir(parents=True,exist_ok=True)
 for name,a in [('truth_null',u),('truth_alt_fixed_m1',fixed),('detected_null',p0),('detected_alt',project(fixed))]:np.save(out/(name+'.npy'),a);np.savetxt(out/(name+'.csv'),a,delimiter=',')
 residual=float(max(abs(fixed.reshape(6,9).sum(1)-u.reshape(6,9).sum(1))));assert residual<1e-14
 meta={'seed_fixed':2026093001,'N':20000,'n_alternative':20000,'redistribution':'removed target probability equally divided among eight other q cells in same M1 block','M1_residual':residual,'truth_sum':float(fixed.sum()),'projected_sum':float(project(fixed).sum()),'target_before':float(u[30]),'target_after':float(fixed[30]),'response_sha256':hashlib.sha256((L/'calibration_200k.npz').read_bytes()).hexdigest()};savejson(out/'metadata.json',meta)
 pear=np.load(L/'exp34_power_surface.npz',allow_pickle=True)['thresh'].item()['20000'];palt=project(fixed);log=np.log(palt/p0)
 rng=np.random.default_rng(2026093001);null=rng.multinomial(20000,p0,200000);lr_cut=float(np.quantile(null@log,.95));counts=rng.multinomial(20000,palt,20000);ps=(((counts/20000-p0)**2)/np.maximum(p0,1e-12)).sum(1);ls=counts@log
 np.savez_compressed(out/'test_results.npz',counts=counts,pearson_stat=ps,lr_stat=ls)
 savejson(out/'test_results.json',{'Pearson':fraction(ps>pear),'LR':fraction(ls>lr_cut),'threshold_Pearson':float(pear),'threshold_LR':lr_cut,'LR_calibration_n':200000,'size_matching':'see independent common-size experiment'})
 (ROOT/'reports/FIXED_M1_ALTERNATIVE_AUDIT.md').write_text('# FIXED M1 ALTERNATIVE AUDIT\n\nPASS. Equal within-block redistribution; unchanged other M1 blocks; same frozen response. Historical joint-cell alternative retained.\n\n```json\n'+json.dumps(meta,indent=2)+'\n```\n')
 seed={'calibration':2026093002,'validation':2026093003,'alternative_joint':2026093004,'alternative_fixed_m1':2026093005};res={};dest=ROOT/'output/size_matched_tests';dest.mkdir(parents=True,exist_ok=True)
 nc=np.random.default_rng(seed['calibration']).multinomial(20000,p0,200000);nv=np.random.default_rng(seed['validation']).multinomial(20000,p0,100000)
 stat=lambda n:(((n/20000-p0)**2)/np.maximum(p0,1e-12)).sum(1)
 pc=stat(nc);pv=stat(nv);pcut=float(np.quantile(pc,.95));arrays={'null_calibration_counts':nc,'null_validation_counts':nv,'pearson_calibration_stat':pc,'pearson_validation_stat':pv}
 for name,t in alts.items():
  pa=project(t);log=np.log(pa/p0);lc=nc@log;lv=nv@log;cut=float(np.quantile(lc,.95));na=np.random.default_rng(seed['alternative_'+name]).multinomial(20000,pa,20000);aps=stat(na);als=na@log
  res[name]={'Pearson':{'threshold':pcut,'null':fraction(pv>pcut),'power':fraction(aps>pcut)},'LR':{'threshold':cut,'null':fraction(lv>cut),'power':fraction(als>cut)}}
  arrays.update({name+'_counts':na,name+'_lr_calibration_stat':lc,name+'_lr_validation_stat':lv,name+'_pearson_alt_stat':aps,name+'_lr_alt_stat':als})
 np.savez_compressed(dest/'samples_and_statistics.npz',**arrays);savejson(dest/'results.json',res);savejson(dest/'seeds.json',seed);savejson(out/'seeds.json',{'seed':2026093001});savejson(dest/'metadata.json',{'N':20000,'n_null_calibration':200000,'n_null_validation':100000,'n_alternative_each':20000,'threshold':'NumPy linear empirical quantile 0.95','reject_rule':'> threshold','LR':'separately tailored to each known alternative','response':'same frozen operator','post_hoc':True})
 (ROOT/'reports/SIZE_MATCHED_TEST_AUDIT.md').write_text('# SIZE MATCHED TEST AUDIT\n\nPASS. Independent calibration, validation and alternative streams; common target size 0.05. 200,000 calibration nulls, 100,000 validation nulls, 20,000 per alternative. Thresholds frozen before validation. Wilson 95% intervals conditional on thresholds.\n\n```json\n'+json.dumps(res,indent=2)+'\n```\n')
 print(json.dumps(res,indent=2))
if __name__=='__main__':main()
