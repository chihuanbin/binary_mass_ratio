from common import *
import csv

def table(name,caption,label,header,rows,notes=''):
 tex='\\begin{table}[ht]\n\\centering\n\\caption{'+caption+'}\n\\label{'+label+'}\n\\begin{tabular}{'+('l'+'c'*(len(header)-1))+'}\n\\toprule\n'+' & '.join(header)+r' \\'+'\n\\midrule\n'
 tex+='\n'.join(' & '.join(map(str,row))+r' \\' for row in rows)+'\n\\bottomrule\n\\end{tabular}\n'
 if notes:tex+='\\par\\smallskip{\\small '+notes+'}\n'
 tex+='\\end{table}\n';(ROOT/'tables'/f'{name}.tex').write_text(tex);savejson(ROOT/'tables'/f'{name}.json',{'caption':caption,'header':header,'rows':rows,'notes':notes})
def main():
 b=json.loads((ROOT/'output/baseline.json').read_text());a=json.loads((L/'exp5_operator_uncertainty.json').read_text());rows=[]
 rows += [['Mean signed $\\Delta_q$',*[f"${b['populations'][p]['mean_q_error_archived']:.6f}$" for p in ['S0','S1']]],['Mean maximum error $E_t$',*[f"{b['populations'][p]['mean_E']:.5f}" for p in ['S0','S1']]],['Fixed-uniform $E_t$', '0','0.004500'],['Fixed-uniform $\\Delta_q$','0','$-0.0004673$']]
 for key,title in [('coverage_fixed_R','Fixed response'),('coverage_joint_R_n','All-cell response perturbation')]:
  rows.append([title+' coverage',*[f"{a[p][key]['successes']}/1080" for p in ['S0','S1']]])
  rows.append([title+' cluster 95\\% interval',*[f"{a[p][key]['cluster_seed_bootstrap_95'][0]:.3f}--{a[p][key]['cluster_seed_bootstrap_95'][1]:.3f}" for p in ['S0','S1']]])
  rows.append([title+' all-54 containment',*[f"{a[p][key]['simultaneous_all_bins']['successes']}/20" for p in ['S0','S1']]])
 table('reconstruction','Historical reconstruction and interval behavior.','tab:reconstruction',['Quantity','S0','S1'],rows,'20 independent realizations per truth, 20,000 systems each; nominal marginal intervals use the 16th/84th percentiles. Historical fits are finite-iteration outputs; convergence was not recorded.')
 res=json.loads((ROOT/'output/size_matched_tests/results.json').read_text());rows=[]
 for name,label in [('joint','Joint-cell'),('fixed_m1','Fixed-$M_1$ $q$ redistribution')]:
  rows.append([label,*[f"{res[name][test][q]['rate']:.5f}" for test in ['Pearson','LR'] for q in ['null','power']]])
 table('detection','Focused statistical-detection benchmark with a common 5\\% target size.','tab:detection',['Alternative','Pearson size','Pearson power','LR size','LR power'],rows,'Thresholds: 200,000 calibration nulls; size: 100,000 independent validation nulls; power: 20,000 alternatives per row. Each LR is tailored to its stated alternative. Wilson 95\\% power intervals: joint-cell Pearson 0.1045--0.1131, LR 0.6588--0.6719; fixed-$M_1$ Pearson 0.0819--0.0897, LR 0.4961--0.5099. Size intervals and counts are supplied in the machine-readable table source.')
 z=json.loads((ROOT/'output/response_prior/results.json').read_text());rows=[]
 for method,title in [('coverage_fixed_R','Fixed (historical)'),('coverage_joint_R_n','All-cell (historical)')]:
  rows.append([title,*[f"{a[p][method]['mean']:.3f}" for p in ['S0','S1']],*[f"{a[p][method]['mean_interval_width']:.5f}" for p in ['S0','S1']]])
 rows.append(['Zero-preserving (new)',*[f"{z['structural_zero_200'][p]['coverage']:.3f}" for p in ['S0','S1']],*[f"{z['structural_zero_200'][p]['mean_width']:.5f}" for p in ['S0','S1']]])
 table('response_prior','Nominal 68\\% marginal interval procedures at the 200-iteration cap.','tab:prior',['Procedure','S0 coverage','S1 coverage','S0 width','S1 width'],rows,'Width is the mean upper-minus-lower latent fraction. The four historical widths are archived aggregate means in the retained interval-summary JSON; historical endpoint arrays were not retained. The new procedure uses independent realizations, so differences from historical procedures are unpaired. Its realization-cluster 95\\% intervals are 0.873--0.925 (S0) and 0.852--0.904 (S1); all-54 containment is 0/20 for both.')
 cat=list(csv.DictReader((ROOT.parent/'arch_code_results/data/hunt24/params.csv').open()));support={'Pleiades':'1--5','Hyades':'1--2.5','Praesepe':'1--3','M67':'1--1.5','NGC188':'None','NGC6791':'None'}
 rows=[[r['cluster'],f"{float(r['log_age']):.4f}",f"{float(r['feh']):.2f}",f"{float(r['dist_pc']):.2f}",f"{float(r['AV']):.4f}",support[r['cluster']]] for r in cat]
 table('catalogue','Centers in the six-cluster master catalogue and nominal support subsets.','tab:catalogue',['Cluster','$\\log_{10}(\\mathrm{age/yr})$','[Fe/H]','$d$ (pc)','$A_V$ (mag)','$M_1/\\msun$'],rows,'Support denotes the mass intervals for which the calibration generator selects that cluster; NGC 188 and NGC 6791 are in the master catalogue but not in any selected subset.')
 savejson(ROOT/'tables/detection_full_source.json',res)
if __name__=='__main__':main()
