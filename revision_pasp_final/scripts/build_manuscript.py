from common import *
import shutil,re,difflib
FINAL=ROOT.parent/'manuscript/pasp_submission_final'
def main():
 source=OLD/'overleaf_upload/manuscript/main.tex';old=source.read_text();a=json.loads((ROOT/'output/response_prior/results.json').read_text());sm=json.loads((ROOT/'output/size_matched_tests/results.json').read_text());ext=json.loads((ROOT/'output/reconstruction_extended/summary.json').read_text())
 head=old[:old.index('\\begin{abstract}')];head=re.sub(r'% Author list.*?\\author',r'\\author',head,flags=re.S);head=head.replace('\\shortauthors{Chi et al.}','\\shortauthors{Chi \\& Wang}').replace('\\submitjournal{PASP}', '% Target journal: PASP')
 intro=r'''\section{Introduction}
Unresolved binaries shift a system's position in a color--magnitude diagram (CMD) through the combined fluxes of their components \citep{Li2020,Li2022,Albrow2022}. Binary fractions and mass-ratio distributions constrain stellar formation and cluster evolution \citep{Duchene2013,Raghavan2010,Moe2017}. Their photometric inference depends on the stellar model, population parameterization, measurement errors, and observational selection. Interpreting inferred structure therefore requires validation matched to the quantity being claimed.

Existing CMD studies use mixture and generative approaches and have performed substantial mock-data validation. Li et al. modeled unresolved binaries in NGC 3532 using Gaia DR2 \citep{Li2020}; Li \& Shao introduced MiMO, validated it on 1,000 mock clusters, and analyzed ten Gaia EDR3 clusters \citep{Li2022}. Generative models have inferred binary populations in M67, Hyades, Praesepe, and six young open clusters \citep{Albrow2022,Albrow2024,Alexander2025}. Broader cluster surveys have measured photometric binary fractions and their variation with cluster properties \citep{Milone2012,Cordoni2023,Donada2023,Jadhav2021}. Gaia DR3 cluster analyses also use joint stellar and cluster-parameter models such as BASE-9 with PARSEC isochrones \citep{Chi2025a}. These approaches have different estimands and validation designs.

Mock recovery comprises at least three distinct operating characteristics. Statistical detection measures rejection probability for a specified population departure at a stated false-positive rate. Reconstruction measures the discrepancy between an estimated latent distribution and the known input distribution. Interval calibration measures repeated-sampling containment of the truth. A sensitive test need not reconstruct the departure accurately, and agreement in a global summary need not establish local accuracy or calibrated intervals.

We construct a controlled finite-state CMD benchmark in which these properties can be measured separately using the same response operator. We retain historical power, reconstruction, and coverage experiments, then add a mass-ratio redistribution that holds the primary-mass marginal fixed, threshold calibration at a common target size, iteration-cap diagnostics, and a response-perturbation comparison preserving empirical zeros. Three analytically specified nonuniform truths extend the reconstruction test beyond the original uniform and near-uniform cases.

The contribution is a reusable validation framework for a discrete population-inference workflow. Its numerical results are conditional on the stated grid, nuisance sampling, support, selection model, and inverse fits supplied the true nuisance parameters. They support methodological conclusions within this benchmark rather than an astrophysical population measurement or a general ranking of CMD algorithms. Throughout, ``statistical detection'' means detecting a distributional departure; photometric detection and selection refer to observational retention.

'''
 benchmark=old[old.index('\\section{Benchmark and forward model}'):old.index('\\section{Response operator and population inference}')]
 benchmark=benchmark.replace('\\section{Benchmark and forward model}','\\section{Benchmark construction}').replace('\\subsection{Latent grid, stellar models, and nuisance sampling}','\\subsection{Latent state space}')
 benchmark=benchmark.replace('The forward code reads', '\\subsection{Stellar and CMD forward model}\nThe forward code reads')
 benchmark=benchmark.replace('Cluster nuisance parameters are drawn','\\subsection{Nuisance sampling}\nCluster nuisance parameters are drawn')
 benchmark=benchmark.replace('For each component,','\\subsection{Flux addition, extinction, and photometric noise}\nFor each component,')
 benchmark=benchmark.replace('For band $b$ and nuisance vector', 'These coefficients are fixed implementation choices in this benchmark; the retained source gives no extinction-law derivation or traceable $R_V$ assumption. They are not fitted or claimed to apply to arbitrary stellar spectra. For band $b$ and nuisance vector')
 benchmark=benchmark.replace('The inverse fit uses a grid','\\subsection{Inverse grid}\nThe inverse fit uses a grid')
 benchmark=benchmark.replace('The master catalogue supplies nuisance-parameter centers;', 'The six entries are Pleiades, Hyades, NGC 188, Praesepe, M67, and NGC 6791; their centers and selected mass intervals are listed in Appendix~\\ref{app:catalogue}. The master catalogue supplies nuisance-parameter centers;')
 response=old[old.index('Let $R_{ij}'):old.index('Figure~\\ref{fig:response} displays')]
 start=response.index('The implementation stops');end=response.index('The zero-efficiency fallback')
 response=response[:start]+r'''The original implementation stops at a maximum component change below $10^{-10}$ or a cap of 200 iterations. Nonpositive detected denominators and efficiencies are replaced by $10^{-300}$ and $10^{-12}$, respectively. A nonpositive or nonfinite updated total returns the last iterate. We treat this as a finite-iteration estimator, since the historical outputs did not store convergence diagnostics. New independent fits record iteration count, final component change, cap hitting, numerical guards, and simplex validity at caps of 200, 500, 1,000, and 5,000. '''+response[end:]
 methods=r'''\section{Response and validation experiments}
\subsection{Response operator and finite-iteration inversion}
'''+response+r'''
The frozen response has unit detected efficiency in every row. It summarizes the within-bin draws, support conditioning, nuisance sampling, photometric-error model, sensitivity cuts, and inverse grid. Figure~\ref{fig:response} displays its probabilities directly.
\begin{figure}[ht]
\centering\includegraphics[width=0.88\columnwidth]{figures/response.pdf}
\caption{Frozen $54\times55$ response. Both axes are flattened in primary-mass-block then $q$-bin order. Within each mass block the nine $q$ bins span 0.1--1.0. The neutral separator isolates the zero-probability NOT\_OBSERVED column. This is a conditional, unit-efficiency calibration response; it is not a Gaia completeness function.}
\label{fig:response}\end{figure}

\subsection{Statistical-detection experiments}
The null is uniform over the 54 latent states. The historical one-cell alternative multiplies the fraction in $2.5\leq M_1/\msun<3.0$, $0.4\leq q<0.5$ by $w=0.75$, then globally normalizes. Its affected primary-mass marginal becomes $8.75/53.75$, instead of $1/6$. Thus this alternative changes the joint $(M_1,q)$ distribution.

For the new fixed-$M_1$ alternative, the same target cell is reduced to $0.75/54$. The removed fraction $0.25/54$ is divided equally among the eight other $q$ cells of that mass block. All other mass blocks are unchanged. Both alternatives are projected through the same frozen response; the new design preserves every primary-mass marginal to machine precision.

The Pearson-distance statistic is $X^2=\sum_j(\hat p_j-p_{0,j})^2/\max(p_{0,j},10^{-12})$, with $\hat p_j=n_j/N$. The simple-versus-simple LR is $T=\sum_j n_j\log[p_{1,j}/p_{0,j}]$ for the stated alternative. These are upper-tail tests, rejecting strictly above an empirical 95th-percentile threshold. The projected probabilities here are positive in all 54 detected categories. Historical thresholds used 2,000 null vectors; the archived power runs used 2,000 alternatives per sampled cell. The LR threshold was specific to each alternative, while Pearson shared a threshold at each $N$. Historical power runs were independent between tests, so their rejection indicators cannot be paired retrospectively.

The new common-size comparison calibrates thresholds on 200,000 null count vectors, freezes them, validates their size on 100,000 independent null vectors, and evaluates 20,000 independent alternatives per design at $N=20,000$. Each test uses the same vectors within a design, but the calibration, validation, and alternative samples have distinct PCG64 seeds. The LR is separately tailored to each known alternative. We report realized sizes, powers, and Wilson 95\% intervals \citep{Brown2001}; intervals condition on the calibrated thresholds. These additions are post hoc benchmark extensions, not preregistered analyses.

\subsection{Reconstruction experiments}
Historical truths are S0 (uniform) and S1 (two adjacent $q$ cells, $[0.4,0.5)$ and $[0.5,0.6)$, downweighted by 0.75 within the same $2.5$--$3.0\,\msun$ block, then globally normalized). S1 differs from the one-cell statistical-detection alternative and is not a fixed-$M_1$ design. Each truth has 20 archived validation realizations at $N=20,000$ from a sequential PCG64 stream seeded with 999. We compare the returned 200-iteration estimate against a fixed-uniform estimator, $\hat t_i=1/54$.

For realization $s$, $E_{t,s}=\max_i|\hat t_{s,i}-t_i|$ measures maximum local error. The discretized mean-$q$ error is $\Delta_{q,s}=\sum_i(\hat t_{s,i}-t_i)q_i$, where $q_i$ is the bin center. These summaries measure different discrepancies. Within-bin uniform sampling makes the design-expected continuous mean equal to this bin-center mean. Historical assets retain every $E_{t,s}$ and aggregate signed $\Delta_q$, but not the original count vectors or fitted state vectors. Exact historical-fit convergence authentication is therefore unavailable. The new convergence experiment uses independently drawn S0/S1 counts and preserves these counts and all fitted vectors.

We also declare three nonuniform truths before generating any realizations: R1 has conditional $q$ weights proportional to bin center; R2 has weights $1+2\exp\{-[(q_i-0.9)/0.1]^2/2\}$, a smooth synthetic high-$q$ excess; R3 uses weights $1.1-q_i$ in the first three mass blocks and $q_i$ in the last three. Each conditional distribution is normalized and each mass block has probability $1/6$. There are 500 independent realizations per truth at $N=20,000$. We archive maximum error, $L_1$ error, RMSE, Jensen--Shannon divergence (natural logarithms), population mean-$q$ error, and conditional mean-$q$ errors for every mass block. The primary estimator remains the 200-iteration procedure; longer caps are sensitivity diagnostics, not replacements for historical results.

\subsection{Interval-calibration procedures}
For each historical realization, 1,000 parametric count-bootstrap samples were drawn around $p(\hat t,R)$ and inverted with either fixed $R$ or an independently perturbed response. The all-cell procedure draws each 55-category row as
\begin{equation}
R^*_{i,\cdot}\sim\operatorname{Dirichlet}(C_{i1}+1,\ldots,C_{i,55}+1),\qquad C_{ij}=3704R_{ij}.
\end{equation}
The full-precision response resolves the calibration counts to integers within $1.14\times10^{-13}$; every row has 3,704 counts. The all-cell procedure assigns a pseudo-count to every zero cell, including NOT\_OBSERVED. Its mean uniform-mixture weight is $55/3759=0.01463$, and its mean missing probability is $1/3759=2.66\times10^{-4}$. We call this a response-perturbation procedure, not a hierarchical posterior.

The new comparison draws Dirichlet$(C_{ij}+1)$ only over cells with $C_{ij}>0$. All empirical zeros, including NOT\_OBSERVED, remain exactly zero. Empirical zeros are an imposed support constraint for this sensitivity test, not proof that the underlying probabilities vanish. The procedure uses 20 independent S0/S1 realizations and 1,000 bootstrap inversions per realization. Nominal fitted probabilities generate bootstrap counts independently of the response used for inversion. It therefore matches the historical resampling construction but uses new realizations; comparisons with historical procedures are unpaired. All new outer and bootstrap counts, response draws, estimates, endpoints, and convergence diagnostics are retained.

The 16th and 84th bootstrap percentiles define nominal 68\% marginal intervals. Coverage pools the 54 state indicators within each realization; a 95\% realization-cluster interval resamples the 20 independent realization-level fractions 20,000 times \citep{Efron1979}. All-54 containment is the fraction of realizations containing every truth in its marginal interval, and is not a nominal simultaneous confidence region. Mean interval width is the average endpoint difference. We repeat inversions at longer caps using the same bootstrap counts and response draws; the bootstrap-generating population remains the 200-iteration fit. This isolates inversion-cap sensitivity rather than redefining a fully recalibrated long-iteration interval method.

'''
 results=r'''\section{Validation results}
\subsection{Detectability of specified population departures}
At the historical one-cell alternative and $N=20,000$, Pearson rejected 205/2,000 alternatives (0.1025), while the tailored LR rejected 1,363/2,000 (0.6815). An independent focused null check at the original thresholds rejected 482/10,000 (0.0482) and 563/10,000 (0.0563), respectively. These raw powers should therefore not be interpreted as an exactly size-matched comparison. Figure~\ref{fig:power} retains the historical sampled-grid results as a description of those frozen thresholds.
\begin{figure}[ht]
\centering\includegraphics[width=\columnwidth]{figures/historical_power.pdf}
\caption{Historical power at $N=20,000$ for one $q$ bin, $0.4\leq q<0.5$. Each sampled cell uses 2,000 alternative draws and the original 2,000-draw null thresholds. The relative weights are categorical sampled levels, followed by global normalization. Both panels share the same color scale; no unsimulated values are interpolated. Annotated values are the historical focus cell. The common-target-size comparison is reported separately in Table~\ref{tab:detection}.}
\label{fig:power}\end{figure}

The new threshold calibration targets 5\% for both tests. Independent validation sizes are 0.05148 for Pearson and 0.05199/0.05098 for the joint/fixed-$M_1$ LR rules. The residual realized-size differences reflect finite calibration and validation sampling and the discrete multinomial count and test-statistic distributions; the comparison uses a common 5\% target size rather than asserting exact equality of realized sizes. Pearson's size interval is 0.05013--0.05287; the LR intervals are 0.05063--0.05338 and 0.04963--0.05236. These exclude threshold-recalibration uncertainty.

For the joint-cell departure, the new powers are 0.10875 (Pearson) and 0.66540 (LR). Holding every $M_1$ marginal fixed gives powers 0.08570 and 0.50300. The largest numerical change in an $M_1$ marginal is below $10^{-14}$. Thus the specified $q$ redistribution remains detectable by its tailored LR, although sensitivity is lower than for the historical joint departure. The result does not identify which unknown $q$ features would be detectable or demonstrate state-by-state reconstruction.
\input{tables/detection.tex}

\subsection{Reconstruction accuracy and iteration-cap sensitivity}
Historical mean maximum errors were 0.01149 (S0) and 0.01152 (S1), compared with fixed-uniform errors of 0 and 0.004500. Mean signed bin-center $q$ errors were $-0.001317$ and $-0.001610$, while the fixed output gives 0 and $-0.0004673$. The true overall means are 0.55 and 0.550467; S1 has conditional mean 0.55294 in the affected mass block and 0.55 elsewhere. These uniform and near-uniform designs deliberately make the fixed estimator a strong baseline. The historical finite-iteration estimator does not improve these error summaries over that baseline (Table~\ref{tab:reconstruction}, Figure~\ref{fig:reconstruction}).
\input{tables/reconstruction.tex}
\begin{figure}[ht]
\centering\includegraphics[width=\columnwidth]{figures/reconstruction.pdf}
\caption{All 20 historical maximum local errors per truth. Dashed lines show the analytic fixed-uniform baseline evaluated against the same truth, including exactly zero for S0. These are finite-iteration estimates. Mean-$q$ signed errors and interval coverage are distinct quantities, reported in Table~\ref{tab:reconstruction}.}
\label{fig:reconstruction}\end{figure}

The independent convergence audit contains 40 outer fits and 40,000 structural-zero bootstrap fits. All 40,040 reach the 200-iteration cap without satisfying $10^{-10}$; all remain above tolerance at 500, 1,000, and 5,000 iterations. No numerical guard exits occur, and every returned vector is finite and normalized. The maximum $\|\hat t_{200}-\hat t_{5000}\|_\infty$ across these fits is 0.07292, and the mean is 0.03312. The optimized update matches the archived implementation within $10^{-10}$ on ten independent test vectors; this checks implementation equivalence, not historical-fit reproduction. The observed-data multinomial log-likelihood was non-decreasing in all 41,540 independent diagnostic fits (including the nonuniform extension), supporting implementation consistency without implying convergence.

For the new S0/S1 outer fits, mean $E_t$ changes from 0.01138/0.01163 at 200 iterations to 0.03936/0.03567 at 5,000. Mean signed $\Delta_q$ changes from $+0.000684/-0.001501$ to $+0.001677/-0.003271$. The longer fits are not declared converged and do not replace any historical estimate. The 200-step procedure is consequently an iteration-limited estimator whose stopping rule materially affects the result.

EXTENDED_RESULTS

\subsection{Interval calibration and empirical-zero preservation}
Historical fixed-response coverage was 913/1,080 (0.845) for S0 and 907/1,080 (0.840) for S1; all-cell response perturbation covered 927/1,080 (0.858) for both. The realization-cluster uncertainty intervals and all-54 containment counts are given in Table~\ref{tab:reconstruction}. Neither procedure attained the nominal 68\% marginal level in these experiments.

With empirical zeros preserved, the independent new coverage is 972/1,080 (0.900) for S0 and 949/1,080 (0.879) for S1 (Table~\ref{tab:prior}). The corresponding realization-cluster 95\% intervals are 0.873--0.925 and 0.852--0.904. All-54 containment is 0/20 for both. The overcoverage therefore persists when zero cells are not filled. Because the historical and new outer realizations differ, their coverage difference does not isolate a response-prior effect. None of the three procedures is validated as a calibrated posterior or confidence construction beyond these experiments.
\input{tables/response_prior.tex}

Using the same new bootstrap draws and changing only the inversion cap to 5,000 increases coverage to 1,067/1,080 (S0) and 1,066/1,080 (S1), with mean widths 0.03006/0.03024 rather than 0.00898/0.00899. All-54 containment becomes 13/20 for both. These are sensitivity diagnostics around the original fitted generating population, not a separately recalibrated interval method. They further show that interval labels and widths depend on the complete numerical procedure.

\subsection{Response projection and proposal diagnostics}
Appendix~\ref{app:projection} projects the separate one-cell statistical-detection and two-cell reconstruction contrasts through the frozen response. The projections preserve normalization and describe those declared contrasts only. Appendix~\ref{app:proposal} accounts for support rejections and observational model outcomes in the proposal replay. Its unit observed retention conditional on support is an outcome of that replay, not a measurement of Gaia completeness.

\section{Discussion}
\subsection{Detectability is alternative-specific}
The LR uses the exact supplied projected alternative. At the common 5\% target size it has higher sensitivity than Pearson for the two declared departures, with closely matched but not identical realized false-positive rates. This is a statement about specified simple alternatives, not an algorithm ranking or a search over unknown structures. The fixed-$M_1$ experiment separates sensitivity to a $q$ redistribution from sensitivity to the globally normalized joint-cell departure. Its lower power also illustrates why the full alternative definition belongs with every detection claim.

\subsection{Reconstruction requires evidence beyond detection}
A test can reject the null without estimating every latent fraction accurately. In S0/S1, small mean-$q$ errors coexist with maximum local errors larger than the fixed-uniform baseline. The nonuniform extension establishes a broader but still finite set of point-estimator comparisons. It does not validate arbitrary distributions. The convergence audit also makes the iteration cap part of the estimand-defining procedure: the published finite-step vectors are not authenticated converged maximum-likelihood solutions. Longer iteration increases local error in the new S0/S1 runs, so satisfying a stricter numerical criterion and improving statistical accuracy are separate objectives.

\subsection{Interval labels require repeated-sampling calibration}
Nominal 68\% marginal intervals overcover on average under all tested 200-step procedures. Preservation of empirical zeros leaves substantial overcoverage, so filling zero cells is not its sole source in this benchmark. The unpaired historical comparison cannot quantify a causal prior effect. The cap sensitivity further shows that resampling, response perturbation, and inversion settings jointly determine interval behavior. Marginal coverage and all-state containment answer different questions, and neither justifies a global calibration claim.

\subsection{Scope and transfer to cluster studies}
The response is conditional on the stellar model, finite bins, uniform within-bin draws, cluster-dependent nuisance sampling, support conditioning, and observational cuts. Nuisance distributions vary with primary mass, so mass-block differences are not isolated physical mass trends. Inverse fits receive true age, distance, extinction, and metallicity; the benchmark does not test their simultaneous inference. It includes no spatial membership-completeness model, crowding, triples, external response mismatch, or Gaia completeness calibration. Tie-breaking rates and forward-interpolated versus inverse-snapped grid effects were not isolated by new ablations. The validation framework can be applied to another cluster only after constructing an appropriate response and repeating the relevant experiments. It does not directly validate or invalidate published real-cluster measurements.

\section{Conclusions}
A CMD population-inference workflow requires separate validation of statistical detection, latent-distribution reconstruction, and interval calibration. A tailored test can be sensitive to a declared feature even when state-by-state reconstruction remains difficult. Agreement in a global mean such as $\langle q\rangle$ does not establish local distribution recovery. Finite iteration affects both estimates and resampling intervals, and a nominal marginal interval label requires empirical repeated-sampling calibration. These conclusions are conditional on the finite response, grid, nuisance distribution, support, and selection model examined here.

\section{Data and code availability}
The existing public repository is \url{https://github.com/chihuanbin/binary_mass_ratio}. The final revision archive is prepared locally for deposit and is not represented as already published. Its release identifiers must be completed at deposit: \texttt{RELEASE\_TAG\_TO\_BE\_FILLED}, \texttt{RELEASE\_COMMIT\_TO\_BE\_FILLED}, \texttt{RELEASE\_DATE\_TO\_BE\_FILLED}, and \texttt{ZENODO\_DOI\_TO\_BE\_FILLED}. The archive contains the frozen response and historical summaries, declared truth vectors, new count samples, thresholds, convergence diagnostics, response draws, interval endpoints, generated figure/table sources, source hashes, environment specification, and executable reproduction commands. Historical count vectors and fitted state vectors were not retained. Rebuilding the raw-to-response calibration additionally requires the versioned MIST and cluster/member inputs documented in the inventory; this revision verifies and uses the archived response without rerunning that chain.

\software{Python, NumPy, SciPy, Matplotlib, Numba, pytest.}

\appendix
\section{Six-cluster master catalogue}\label{app:catalogue}
The parameter centers in Table~\ref{tab:catalogue} are read from the retained master-catalogue CSV. NGC 188 and NGC 6791 are master-catalogue entries but do not enter the calibration generator's selected mass-support subsets. This distinction precedes the system-level component-support rejection described in Section~2.
\input{tables/catalogue.tex}

\section{Projection of two specified contrasts}\label{app:projection}
The one-bin statistical-detection contrast and the distinct two-bin S0-to-S1 reconstruction contrast are projected deterministically. Since every frozen row has unit detected efficiency, $\Delta p=R_{\rm det}^{\mathsf T}\Delta t$; the signed changes sum to zero in both spaces. Figure~\ref{fig:projection} displays these two separate comparisons. They are not sequential stages or an identifiability test.
\begin{figure}[ht]
\centering\includegraphics[width=0.94\columnwidth]{figures/projection.pdf}
\caption{Two separate contrasts, each shown as latent and projected probability changes. The one-bin statistical-detection contrast uses the historical globally normalized alternative; the two-bin reconstruction contrast compares S0 and S1. Values come directly from the frozen response and specified truths. They are not sequential stages of a single experiment.}
\label{fig:projection}\end{figure}

\section{Proposal and first-failure accounting}\label{app:proposal}
The restricted replay contains 296,000 proposals: 50,386 fail the adopted interpolation-support rule and 245,614 are in support. All 245,614 pass the subsequent listed model checks, with zero first failures for interpolation domain, inverse fit, $G$ cut, uncertainty cut, or nonfinite photometry. The separately retained overlapping flags for $G$ and uncertainty cuts are also zero and are not additive first-failure counts. These outcomes describe the replay, not Gaia completeness (Figure~\ref{fig:accounting}). Event identifiers linking these proposals to the calibration quota of 200,016 rows were not retained. The arithmetic difference, 45,598, therefore remains unlinked and is not assigned to a loss or surplus category.
\begin{figure}[ht]
\centering\includegraphics[width=\columnwidth]{figures/accounting.pdf}
\caption{Proposal/support accounting (left) and first-failure outcomes conditional on support (right). Support rejects and in-support counts sum to 296,000. All 245,614 in-support events pass the subsequent listed checks in the replay. No quota stage is shown because calibration-row linkage is absent. These are simulation model outcomes, not measured catalogue completeness.}
\label{fig:accounting}\end{figure}

\section{Nonuniform reconstruction extension}\label{app:extended}
The three truths are analytically declared in Section~3 and retain uniform primary-mass marginals. Figure~\ref{fig:extended} compares maximum local error with the fixed-uniform baseline across all 500 realizations per truth. The complete error metrics, conditional mean-$q$ errors, counts, fitted vectors, and iteration-cap diagnostics are supplied in the archive. The synthetic high-$q$ and mass-dependent shapes are benchmark designs rather than inferred astrophysical populations.
\begin{figure}[ht]\centering\includegraphics[width=0.88\columnwidth]{figures/extended.pdf}
\caption{Nonuniform reconstruction extension at $N=20,000$. Points show the mean maximum local error over 500 realizations; bars span the 16th--84th percentiles of realization errors, not confidence intervals on the mean. Dashed segments show analytic fixed-uniform errors. The 200-iteration and 5,000-iteration estimates use identical counts and initialization; neither cap is claimed to achieve convergence. All realizations are included.}
\label{fig:extended}\end{figure}

\clearpage
\bibliographystyle{aasjournalv7}
\bibliography{references}
\end{document}
'''
 ext_sentence='For the three nonuniform truths, the mean maximum errors at 200 iterations are '+', '.join(f"{ext[p]['EM200']['E']:.5f}" for p in ['R1','R2','R3'])+' (R1--R3), compared with fixed-uniform errors '+', '.join(f"{ext[p]['uniform']['E']:.5f}" for p in ['R1','R2','R3'])+r'''. These shapes broaden the baseline comparison without constituting a population discovery. Figure~\ref{fig:extended} and the full metrics in Appendix~\ref{app:extended} show where the finite-iteration inversion adds information for these declared distributions.'''
 results=results.replace('EXTENDED_RESULTS',ext_sentence)
 abstract=r'''\begin{abstract}
Photometric binary mass-ratio structure requires validation of sensitivity, reconstruction accuracy, and uncertainty calibration as distinct empirical properties. We construct a controlled CMD benchmark with a discrete joint primary-mass--mass-ratio response and evaluate these three operating characteristics separately. For a specified historical one-cell departure at $N=20,000$, Pearson and tailored simple-versus-simple likelihood-ratio tests have powers 0.1025 and 0.6815 at their original thresholds; with independently calibrated common 5\% target size and primary-mass marginals held fixed, their powers are 0.0857 and 0.5030. Historical uniform and near-uniform reconstruction errors exceed a fixed-uniform baseline, while three declared nonuniform truths extend the comparison; new fits exhibit substantial iteration-cap sensitivity. Nominal 68\% marginal interval coverage is 84.0\%--85.8\% historically and 87.9\%--90.0\% under an independent empirical-zero-preserving response perturbation. Sensitivity to a specified departure, local latent recovery, and interval calibration therefore support different claims even within one frozen response. The results are conditional on the finite grid, nuisance sampling, support, selection model, and inverse fits supplied the true nuisance parameters.
\end{abstract}
\keywords{Binary stars --- Color-magnitude diagrams --- Stellar populations}

'''
 text=head+abstract+intro+benchmark+methods+results
 text=text.replace('truncated at zero','clipped at zero').replace('the iteration cap part of the estimand-defining procedure','the iteration cap an explicit part of the estimation procedure').replace('the published finite-step vectors are not authenticated converged maximum-likelihood solutions','the original finite-step outputs are not documented converged maximum-likelihood solutions')
 release_file=ROOT/'provenance/public_release.json'
 release=None
 if release_file.exists():
  release=json.loads(release_file.read_text())
  if release.get('doi') and release.get('commit'):
   availability=(r'The reproducibility and benchmark package is publicly available from the versioned GitHub release '+r'\href{'+release['github_release_url']+r'}{\texttt{'+release['tag'].replace('_',r'\_')+r'}} and is archived on Zenodo at \url{https://doi.org/'+release['doi']+r'} \citep{BenchmarkRelease}. The release is tied to commit \texttt{'+release['commit']+r'} and was published on '+release['release_date']+r'. The archive contains the frozen response, historical aggregate summaries, declared truth vectors, independent count samples, test thresholds, reconstruction and likelihood diagnostics, response draws, interval endpoints for the new experiments, generated figure/table sources, source hashes, environment specification, and executable reproduction commands. The earlier files \texttt{reproducibility\_package.zip} and \texttt{CMD Binary-Mass-Ratio Benchmark.zip} remain in the repository as historical packages. Historical count vectors, fitted state vectors, and interval endpoint arrays were not retained. Rebuilding the raw-to-response calibration additionally requires the versioned MIST and cluster/member inputs documented in the inventory; this release verifies and uses the archived response without rerunning that chain.')
   a=text.index('The existing public repository is ');b=text.index('\n',a);text=text[:a]+availability+text[b:]
 FINAL.mkdir(parents=True,exist_ok=True)
 (FINAL/'main.tex').write_text(text)
 for f in ['references.bib','aastex701.cls','aasjournalv7.bst']:shutil.copy2(OLD/'overleaf_upload/manuscript'/f,FINAL/f)
 bib=(FINAL/'references.bib').read_text();bib=bib.replace('{Alexander}, Jason S. and {Albrow}, Michael D. and others','{Alexander}, Jason S. and {Albrow}, Michael D.');bib=bib.replace('{Li}, Lu and {Shao}, Zhengyi and others','{Li}, Lu and {Shao}, Zhengyi and {Li}, Zhao-Zhou and {Yu}, Jincheng and {Zhong}, Jing and {Chen}, Li');bib=bib.replace('1050-1064','1050--1064')
 if release and release.get('doi') and release.get('commit'):
  bib+='\n@MISC{BenchmarkRelease,\n author = {{Chi}, Huanbin and {Wang}, Feng},\n title = {{CMD Binary-Mass-Ratio Benchmark: PASP final release}},\n year = {2026},\n publisher = {Zenodo},\n version = {'+release['tag']+'},\n doi = {'+release['doi']+'},\n url = {https://doi.org/'+release['doi']+'}\n}\n'
 (FINAL/'references.bib').write_text(bib)
 shutil.copytree(ROOT/'tables',FINAL/'tables',dirs_exist_ok=True);(FINAL/'figures').mkdir(exist_ok=True)
 for name in ['response','historical_power','reconstruction','projection','accounting','extended']:shutil.copy2(ROOT/'figures'/f'{name}.pdf',FINAL/'figures'/f'{name}.pdf')
 (ROOT/'reports/manuscript_changes.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),text.splitlines(True),fromfile=str(source),tofile=str(FINAL/'main.tex'))))
 savejson(ROOT/'provenance/manuscript_build.json',{'baseline_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'final_sha256':hashlib.sha256((FINAL/'main.tex').read_bytes()).hexdigest(),'title_decision':'retain original title: fixed-M1 experiment completed and still sensitive','authors':['Huanbin Chi','Feng Wang']})
 print(FINAL/'main.tex')
if __name__=='__main__':main()
