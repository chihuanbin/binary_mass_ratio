# CMD Binary-Mass-Ratio Benchmark

Authors: Huanbin Chi and Feng Wang.

A controlled simulation benchmark for statistical detection of specified population departures, finite-iteration latent-distribution reconstruction, and empirical interval calibration. These operating characteristics are validated separately. The benchmark does not validate Gaia completeness or arbitrary real-cluster distributions.

The original `reproducibility_package.zip` and `CMD Binary-Mass-Ratio Benchmark.zip` are retained as historical packages. The complete PASP final source/result archive is distributed through versioned GitHub releases. Analysis and audit scripts are versioned under `revision_pasp_final/`; execute them from the complete release archive, which also contains their frozen inputs and independently generated count/response arrays.

The final release patch authenticates four historical width aggregates and replays 41,540 retained independent fits for likelihood monotonicity, without new simulation draws. All likelihoods are non-decreasing; the fit stopping tolerance is not satisfied at the recorded iteration caps.
