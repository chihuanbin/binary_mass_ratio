# Reproducing the final PASP benchmark

Use the recorded environment in `environment.yml` and `provenance/environment.json`. Commands below run from the workspace/release root (which contains both `revision_pasp_final` and the frozen `pasp_cmd_recoverability_revision` source tree). `python` denotes the recorded Python environment. TeX Live/latexmk and Poppler are needed for compilation/rendering.

```sh
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/pasp_mpl
python revision_pasp_final/scripts/verify_baseline.py
python revision_pasp_final/scripts/run_detection.py
python revision_pasp_final/scripts/run_convergence_response.py
python revision_pasp_final/scripts/run_extended.py
python revision_pasp_final/scripts/build_tables.py
python revision_pasp_final/scripts/build_figures.py
python revision_pasp_final/scripts/build_manuscript.py
python revision_pasp_final/scripts/write_audits.py
python revision_pasp_final/scripts/audit_numerical_text.py
python -m pytest revision_pasp_final/tests -q
latexmk -cd -pdf -interaction=nonstopmode -halt-on-error manuscript/pasp_submission_final/main.tex
python revision_pasp_final/scripts/render_pages.py
```

The combined detection script performs both fixed-M1 and common-target-size experiments. The combined convergence/response script performs all cap reruns and empirical-zero-preserving interval experiments. Historical count vectors, fitted vectors and bootstrap endpoints are absent; those fits cannot be authenticated or recreated as historical values. Newly seeded outputs are independent extensions.

Reference metadata are cached. Refreshing them separately requires network access:

```sh
python revision_pasp_final/scripts/audit_references.py
```

All scientific experiments use the same frozen `calibration_200k.npz`. Seeds are saved in output/configuration files and PCG64 is used. Numba is restricted to four threads; floating differences across platforms are checked with explicit tolerances. Reproduction scripts write only new namespaces; baseline verification compares original hashes before analysis. PDF creation dates may differ; exact scientific array and summary reproduction is the required criterion.

`python revision_pasp_final/scripts/build_release.py` creates the portable source/result release and clean Overleaf archive. Raw MIST and cluster/member tables are available in the original local project but excluded from this portable release; source inventory includes their hashes. Rebuilding the raw-to-response chain is outside these commands.

## Deposit status

Public repository verified: https://github.com/chihuanbin/binary_mass_ratio
Current remote commit inspected: `10c907f21bb9f43ef385372035eecddf45761cce` (2026-09-28). This is the existing remote snapshot, NOT the final revision release. No tagged GitHub release was found during this audit.

Fill only after publishing the revision archive:

- `RELEASE_TAG_TO_BE_FILLED`
- `RELEASE_COMMIT_TO_BE_FILLED`
- `RELEASE_DATE_TO_BE_FILLED`
- `ZENODO_DOI_TO_BE_FILLED`

No remote publication was performed. The final PDF retains explicit placeholders rather than claiming that unpublished files are public.
