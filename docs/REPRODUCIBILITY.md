# Reproducibility

## Environment setup

~~~text
python -m venv .venv
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
~~~

On Windows, activate the environment with .\.venv\Scripts\Activate.ps1. On macOS or Linux, use source .venv/bin/activate.

## Authoritative manuscript run

~~~text
python -m correlated_abiogenesis.appendix_b --seed 7 --replicates 10000 --output-dir outputs
~~~

The equivalent Windows wrapper is:

~~~powershell
.\scripts\run_all.ps1 -Seed 7 -Replicates 10000 -DoNotOpenOutput
~~~

The run uses master seed 7. NumPy SeedSequence creates independent streams for the estimator, propagation, and correlated-origin realizations. Appendix B.2 uses 10,000 Monte Carlo replicates for every sample size.

## Expected outputs

The command creates:

- figure_9_appendix_b.png and .pdf;
- figure_10_matched_constructions.png and .pdf;
- appendix_b1_crossover.csv;
- appendix_b2_estimator_summary.csv;
- appendix_b3_pair_correlation.csv;
- appendix_b_parameters.json; and
- appendix_b_results.md.

Reference PNGs are stored in figures/. Reference tables and run metadata are stored in results/.

## Numerical validation

Expected principal values:

- p_I = 0.248685199098
- T_I = 0.285930539413
- R_x(xi = 1) = 11.843345
- R_x(xi = 2) = 19.527808
- R_x(xi = 5) = 35.075159

The Appendix B.2 medians and empirical 95% interval widths are:

| n | median xi_hat | interval width |
|---:|---:|---:|
| 1 | 1.769 | 4.379 |
| 2 | 1.894 | 3.105 |
| 3 | 1.923 | 2.642 |
| 5 | 1.962 | 2.030 |
| 10 | 1.978 | 1.412 |
| 30 | 1.994 | 0.811 |
| 100 | 1.998 | 0.454 |

The seeded Appendix B.3 propagation realization contains 834 latent parents and 3,625 observed points. The correlated-origin realization contains 3,773 independent-origin points.

## Validation commands

~~~text
python -m compileall -q src scripts tests
python -m pytest
python -m correlated_abiogenesis.appendix_b --quick --output-dir outputs-smoke
~~~

The full manuscript command, rather than the smoke run, must be used to reproduce the committed reference artifacts.

## Platform notes

The tested environment is recorded in SOFTWARE_ENVIRONMENT.md and requirements-lock.txt. Identical versions and seed reproduce the numerical tables and PNGs byte-for-byte in the validated Windows environment. PDF files can contain time-dependent metadata even when their rendered contents are unchanged.
