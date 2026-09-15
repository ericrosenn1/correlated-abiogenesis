# Correlated Abiogenesis

[![Tests](https://github.com/ericrosenn1/correlated-abiogenesis/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/ericrosenn1/correlated-abiogenesis/actions/workflows/tests.yml) [![Codecov](https://codecov.io/gh/ericrosenn1/correlated-abiogenesis/graph/badge.svg)](https://app.codecov.io/gh/ericrosenn1/correlated-abiogenesis) [![Python 3.10-3.12](https://img.shields.io/badge/Python-3.10--3.12-3776AB?logo=python&logoColor=white)](https://github.com/ericrosenn1/correlated-abiogenesis/blob/main/pyproject.toml) [![MIT License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/ericrosenn1/correlated-abiogenesis/blob/main/LICENSE) [![Release](https://img.shields.io/github/v/release/ericrosenn1/correlated-abiogenesis?display_name=tag&sort=semver)](https://github.com/ericrosenn1/correlated-abiogenesis/releases)

[![CITATION.cff](https://img.shields.io/badge/CITATION.cff-cite-555.svg)](https://github.com/ericrosenn1/correlated-abiogenesis/blob/main/CITATION.cff) [![Last commit](https://img.shields.io/github/last-commit/ericrosenn1/correlated-abiogenesis/main?label=last%20commit)](https://github.com/ericrosenn1/correlated-abiogenesis/commits/main) [![Contributors](https://img.shields.io/github/contributors/ericrosenn1/correlated-abiogenesis)](https://github.com/ericrosenn1/correlated-abiogenesis/graphs/contributors) [![Repository size](https://img.shields.io/github/repo-size/ericrosenn1/correlated-abiogenesis)](https://github.com/ericrosenn1/correlated-abiogenesis) [![CodeQL](https://github.com/ericrosenn1/correlated-abiogenesis/actions/workflows/codeql.yml/badge.svg?branch=main)](https://github.com/ericrosenn1/correlated-abiogenesis/actions/workflows/codeql.yml)

This repository contains the calculations and simulations supporting Appendix B and Figures 9–10 of the manuscript “Correlated Abiogenesis and the Sign of the Bayesian Update: Why a Second Origin of Life Need Not Imply an Abundant Galaxy.”

## Scope

The implementation covers:

1. the zero-background rare-basin between-model crossover calculation;
2. the known-centre estimator simulation for the decay scale xi; and
3. the matched Neyman–Scott / shot-noise Cox positional-degeneracy simulation.

The simulations are illustrative. The positional equivalence in Proposition 5 is analytical and is not inferred from visual agreement between simulated realizations.

## Requirements

Python 3.10–3.12 is supported. The recommended installation is an editable install with the development dependencies.

## Installation

Create and activate a virtual environment:

~~~powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
~~~

On macOS or Linux, activate with source .venv/bin/activate. A requirements-file installation is also available:

~~~text
python -m pip install -r requirements.txt
~~~

## Reproduce Appendix B

The authoritative manuscript-associated run is:

~~~text
python -m correlated_abiogenesis.appendix_b --seed 7 --replicates 10000 --output-dir outputs
~~~

Windows users may instead run:

~~~powershell
.\scripts\run_all.ps1 -Seed 7 -Replicates 10000 -DoNotOpenOutput
~~~

## Expected numerical checks

- p_I = 0.248685199098
- T_I = 0.285930539413
- crossover radius R_x: 11.843345 for xi = 1, 19.527808 for xi = 2, and 35.075159 for xi = 5

For Appendix B.2, with true xi = 2 and 10,000 Monte Carlo replicates:

| n | median xi_hat | empirical 95% interval width |
|---:|---:|---:|
| 1 | 1.769 | 4.379 |
| 2 | 1.894 | 3.105 |
| 3 | 1.923 | 2.642 |
| 5 | 1.962 | 2.030 |
| 10 | 1.978 | 1.412 |
| 30 | 1.994 | 0.811 |
| 100 | 1.998 | 0.454 |

The seeded Appendix B.3 propagation realization contains 834 latent parent origins and 3,625 observed points. The correlated-origin realization contains 3,773 independent-origin points.

## Output files

The reproduction command writes runtime artifacts to outputs/:

- figure_9_appendix_b.png and figure_9_appendix_b.pdf: the three-panel Appendix B summary;
- figure_10_matched_constructions.png and figure_10_matched_constructions.pdf: the matched positional constructions;
- appendix_b2_estimator_summary.csv: the Appendix B.2 estimator table;
- appendix_b1_crossover.csv and appendix_b3_pair_correlation.csv: crossover and pair-correlation results;
- appendix_b_parameters.json and appendix_b_results.md: run parameters and a readable numerical summary.

The manuscript-associated PNGs are committed in figures/, and the small reference tables and parameter records are committed in results/.

## Tests

~~~text
python -m pytest
~~~

The repository contains numerical and file-level regression checks. The current publication preparation run passes all 6 tests.

## Repository organization

- src/ contains the authoritative Python package.
- scripts/ contains the command-line and PowerShell wrappers.
- tests/ contains numerical and repository-file checks.
- figures/ contains the committed manuscript-associated reference images.
- results/ contains the committed reference numerical outputs.
- docs/ contains methods, reproduction instructions, and the figure map.

## Reproducibility

The manuscript-associated run uses master random seed 7 and 10,000 replicates for Appendix B.2. Package and environment specifications are included, together with the tested exact environment. Generated reference figures and numerical summaries are committed for comparison with reruns.

## Citation

Repository citation information is available in CITATION.cff.

## License

This software is released under the permissive MIT License. See LICENSE.
