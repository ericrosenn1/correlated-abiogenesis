# Contributing

Changes should preserve the mathematical definitions used in the manuscript and include tests for altered calculations. Before opening a pull request, run:

~~~text
python -m pip install -e ".[dev]"
python -m compileall -q src scripts tests
python -m pytest
python -m correlated_abiogenesis.appendix_b --quick --output-dir outputs-smoke
~~~

Commit generated figures or tables only when they are deliberate manuscript-associated reference outputs.
