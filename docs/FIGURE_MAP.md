# Appendix B figure and result map

| Manuscript item | Implementation | Reproduction command | Committed reference |
|---|---|---|---|
| Figure 9 | run_appendix_b and make_figure_9 in src/correlated_abiogenesis/appendix_b.py | python -m correlated_abiogenesis.appendix_b --seed 7 --replicates 10000 --output-dir outputs | figures/figure_9_appendix_b.png |
| Figure 10 | run_appendix_b and make_figure_10 in src/correlated_abiogenesis/appendix_b.py | python -m correlated_abiogenesis.appendix_b --seed 7 --replicates 10000 --output-dir outputs | figures/figure_10_matched_constructions.png |
| Appendix B.2 estimator table | simulate_known_center_estimator and summarize_estimator_simulation in src/correlated_abiogenesis/appendix_b.py | same command | results/appendix_b2_estimator_summary.csv |

The same run also produces the B.1 crossover table and the B.3 pair-correlation table. Their committed references are results/appendix_b1_crossover.csv and results/appendix_b3_pair_correlation.csv.

The simulations provide computational illustrations. Proposition 5 establishes the positional equivalence analytically.
