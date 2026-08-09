# Appendix B computational results

Master random seed: `7`

## B.1 Between-model crossover

- Homogeneous predictive probability: `0.248685199098`
- Probability-equivalent Poisson mean: `0.285930539413`

| xi | crossover radius |
|---:|---:|
| 1 | 11.843345497280 |
| 2 | 19.527807911201 |
| 5 | 35.075158799890 |

## B.2 Known-centre estimator

| n | median xi_hat | 2.5% quantile | 97.5% quantile | interval width |
|---:|---:|---:|---:|---:|
| 1 | 1.768722 | 0.401092 | 4.779715 | 4.378624 |
| 2 | 1.894207 | 0.748119 | 3.852960 | 3.104842 |
| 3 | 1.923197 | 0.924736 | 3.566874 | 2.642139 |
| 5 | 1.961948 | 1.115050 | 3.144856 | 2.029806 |
| 10 | 1.977752 | 1.357777 | 2.769688 | 1.411912 |
| 30 | 1.994431 | 1.616430 | 2.427363 | 0.810933 |
| 100 | 1.998356 | 1.781127 | 2.234718 | 0.453590 |

## B.3 Matched cluster constructions

- Window side length: `60.0`
- Parent intensity: `0.23222222222222222`
- Mean offspring count: `4.4`
- Gaussian dispersal sigma: `0.65`
- Propagation realization: `834` latent parents and `3625` observed points
- Correlated-origin realization: `3773` independent-origin points

The two observed point patterns are independent realizations of the same matched Thomas/Neyman-Scott positional law. Their finite-sample pair-correlation estimates are not expected to coincide pointwise.
