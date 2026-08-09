# Computational methods

## B.1 Between-model predictive crossover

The homogeneous model uses a Gamma prior with shape `alpha` and rate `beta`. After observing `n` independent origins in searched volume `V`, the predictive probability of at least one independent origin in target volume `|F|` is

```text
p_I = 1 - [(beta + V) / (beta + V + |F|)]^(alpha + n).
```

The probability-equivalent Poisson mean is

```text
T_I = -log(1 - p_I).
```

Under the zero-background rare-basin model and the small-target approximation,

```text
delta_R = [mu |F| / (8 pi xi^3)] exp(-R / xi),
p_C(R) = 1 - exp(-delta_R).
```

The between-model crossover is the solution to `delta_R = T_I`:

```text
R_x = xi log{mu |F| / [8 pi xi^3 T_I]}.
```

The default illustrative parameters are stored in `CrossoverParameters` and written to the run-specific JSON record.

## B.2 Known-centre decay-scale estimator

For the normalized kernel

```text
g_xi(r) = exp(-r / xi) / (8 pi xi^3),
```

the radial density in three dimensions is

```text
f_R(r) = r^2 exp(-r / xi) / (2 xi^3),
```

which is Gamma with shape `3` and scale `xi`. For `n` independent radii and a known centre,

```text
xi_hat = (1 / 3n) sum_i r_i = mean(r_i) / 3.
```

The simulation draws the radial sum directly from Gamma with shape `3n` and scale `xi`, then divides by `3n`. The reported interval is the empirical interval between the 2.5th and 97.5th percentiles.

## B.3 Matched Neyman-Scott and shot-noise Cox constructions

The two-dimensional simulation is used for visualization. Parent centres form a homogeneous Poisson process with intensity `kappa` on a periodic square of side length `L`. Each parent produces a Poisson number of observed points with mean `alpha`. Offspring displacements are isotropic Gaussian with standard deviation `sigma` in each coordinate and are wrapped periodically.

Conditional on the same centre process, the observed point pattern is an inhomogeneous Poisson process with intensity

```text
Lambda(x) = alpha sum_j k_sigma(x - c_j).
```

This is both the conditional offspring process of the Neyman-Scott construction and the conditional point process directed by the associated shot-noise Cox intensity. The equality of the unmarked positional laws is analytical; the simulations are illustrative.

The pair-correlation estimate uses minimum-image pair distances. For annulus area `A_r` in a window of area `A`, the homogeneous-Poisson expected unordered pair count is

```text
[n(n - 1) / 2] A_r / A.
```

The observed annulus pair count divided by this expectation estimates `g(r)`. The theoretical Thomas-process pair correlation is

```text
g(r) = 1 + exp[-r^2 / (4 sigma^2)] / (4 pi kappa sigma^2).
```
