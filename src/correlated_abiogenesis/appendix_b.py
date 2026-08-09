"""Reproduce the calculations and simulations reported in Appendix B.

The module implements three independent components:

1. The zero-background rare-basin crossover calculation.
2. The known-centre estimator simulation for the decay scale xi.
3. The matched Neyman-Scott / shot-noise Cox illustration.

Run as a module:

    python -m correlated_abiogenesis.appendix_b --output-dir outputs
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy import __version__ as scipy_version
from scipy.spatial import cKDTree


MASTER_SEED = 7


@dataclass(frozen=True)
class CrossoverParameters:
    """Parameters for the Appendix B.1 crossover calculation."""

    alpha: float = 2.0
    n_origins: int = 1
    beta: float = 3.0
    searched_volume: float = 7.0
    target_volume: float = 1.0
    basin_mass: float = 1_000_000.0
    xi_values: tuple[float, ...] = (1.0, 2.0, 5.0)
    radius_min: float = 0.0
    radius_max: float = 38.0
    radius_points: int = 700


@dataclass(frozen=True)
class EstimatorParameters:
    """Parameters for the Appendix B.2 known-centre simulation."""

    true_xi: float = 2.0
    sample_sizes: tuple[int, ...] = (1, 2, 3, 5, 10, 30, 100)
    replicates: int = 10_000
    lower_quantile: float = 0.025
    upper_quantile: float = 0.975


@dataclass(frozen=True)
class ClusterParameters:
    """Parameters for the Appendix B.3 matched cluster-process simulation."""

    window_size: float = 60.0
    parent_intensity: float = 836.0 / 3600.0
    mean_offspring: float = 4.4
    dispersal_sigma: float = 0.65
    maximum_pair_distance: float = 12.0
    pair_distance_bins: int = 30


@dataclass(frozen=True)
class ClusterRealization:
    """One realization of a periodic Thomas/Neyman-Scott process."""

    parents: NDArray[np.float64]
    observed_points: NDArray[np.float64]
    offspring_counts: NDArray[np.int64]


@dataclass(frozen=True)
class PairCorrelationEstimate:
    """Estimated pair-correlation function on a periodic square window."""

    radius: NDArray[np.float64]
    estimate: NDArray[np.float64]
    pair_counts: NDArray[np.int64]
    poisson_expected_counts: NDArray[np.float64]


def gamma_poisson_predictive_probability(
    alpha: float,
    n_origins: int,
    beta: float,
    searched_volume: float,
    target_volume: float,
) -> float:
    """Return the Gamma-Poisson posterior predictive probability of >=1 event.

    Gamma(alpha, beta) uses shape ``alpha`` and rate ``beta``.
    """

    if alpha <= 0 or beta <= 0:
        raise ValueError("alpha and beta must be positive")
    if n_origins < 0:
        raise ValueError("n_origins must be non-negative")
    if searched_volume < 0 or target_volume < 0:
        raise ValueError("volumes must be non-negative")

    posterior_shape = alpha + n_origins
    posterior_rate = beta + searched_volume
    void_probability = (posterior_rate / (posterior_rate + target_volume)) ** posterior_shape
    return float(1.0 - void_probability)


def probability_equivalent_poisson_mean(probability: float) -> float:
    """Map a probability in [0, 1) to the equivalent Poisson mean."""

    if not 0.0 <= probability < 1.0:
        raise ValueError("probability must lie in [0, 1)")
    return float(-math.log1p(-probability))


def rare_basin_integrated_mean(
    radius: NDArray[np.float64] | float,
    basin_mass: float,
    target_volume: float,
    xi: float,
) -> NDArray[np.float64] | float:
    """Small-target approximation to the basin-contributed Poisson mean."""

    if basin_mass <= 0 or target_volume <= 0 or xi <= 0:
        raise ValueError("basin_mass, target_volume, and xi must be positive")

    prefactor = basin_mass * target_volume / (8.0 * math.pi * xi**3)
    radius_array = np.asarray(radius, dtype=float)
    if np.any(radius_array < 0):
        raise ValueError("radius must be non-negative")
    result = prefactor * np.exp(-radius_array / xi)
    if np.ndim(radius) == 0:
        return float(result)
    return result


def rare_basin_predictive_probability(integrated_mean: NDArray[np.float64] | float):
    """Return 1 - exp(-delta) for a non-negative integrated mean delta."""

    delta = np.asarray(integrated_mean, dtype=float)
    if np.any(delta < 0):
        raise ValueError("integrated_mean must be non-negative")
    probability = -np.expm1(-delta)
    if np.ndim(integrated_mean) == 0:
        return float(probability)
    return probability


def crossover_radius(
    basin_mass: float,
    target_volume: float,
    xi: float,
    homogeneous_probability: float,
) -> float:
    """Return the small-target between-model crossover radius.

    The returned value may be non-positive. A non-positive value means that the
    formal crossover lies at or outside the physical boundary R >= 0.
    """

    threshold = probability_equivalent_poisson_mean(homogeneous_probability)
    amplitude = basin_mass * target_volume / (8.0 * math.pi * xi**3)
    if amplitude <= 0 or threshold <= 0:
        raise ValueError("amplitude and probability-equivalent threshold must be positive")
    return float(xi * math.log(amplitude / threshold))


def simulate_known_center_estimator(
    rng: np.random.Generator,
    parameters: EstimatorParameters,
) -> dict[int, NDArray[np.float64]]:
    """Simulate the known-centre maximum-likelihood estimator of xi.

    For the normalized three-dimensional exponential kernel,

        g_xi(r) = exp(-r/xi) / (8*pi*xi^3),

    the radial distance has Gamma(shape=3, scale=xi) distribution. For n
    independent points and a known centre, the maximum-likelihood estimator is

        xi_hat = mean(r_i) / 3.

    The sum of n radial distances is Gamma(shape=3*n, scale=xi), which permits
    direct and numerically efficient simulation of xi_hat.
    """

    if parameters.true_xi <= 0:
        raise ValueError("true_xi must be positive")
    if parameters.replicates <= 0:
        raise ValueError("replicates must be positive")
    if not 0 <= parameters.lower_quantile < parameters.upper_quantile <= 1:
        raise ValueError("quantiles must satisfy 0 <= lower < upper <= 1")

    estimates: dict[int, NDArray[np.float64]] = {}
    for n in parameters.sample_sizes:
        if n <= 0:
            raise ValueError("all sample sizes must be positive")
        radial_sums = rng.gamma(
            shape=3.0 * n,
            scale=parameters.true_xi,
            size=parameters.replicates,
        )
        estimates[n] = radial_sums / (3.0 * n)
    return estimates


def summarize_estimator_simulation(
    estimates: dict[int, NDArray[np.float64]],
    parameters: EstimatorParameters,
) -> list[dict[str, float | int]]:
    """Summarize the Monte Carlo estimator distributions."""

    rows: list[dict[str, float | int]] = []
    for n in parameters.sample_sizes:
        values = estimates[n]
        q_low, median, q_high = np.quantile(
            values,
            [parameters.lower_quantile, 0.5, parameters.upper_quantile],
        )
        rows.append(
            {
                "n": n,
                "median_xi_hat": float(median),
                "q_lower": float(q_low),
                "q_upper": float(q_high),
                "interval_width": float(q_high - q_low),
                "mean_xi_hat": float(np.mean(values)),
                "standard_deviation": float(np.std(values, ddof=1)),
            }
        )
    return rows


def simulate_periodic_thomas_process(
    rng: np.random.Generator,
    parameters: ClusterParameters,
) -> ClusterRealization:
    """Simulate a Thomas/Neyman-Scott process on a periodic square window.

    Parent locations form a homogeneous Poisson process. Each parent generates a
    Poisson number of offspring with isotropic Gaussian displacement. Coordinates
    are wrapped onto the square, eliminating boundary loss and making the pair-
    correlation calculation use a transparent toroidal edge treatment.
    """

    if parameters.window_size <= 0:
        raise ValueError("window_size must be positive")
    if parameters.parent_intensity <= 0:
        raise ValueError("parent_intensity must be positive")
    if parameters.mean_offspring <= 0:
        raise ValueError("mean_offspring must be positive")
    if parameters.dispersal_sigma <= 0:
        raise ValueError("dispersal_sigma must be positive")

    area = parameters.window_size**2
    parent_count = int(rng.poisson(parameters.parent_intensity * area))
    if parent_count == 0:
        raise RuntimeError("simulation generated zero parents; choose a larger expected count")

    parents = rng.uniform(0.0, parameters.window_size, size=(parent_count, 2))
    offspring_counts = rng.poisson(parameters.mean_offspring, size=parent_count).astype(np.int64)
    total_offspring = int(offspring_counts.sum())
    if total_offspring == 0:
        raise RuntimeError("simulation generated zero offspring")

    parent_index = np.repeat(np.arange(parent_count), offspring_counts)
    displacements = rng.normal(
        loc=0.0,
        scale=parameters.dispersal_sigma,
        size=(total_offspring, 2),
    )
    observed_points = (parents[parent_index] + displacements) % parameters.window_size

    return ClusterRealization(
        parents=np.asarray(parents, dtype=float),
        observed_points=np.asarray(observed_points, dtype=float),
        offspring_counts=offspring_counts,
    )


def estimate_pair_correlation_periodic(
    points: NDArray[np.float64],
    window_size: float,
    bin_edges: NDArray[np.float64],
) -> PairCorrelationEstimate:
    """Estimate g(r) with a periodic-window annulus-count estimator.

    Pair distances use the minimum-image convention. Under complete spatial
    randomness, the expected unordered pair count in an annulus is

        n(n-1)/2 * annulus_area / window_area.

    Dividing observed counts by this expectation yields an estimator centered on
    one for a homogeneous Poisson process.
    """

    points = np.asarray(points, dtype=float)
    bin_edges = np.asarray(bin_edges, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (n, 2)")
    if len(points) < 2:
        raise ValueError("at least two points are required")
    if window_size <= 0:
        raise ValueError("window_size must be positive")
    if bin_edges.ndim != 1 or len(bin_edges) < 2 or np.any(np.diff(bin_edges) <= 0):
        raise ValueError("bin_edges must be a strictly increasing one-dimensional array")
    if bin_edges[0] < 0 or bin_edges[-1] > window_size / 2:
        raise ValueError("pair-distance bins must lie within [0, window_size/2]")

    wrapped = points % window_size
    tree = cKDTree(wrapped, boxsize=window_size)
    pairs = tree.query_pairs(float(bin_edges[-1]), output_type="ndarray")

    if pairs.size == 0:
        pair_distances = np.empty(0, dtype=float)
    else:
        differences = np.abs(wrapped[pairs[:, 0]] - wrapped[pairs[:, 1]])
        minimum_image = np.minimum(differences, window_size - differences)
        pair_distances = np.sqrt(np.sum(minimum_image**2, axis=1))

    pair_counts, _ = np.histogram(pair_distances, bins=bin_edges)
    annulus_area = math.pi * (bin_edges[1:] ** 2 - bin_edges[:-1] ** 2)
    window_area = window_size**2
    unordered_pairs = len(points) * (len(points) - 1) / 2.0
    expected_counts = unordered_pairs * annulus_area / window_area
    estimate = np.divide(
        pair_counts,
        expected_counts,
        out=np.full_like(expected_counts, np.nan, dtype=float),
        where=expected_counts > 0,
    )
    radius = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    return PairCorrelationEstimate(
        radius=radius,
        estimate=estimate,
        pair_counts=pair_counts.astype(np.int64),
        poisson_expected_counts=expected_counts,
    )


def theoretical_thomas_pair_correlation(
    radius: NDArray[np.float64],
    parent_intensity: float,
    dispersal_sigma: float,
) -> NDArray[np.float64]:
    """Return the pair-correlation function of a two-dimensional Thomas process."""

    radius = np.asarray(radius, dtype=float)
    if parent_intensity <= 0 or dispersal_sigma <= 0:
        raise ValueError("parent_intensity and dispersal_sigma must be positive")
    coefficient = 1.0 / (4.0 * math.pi * parent_intensity * dispersal_sigma**2)
    return 1.0 + coefficient * np.exp(-(radius**2) / (4.0 * dispersal_sigma**2))


def _write_csv(path: Path, rows: Iterable[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _save_figure(fig: plt.Figure, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_figure_9(
    output_dir: Path,
    crossover_parameters: CrossoverParameters,
    estimator_parameters: EstimatorParameters,
    estimator_estimates: dict[int, NDArray[np.float64]],
    propagation_pair_correlation: PairCorrelationEstimate,
    correlated_pair_correlation: PairCorrelationEstimate,
    homogeneous_probability: float,
    cluster_parameters: ClusterParameters,
) -> None:
    """Generate the three-panel Appendix Figure 9."""

    threshold = probability_equivalent_poisson_mean(homogeneous_probability)
    radius = np.linspace(
        crossover_parameters.radius_min,
        crossover_parameters.radius_max,
        crossover_parameters.radius_points,
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.6))

    ax = axes[0]
    for xi in crossover_parameters.xi_values:
        delta = rare_basin_integrated_mean(
            radius,
            crossover_parameters.basin_mass,
            crossover_parameters.target_volume,
            xi,
        )
        r_cross = crossover_radius(
            crossover_parameters.basin_mass,
            crossover_parameters.target_volume,
            xi,
            homogeneous_probability,
        )
        ax.plot(radius, delta, linewidth=1.8, label=fr"$\xi={xi:g}$")
        if 0.0 <= r_cross <= crossover_parameters.radius_max:
            ax.axvline(r_cross, linestyle=":", linewidth=1.0, alpha=0.8)
    ax.axhline(threshold, linestyle="--", linewidth=1.2, label=fr"$-\log(1-p_I)={threshold:.3f}$")
    ax.set_yscale("log")
    ax.set_xlabel(r"Distance from basin centre, $R$")
    ax.set_ylabel(r"Integrated expected count, $\delta_R$")
    ax.set_title("(A) Between-model predictive crossover")
    ax.text(
        0.02,
        0.04,
        fr"$\mu={crossover_parameters.basin_mass:.0e},\ |F|={crossover_parameters.target_volume:g}$",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
    )
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.25)

    ax = axes[1]
    box_data = [estimator_estimates[n] for n in estimator_parameters.sample_sizes]
    ax.boxplot(
        box_data,
        tick_labels=[str(n) for n in estimator_parameters.sample_sizes],
        whis=(
            100.0 * estimator_parameters.lower_quantile,
            100.0 * estimator_parameters.upper_quantile,
        ),
        showfliers=False,
        widths=0.62,
    )
    ax.axhline(
        estimator_parameters.true_xi,
        linestyle="--",
        linewidth=1.2,
        label=fr"True $\xi={estimator_parameters.true_xi:g}$",
    )
    ax.set_xlabel(r"Number of independent origins, $n$")
    ax.set_ylabel(r"Known-centre estimator, $\widehat{\xi}$")
    ax.set_title("(B) Precision improves continuously")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", alpha=0.25)

    ax = axes[2]
    ax.plot(
        propagation_pair_correlation.radius,
        propagation_pair_correlation.estimate,
        marker="o",
        markersize=3,
        linewidth=1.2,
        label="Propagation interpretation",
    )
    ax.plot(
        correlated_pair_correlation.radius,
        correlated_pair_correlation.estimate,
        marker="s",
        markersize=3,
        linewidth=1.2,
        label="Correlated-origin interpretation",
    )
    theory = theoretical_thomas_pair_correlation(
        propagation_pair_correlation.radius,
        cluster_parameters.parent_intensity,
        cluster_parameters.dispersal_sigma,
    )
    ax.plot(
        propagation_pair_correlation.radius,
        theory,
        linestyle="--",
        linewidth=1.2,
        label="Matched-process expectation",
    )
    ax.axhline(1.0, linestyle=":", linewidth=1.0, label="Poisson reference")
    ax.set_xlabel(r"Separation, $r$")
    ax.set_ylabel(r"Pair correlation, $\widehat{g}(r)$")
    ax.set_title("(C) Same clustering law, different interpretation")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(alpha=0.25)

    fig.tight_layout()
    _save_figure(fig, output_dir / "figure_9_appendix_b")


def make_figure_10(
    output_dir: Path,
    parameters: ClusterParameters,
    propagation: ClusterRealization,
    correlated: ClusterRealization,
    propagation_pair_correlation: PairCorrelationEstimate,
    correlated_pair_correlation: PairCorrelationEstimate,
) -> None:
    """Generate the detailed matched-construction illustration in Figure 10."""

    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.9))

    ax = axes[0]
    ax.scatter(
        propagation.observed_points[:, 0],
        propagation.observed_points[:, 1],
        s=2.2,
        alpha=0.65,
        linewidths=0,
    )
    ax.set_xlim(0, parameters.window_size)
    ax.set_ylim(0, parameters.window_size)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        "(A) Propagation interpretation\n"
        f"{len(propagation.parents):,} latent parent origins; "
        f"{len(propagation.observed_points):,} observed worlds"
    )
    ax.grid(False)

    ax = axes[1]
    ax.scatter(
        correlated.observed_points[:, 0],
        correlated.observed_points[:, 1],
        s=2.2,
        alpha=0.65,
        linewidths=0,
    )
    ax.set_xlim(0, parameters.window_size)
    ax.set_ylim(0, parameters.window_size)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        "(B) Correlated-origin interpretation\n"
        f"{len(correlated.observed_points):,} independent origins; no transport"
    )
    ax.grid(False)

    ax = axes[2]
    ax.plot(
        propagation_pair_correlation.radius,
        propagation_pair_correlation.estimate,
        marker="o",
        markersize=3,
        linewidth=1.3,
        label="Propagation realization",
    )
    ax.plot(
        correlated_pair_correlation.radius,
        correlated_pair_correlation.estimate,
        marker="s",
        markersize=3,
        linewidth=1.3,
        label="Correlated-origin realization",
    )
    theory = theoretical_thomas_pair_correlation(
        propagation_pair_correlation.radius,
        parameters.parent_intensity,
        parameters.dispersal_sigma,
    )
    ax.plot(
        propagation_pair_correlation.radius,
        theory,
        linestyle="--",
        linewidth=1.3,
        label="Matched-process expectation",
    )
    ax.axhline(1.0, linestyle=":", linewidth=1.0, label="Poisson reference")
    ax.set_xlabel(r"Separation, $r$")
    ax.set_ylabel(r"Pair correlation, $\widehat{g}(r)$")
    ax.set_title("(C) Unmarked positions do not identify the mechanism")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(alpha=0.25)

    fig.suptitle(
        "Matched Neyman-Scott and shot-noise Cox constructions",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    _save_figure(fig, output_dir / "figure_10_matched_constructions")


def _write_results_markdown(
    path: Path,
    homogeneous_probability: float,
    threshold: float,
    crossover_rows: list[dict[str, float]],
    estimator_rows: list[dict[str, float | int]],
    propagation: ClusterRealization,
    correlated: ClusterRealization,
    cluster_parameters: ClusterParameters,
    master_seed: int,
) -> None:
    lines = [
        "# Appendix B computational results",
        "",
        f"Master random seed: `{master_seed}`",
        "",
        "## B.1 Between-model crossover",
        "",
        f"- Homogeneous predictive probability: `{homogeneous_probability:.12f}`",
        f"- Probability-equivalent Poisson mean: `{threshold:.12f}`",
        "",
        "| xi | crossover radius |",
        "|---:|---:|",
    ]
    for row in crossover_rows:
        lines.append(f"| {row['xi']:.6g} | {row['crossover_radius']:.12f} |")

    lines.extend(
        [
            "",
            "## B.2 Known-centre estimator",
            "",
            "| n | median xi_hat | 2.5% quantile | 97.5% quantile | interval width |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in estimator_rows:
        lines.append(
            "| {n:d} | {median_xi_hat:.6f} | {q_lower:.6f} | "
            "{q_upper:.6f} | {interval_width:.6f} |".format(**row)
        )

    lines.extend(
        [
            "",
            "## B.3 Matched cluster constructions",
            "",
            f"- Window side length: `{cluster_parameters.window_size}`",
            f"- Parent intensity: `{cluster_parameters.parent_intensity}`",
            f"- Mean offspring count: `{cluster_parameters.mean_offspring}`",
            f"- Gaussian dispersal sigma: `{cluster_parameters.dispersal_sigma}`",
            f"- Propagation realization: `{len(propagation.parents)}` latent parents and "
            f"`{len(propagation.observed_points)}` observed points",
            f"- Correlated-origin realization: `{len(correlated.observed_points)}` independent-origin points",
            "",
            "The two observed point patterns are independent realizations of the same matched "
            "Thomas/Neyman-Scott positional law. Their finite-sample pair-correlation estimates "
            "are not expected to coincide pointwise.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def run_appendix_b(
    output_dir: Path,
    master_seed: int = MASTER_SEED,
    estimator_replicates: int | None = None,
) -> dict[str, object]:
    """Run all Appendix B calculations, simulations, tables, and figures."""

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    crossover_parameters = CrossoverParameters()
    estimator_parameters = EstimatorParameters(
        replicates=(
            estimator_replicates
            if estimator_replicates is not None
            else EstimatorParameters().replicates
        )
    )
    cluster_parameters = ClusterParameters()

    seed_sequence = np.random.SeedSequence(master_seed)
    estimator_seed, propagation_seed, correlated_seed = seed_sequence.spawn(3)
    estimator_rng = np.random.default_rng(estimator_seed)
    propagation_rng = np.random.default_rng(propagation_seed)
    correlated_rng = np.random.default_rng(correlated_seed)

    homogeneous_probability = gamma_poisson_predictive_probability(
        crossover_parameters.alpha,
        crossover_parameters.n_origins,
        crossover_parameters.beta,
        crossover_parameters.searched_volume,
        crossover_parameters.target_volume,
    )
    threshold = probability_equivalent_poisson_mean(homogeneous_probability)

    crossover_rows: list[dict[str, float]] = []
    for xi in crossover_parameters.xi_values:
        r_cross = crossover_radius(
            crossover_parameters.basin_mass,
            crossover_parameters.target_volume,
            xi,
            homogeneous_probability,
        )
        crossover_rows.append(
            {
                "xi": float(xi),
                "crossover_radius": float(r_cross),
                "integrated_mean_at_crossover": float(
                    rare_basin_integrated_mean(
                        r_cross,
                        crossover_parameters.basin_mass,
                        crossover_parameters.target_volume,
                        xi,
                    )
                ),
            }
        )

    estimator_estimates = simulate_known_center_estimator(estimator_rng, estimator_parameters)
    estimator_rows = summarize_estimator_simulation(estimator_estimates, estimator_parameters)

    propagation = simulate_periodic_thomas_process(propagation_rng, cluster_parameters)
    correlated = simulate_periodic_thomas_process(correlated_rng, cluster_parameters)
    bin_edges = np.linspace(
        0.0,
        cluster_parameters.maximum_pair_distance,
        cluster_parameters.pair_distance_bins + 1,
    )
    propagation_pair_correlation = estimate_pair_correlation_periodic(
        propagation.observed_points,
        cluster_parameters.window_size,
        bin_edges,
    )
    correlated_pair_correlation = estimate_pair_correlation_periodic(
        correlated.observed_points,
        cluster_parameters.window_size,
        bin_edges,
    )

    _write_csv(
        output_dir / "appendix_b1_crossover.csv",
        crossover_rows,
        ["xi", "crossover_radius", "integrated_mean_at_crossover"],
    )
    _write_csv(
        output_dir / "appendix_b2_estimator_summary.csv",
        estimator_rows,
        [
            "n",
            "median_xi_hat",
            "q_lower",
            "q_upper",
            "interval_width",
            "mean_xi_hat",
            "standard_deviation",
        ],
    )

    pair_rows: list[dict[str, float | int]] = []
    theory = theoretical_thomas_pair_correlation(
        propagation_pair_correlation.radius,
        cluster_parameters.parent_intensity,
        cluster_parameters.dispersal_sigma,
    )
    for i, radius_value in enumerate(propagation_pair_correlation.radius):
        pair_rows.append(
            {
                "radius": float(radius_value),
                "propagation_g_hat": float(propagation_pair_correlation.estimate[i]),
                "correlated_origin_g_hat": float(correlated_pair_correlation.estimate[i]),
                "theoretical_g": float(theory[i]),
                "propagation_pair_count": int(propagation_pair_correlation.pair_counts[i]),
                "correlated_origin_pair_count": int(correlated_pair_correlation.pair_counts[i]),
            }
        )
    _write_csv(
        output_dir / "appendix_b3_pair_correlation.csv",
        pair_rows,
        [
            "radius",
            "propagation_g_hat",
            "correlated_origin_g_hat",
            "theoretical_g",
            "propagation_pair_count",
            "correlated_origin_pair_count",
        ],
    )

    make_figure_9(
        output_dir,
        crossover_parameters,
        estimator_parameters,
        estimator_estimates,
        propagation_pair_correlation,
        correlated_pair_correlation,
        homogeneous_probability,
        cluster_parameters,
    )
    make_figure_10(
        output_dir,
        cluster_parameters,
        propagation,
        correlated,
        propagation_pair_correlation,
        correlated_pair_correlation,
    )

    parameter_record = {
        "master_seed": master_seed,
        "seed_strategy": "numpy.random.SeedSequence(master_seed).spawn(3)",
        "crossover": asdict(crossover_parameters),
        "known_center_estimator": asdict(estimator_parameters),
        "cluster_process": asdict(cluster_parameters),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy_version,
            "matplotlib": matplotlib.__version__,
        },
        "edge_treatment": "periodic square window with minimum-image distances",
        "pair_correlation_estimator": (
            "unordered annulus pair counts divided by homogeneous-Poisson expected counts"
        ),
    }
    (output_dir / "appendix_b_parameters.json").write_text(
        json.dumps(parameter_record, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    _write_results_markdown(
        output_dir / "appendix_b_results.md",
        homogeneous_probability,
        threshold,
        crossover_rows,
        estimator_rows,
        propagation,
        correlated,
        cluster_parameters,
        master_seed,
    )

    summary = {
        "output_dir": str(output_dir),
        "homogeneous_probability": homogeneous_probability,
        "probability_equivalent_threshold": threshold,
        "crossover_rows": crossover_rows,
        "estimator_rows": estimator_rows,
        "propagation_parent_count": len(propagation.parents),
        "propagation_observed_count": len(propagation.observed_points),
        "correlated_origin_count": len(correlated.observed_points),
    }
    return summary


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Reproduce Appendix B calculations, simulations, tables, and figures."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs"),
        help="Directory for generated figures, tables, and parameter records.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=MASTER_SEED,
        help=f"Master random seed (default: {MASTER_SEED}).",
    )
    parser.add_argument(
        "--replicates",
        type=int,
        default=EstimatorParameters().replicates,
        help="Monte Carlo replicates per sample size for Appendix B.2.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use 1,000 Appendix B.2 replicates for a fast smoke test.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    replicates = 1_000 if args.quick else args.replicates
    if replicates <= 0:
        parser.error("--replicates must be positive")

    summary = run_appendix_b(
        output_dir=args.output_dir,
        master_seed=args.seed,
        estimator_replicates=replicates,
    )

    print("Appendix B reproduction completed")
    print(f"Output directory: {summary['output_dir']}")
    print(f"Homogeneous predictive p_I: {summary['homogeneous_probability']:.12f}")
    print(
        "Probability-equivalent threshold T_I: "
        f"{summary['probability_equivalent_threshold']:.12f}"
    )
    for row in summary["crossover_rows"]:
        print(
            f"xi={row['xi']:g}: crossover radius="
            f"{row['crossover_radius']:.6f}"
        )
    print(
        "Propagation realization: "
        f"{summary['propagation_parent_count']} latent parents, "
        f"{summary['propagation_observed_count']} observed points"
    )
    print(
        "Correlated-origin realization: "
        f"{summary['correlated_origin_count']} independent-origin points"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
