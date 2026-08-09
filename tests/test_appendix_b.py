from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from correlated_abiogenesis.appendix_b import (
    ClusterParameters,
    EstimatorParameters,
    crossover_radius,
    estimate_pair_correlation_periodic,
    gamma_poisson_predictive_probability,
    probability_equivalent_poisson_mean,
    rare_basin_integrated_mean,
    run_appendix_b,
    simulate_known_center_estimator,
    simulate_periodic_thomas_process,
    summarize_estimator_simulation,
)


def test_gamma_poisson_predictive_matches_appendix_value() -> None:
    p_i = gamma_poisson_predictive_probability(2.0, 1, 3.0, 7.0, 1.0)
    expected = 1.0 - (10.0 / 11.0) ** 3
    assert math.isclose(p_i, expected, rel_tol=0.0, abs_tol=1e-14)


def test_crossover_satisfies_equal_mean_threshold() -> None:
    p_i = gamma_poisson_predictive_probability(2.0, 1, 3.0, 7.0, 1.0)
    threshold = probability_equivalent_poisson_mean(p_i)
    for xi in (1.0, 2.0, 5.0):
        radius = crossover_radius(1_000_000.0, 1.0, xi, p_i)
        delta = rare_basin_integrated_mean(radius, 1_000_000.0, 1.0, xi)
        assert math.isclose(delta, threshold, rel_tol=1e-12, abs_tol=1e-12)


def test_known_center_estimator_is_nearly_unbiased() -> None:
    parameters = EstimatorParameters(true_xi=2.0, sample_sizes=(30,), replicates=20_000)
    estimates = simulate_known_center_estimator(np.random.default_rng(123), parameters)
    summary = summarize_estimator_simulation(estimates, parameters)[0]
    assert abs(summary["mean_xi_hat"] - 2.0) < 0.03
    assert summary["interval_width"] > 0.0


def test_periodic_cluster_simulation_and_pair_correlation() -> None:
    parameters = ClusterParameters(
        window_size=20.0,
        parent_intensity=0.15,
        mean_offspring=3.0,
        dispersal_sigma=0.7,
        maximum_pair_distance=5.0,
        pair_distance_bins=10,
    )
    realization = simulate_periodic_thomas_process(np.random.default_rng(321), parameters)
    assert len(realization.parents) > 0
    assert len(realization.observed_points) > 0
    assert np.all(realization.observed_points >= 0.0)
    assert np.all(realization.observed_points < parameters.window_size)

    edges = np.linspace(0.0, parameters.maximum_pair_distance, parameters.pair_distance_bins + 1)
    estimate = estimate_pair_correlation_periodic(
        realization.observed_points,
        parameters.window_size,
        edges,
    )
    assert len(estimate.radius) == parameters.pair_distance_bins
    assert np.all(np.isfinite(estimate.estimate))
    assert np.all(estimate.poisson_expected_counts > 0.0)


def test_quick_end_to_end_run(tmp_path: Path) -> None:
    output_dir = tmp_path / "appendix_b"
    summary = run_appendix_b(output_dir, master_seed=7, estimator_replicates=500)
    assert summary["propagation_parent_count"] > 0
    assert summary["correlated_origin_count"] > 0

    expected_files = {
        "appendix_b1_crossover.csv",
        "appendix_b2_estimator_summary.csv",
        "appendix_b3_pair_correlation.csv",
        "appendix_b_parameters.json",
        "appendix_b_results.md",
        "figure_9_appendix_b.png",
        "figure_9_appendix_b.pdf",
        "figure_10_matched_constructions.png",
        "figure_10_matched_constructions.pdf",
    }
    assert expected_files.issubset({path.name for path in output_dir.iterdir()})
