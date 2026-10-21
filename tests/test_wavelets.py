"""Component recovery exercises the wavelet route with the fallback disabled."""

from pathlib import Path

import numpy as np
import pytest

from laplacex import DecomposeConfig, decompose


@pytest.mark.parametrize("rates,amplitudes,grid", [
    ([2.], [1.5], "uniform"),
    ([.3, 3.], [2., .5], "uniform"),
    ([.01, 1., 100.], [1., 1., 1.], "log"),
    ([.01, 100.], [1., 1.], "log"),
    ([.03, 2., 60.], [.7, 1.3, .5], "log"),
    ([.5, .7], [2., .5], "uniform"),
    ([.02, 1., 50.], [1., 1., .001], "log"),
    ([.02, 1., 50.], [.001, 1., 1.], "log"),
    ([.02, 1., 50.], [1., .001, 1.], "log"),
    ([.5, 5.], [2., -.5], "uniform"),
])
def test_exact_component_identification(rates, amplitudes, grid):
    t = np.geomspace(.001, 1000, 600) if grid == "log" else np.linspace(.01, 20, 600)
    f = np.exp(-t[:, None]*rates)@amplitudes
    components = decompose(t, f, auto_deflate=False)
    assert len(components) == len(rates)
    np.testing.assert_allclose(sorted(components.rates), sorted(rates), rtol=2e-7)
    np.testing.assert_allclose(components.reconstruct(t), f, rtol=1e-12, atol=2e-14)
    assert all(c.frequency == 0 for c in components)


@pytest.mark.parametrize("seed", [7, 19, 41, 83, 127, 313])
def test_noisy_component_count(seed):
    t = np.geomspace(.002, 400, 600)
    rates = np.array([.03, 2., 60.])
    clean = np.exp(-t[:, None]*rates)@np.array([.7, 1.3, .5])
    f = clean+np.random.default_rng(seed).normal(0, .001, len(t))
    components = decompose(t, f, auto_deflate=False)
    assert len(components) == 3
    np.testing.assert_allclose(sorted(components.rates), rates, rtol=.015)
    assert np.sqrt(np.mean((components.reconstruct(t)-clean)**2)) < .0006


@pytest.mark.parametrize("rates,amplitudes", [
    ([.42, .58], [1.2, .8]),
    ([.9, 1.3], [.4, 2.]),
    ([.13, .18], [2., .2]),
    ([.7, 1.1], [1., -.3]),
])
@pytest.mark.parametrize("sigma", [0., 1e-6])
def test_unresolved_pairs_with_fallback_disabled(rates, amplitudes, sigma):
    t = np.linspace(.01, 20/min(rates), 600)
    clean = np.exp(-t[:, None]*rates)@amplitudes
    observed = clean+np.random.default_rng(20261031).normal(0, sigma, len(t))
    components = decompose(t, observed, auto_deflate=False)
    assert len(components) == 2
    np.testing.assert_allclose(sorted(components.rates), rates, rtol=2e-4 if sigma else 2e-10)
    assert np.sqrt(np.mean((components.reconstruct(t)-clean)**2)) < (sigma/2 if sigma else 1e-13)


@pytest.mark.parametrize("name", ["log41", "uniform47"])
@pytest.mark.parametrize("known_noise", [False, True])
def test_noisy_boundary_and_missing_component_regressions(name, known_noise):
    # Exact observations from seed 20261107, grid/index log/41 and uniform/47
    # in benchmarks/validate_wavelets.py. Preserve the actual noise realization.
    with np.load(Path(__file__).parent/"data/noisy_detection.npz") as fixture:
        t, observed = fixture[name+"_t"], fixture[name+"_y"]
        clean, rates = fixture[name+"_clean"], fixture[name+"_rates"]
    cfg = DecomposeConfig()
    if known_noise:
        cfg.noise_sigma = .001
    components = decompose(t, observed, cfg=cfg, auto_deflate=False)
    assert len(components) == len(rates)
    np.testing.assert_allclose(sorted(components.rates), rates, rtol=.04)
    assert np.sqrt(np.mean((components.reconstruct(t)-clean)**2)) < .0002


@pytest.mark.parametrize("sigma", [0., 1e-4])
def test_identifiable_fast_boundary_mode_is_retained(sigma):
    t = np.geomspace(.001, 1000, 600)
    rates = np.array([.2, 20., 2000.])
    clean = np.exp(-t[:, None]*rates)@np.array([.5, .4, .3])
    observed = clean+np.random.default_rng(20261201).normal(0, sigma, len(t))
    components = decompose(t, observed, auto_deflate=False)
    assert len(components) == 3
    np.testing.assert_allclose(sorted(components.rates), rates, rtol=.002 if sigma else 1e-10)


def test_mode_budget_applies_to_wavelets():
    t = np.geomspace(.001, 1000, 600)
    f = np.exp(-t[:, None]*[.01, 1., 100.]).sum(axis=1)
    assert len(decompose(t, f, max_components=2, auto_deflate=False)) <= 2
    with pytest.raises(ValueError):
        decompose(t, f, max_components=33)


@pytest.mark.parametrize("scale", [1e-200, 1e-100, 1., 1e100, 1e200])
def test_extreme_amplitude_scaling(scale):
    t = np.geomspace(.002, 400, 600)
    clean = np.exp(-t[:, None]*[.03, 2., 60.])@np.array([.7, 1.3, .5])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        components = decompose(t, clean*scale)
        assert len(components) == 3
        np.testing.assert_allclose(components.reconstruct(t)/scale, clean, atol=1e-13)


@pytest.mark.parametrize("seed", range(10))
def test_noise_and_constant_do_not_generate_decay_modes(seed):
    t = np.geomspace(.001, 1000, 600)
    noise = np.random.default_rng(seed).normal(0, .001, len(t))
    assert len(decompose(t, noise)) == 0
    components = decompose(t, 1+noise)
    assert len(components) == 1
    assert components[0].rate == 0
    assert abs(components[0].amplitude-1) < .0002


def test_broad_wavelet_configuration_stays_finite():
    cfg = DecomposeConfig()
    cfg.scale_min, cfg.scale_max = 1e100, 1e200
    t = np.linspace(.01, 10, 100)
    components = decompose(t, np.exp(-t), cfg=cfg)
    assert np.all(np.isfinite(components.reconstruct(t)))


@pytest.mark.parametrize("criterion", ["bic", "aicc"])
def test_information_criterion_configuration(criterion):
    cfg = DecomposeConfig()
    cfg.model_selection = criterion
    t = np.geomspace(.001, 1000, 600)
    f = np.exp(-t[:, None]*[.01, 1., 100.]).sum(axis=1)
    assert len(decompose(t, f, cfg=cfg, auto_deflate=False)) == 3
    cfg.model_selection = "invalid"
    with pytest.raises(ValueError, match="model_selection"):
        decompose(t, f, cfg=cfg)
