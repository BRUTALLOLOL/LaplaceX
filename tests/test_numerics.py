"""Accuracy, input safety, and API regressions with independent analytic references."""

import math
from importlib.metadata import version

import numpy as np
import pytest
from laplacex import (
    DecomposeConfig,
    ExponentialComponent,
    ExponentialComponents,
    LaplaceTransform,
    TransformConfig,
    decompose,
    inverse,
    laplace,
    utils,
)


def test_release_version_matches_installed_metadata_and_native_module():
    import laplacex
    from laplacex import _laplacex

    assert laplacex.__version__ == _laplacex.__version__ == version("LaplaceX") == "1.0.0"


@pytest.mark.parametrize("rate", [1e-12, 1., 1e12])
def test_forward_across_scales(rate):
    cfg = TransformConfig()
    cfg.t_max = 100 / rate
    cfg.tol = 1e-12 / rate
    value = laplace(lambda t: math.exp(-rate*t), rate, cfg=cfg)
    assert value == pytest.approx(.5/rate, rel=2e-10, abs=0)


def test_forward_respects_t_max():
    cfg = TransformConfig()
    cfg.t_max = .1
    s = .001
    expected = -math.expm1(-s*cfg.t_max)/s
    assert laplace(lambda t: 1, s, cfg=cfg) == pytest.approx(expected, rel=1e-12)


def test_forward_integrable_endpoint_singularity():
    assert laplace(lambda t: 1/math.sqrt(t), 2) == pytest.approx(
        math.sqrt(math.pi/2), rel=1e-10
    )


def test_forward_callback_cost():
    calls = 0

    def f(t):
        nonlocal calls
        calls += 1
        return math.exp(-2*t)

    assert laplace(f, 1) == pytest.approx(1/3, rel=1e-12)
    assert calls < 2000  # Broken embedded Gauss weights caused excessive subdivision.


@pytest.mark.parametrize("F, exact", [
    (lambda s: 1/(s+2), lambda t: np.exp(-2*t)),
    (lambda s: 1/s, lambda t: np.ones_like(t)),
    (lambda s: 1/s**2, lambda t: t),
    (lambda s: (s+.5)/((s+.5)**2+9), lambda t: np.exp(-.5*t)*np.cos(3*t)),
])
def test_inverse_analytic_pairs(F, exact):
    t = np.linspace(.1, 5, 100)
    np.testing.assert_allclose(inverse(F, t), exact(t), atol=1e-9, rtol=1e-9)


def test_inverse_oscillations_and_order_convergence():
    t = np.linspace(.1, 5, 100)
    F = lambda s: 10/(s*s+100)
    low_error = np.max(abs(inverse(F, t, N=128) - np.sin(10*t)))
    high_error = np.max(abs(inverse(F, t, N=256) - np.sin(10*t)))
    assert high_error < 1e-7
    assert high_error < low_error/100


def test_inverse_auto_scale_does_not_prefer_missing_oscillations():
    t = np.linspace(.1, 5, 100)
    actual = inverse(lambda s: 10/(s*s+100), t)
    assert np.max(abs(actual - np.sin(10*t))) < .1


def test_inverse_honors_sigma_and_b():
    # With sigma=1 and b=3, exp(-2t) is exactly the zeroth Laguerre function.
    points = np.array([0., .1, 1., 5.])
    np.testing.assert_allclose(inverse(lambda s: 1/(s+2), points, sigma=1, b=3, N=2),
                               np.exp(-2*points), atol=2e-13, rtol=2e-13)


def test_sampled_forward_origin_and_finite_window():
    t = np.linspace(.001, 20, 2000)
    np.testing.assert_allclose(laplace(np.exp(-t), [1, 2+1j], t=t),
                               1/(np.array([1, 2+1j])+1), atol=3e-8, rtol=0)
    t = np.linspace(0, 2, 20)
    assert laplace(np.ones_like(t), 1, t=t) == pytest.approx(-math.expm1(-2), rel=1e-12)


def test_strided_numpy_arrays():
    t = np.linspace(.01, 10, 200)
    f = np.exp(-t)
    s = np.array([1, 99, 2+1j, 99, 3+2j, 99], dtype=complex)[::2]
    np.testing.assert_allclose(laplace(lambda t: math.exp(-t), s), 1/(s+1), atol=1e-12)
    np.testing.assert_allclose(laplace(f[::2], s, t=t[::2]),
                               laplace(f[::2].copy(), s.copy(), t=t[::2].copy()), atol=1e-14)
    np.testing.assert_allclose(inverse(lambda s: 1/(s+1), t[::3]), np.exp(-t[::3]),
                               atol=1e-10)
    np.testing.assert_array_equal(utils.gen_step(t[::-1]), (t[::-1] >= 1).astype(float))


@pytest.mark.parametrize("engine", [None, LaplaceTransform()])
def test_shapes(engine):
    forward = laplace if engine is None else engine.laplace
    backward = inverse if engine is None else engine.inverse
    s = np.array([[1, 2+1j], [3, 4]])
    np.testing.assert_allclose(forward(lambda t: math.exp(-t), s), 1/(s+1), atol=1e-12)
    t = np.array([[.1, .2], [.5, 1]])
    np.testing.assert_allclose(backward(lambda s: 1/(s+1), t), np.exp(-t), atol=1e-10)
    assert backward(lambda s: 1/(s+1), 1.).shape == ()
    assert forward(lambda t: 1, np.empty((0, 2))).shape == (0, 2)


def test_cached_components_do_not_override_explicit_image():
    engine = LaplaceTransform()
    t = np.linspace(.01, 10, 200)
    engine.decompose(t, np.exp(-t))
    np.testing.assert_allclose(engine.inverse(lambda s: 1/(s+2), t), np.exp(-2*t),
                               atol=1e-10)
    with pytest.raises(ValueError, match="method"):
        engine.inverse(lambda s: 1/s, t, method="typo")


def test_oscillatory_component_image_and_phase():
    c = ExponentialComponent()
    c.amplitude, c.rate, c.frequency, c.phase = 2, .5, 3, .7
    components = ExponentialComponents([c])
    s = np.array([[1, 2+1j], [3-2j, 4]])
    expected = 2*((s+.5)*math.cos(.7)-3*math.sin(.7))/((s+.5)**2+9)
    np.testing.assert_allclose(components.reconstruct_laplace(s), expected, atol=2e-15)
    # Independent numerical integration of the time-domain component.
    numerical = laplace(lambda t: 2*math.exp(-.5*t)*math.cos(3*t+.7), s)
    np.testing.assert_allclose(components.reconstruct_laplace(s), numerical, atol=1e-10)


def test_decomposition_recovers_clean_rates():
    t = np.linspace(.01, 15, 600)
    f = 2*np.exp(-.3*t)+.5*np.exp(-3*t)
    components = decompose(t, f)
    np.testing.assert_allclose(components.reconstruct(t), f, atol=1e-7, rtol=1e-7)
    np.testing.assert_allclose(sorted(components.rates), [.3, 3], rtol=1e-5)
    assert sum(components.significances) == pytest.approx(1.)


@pytest.mark.parametrize("scale", [1e-12, 1., 1e12])
def test_decomposition_amplitude_scale(scale):
    t = np.linspace(.01, 20, 400)
    f = scale*(2*np.exp(-.5*t)+.5*np.exp(-5*t))
    components = decompose(t, f)
    np.testing.assert_allclose(components.reconstruct(t)/scale, f/scale, atol=1e-7, rtol=1e-7)


def test_decomposition_signed_amplitudes():
    t = np.linspace(.01, 20, 400)
    f = 2*np.exp(-.5*t)-.5*np.exp(-5*t)
    components = decompose(t, f)
    np.testing.assert_allclose(components.reconstruct(t), f, atol=1e-7, rtol=1e-7)


def test_decomposition_denoises_seeded_signal():
    t = np.linspace(.01, 15, 400)
    clean = np.exp(-t)
    noisy = clean + np.random.default_rng(42).normal(0, .01, len(t))
    prediction = decompose(t, noisy).reconstruct(t)
    assert np.sqrt(np.mean((prediction-clean)**2)) < .002


@pytest.mark.parametrize("t,f", [
    ([1, 2, 3, 4], [1, 2]),
    ([1, 3, 2, 4], [1, 2, 3, 4]),
    ([1, 2, 2, 4], [1, 2, 3, 4]),
    ([1, 2, float('nan'), 4], [1, 2, 3, 4]),
    ([1, 2, 3, 4], [1, float('inf'), 3, 4]),
    ([[1, 2], [3, 4]], [[1, 2], [3, 4]]),
])
def test_invalid_samples_raise(t, f):
    with pytest.raises(ValueError):
        laplace(f, 1, t=t)
    with pytest.raises(ValueError):
        decompose(t, f)


@pytest.mark.parametrize("field,value", [("tol", 0), ("rel_tol", -1),
                                         ("max_levels", -1), ("t_max", float('nan'))])
def test_invalid_transform_config(field, value):
    cfg = TransformConfig()
    setattr(cfg, field, value)
    with pytest.raises(ValueError):
        laplace(lambda t: 1, 1, cfg=cfg)


@pytest.mark.parametrize("field,value", [("n_scales", 1), ("n_log_points", 0),
                                         ("scale_min", 0), ("wavelet", "unknown")])
def test_invalid_decompose_config(field, value):
    cfg = DecomposeConfig()
    setattr(cfg, field, value)
    with pytest.raises(ValueError):
        decompose([1, 2, 3, 4], [4, 3, 2, 1], cfg=cfg)


@pytest.mark.parametrize("times", [[-1], [float('nan')], [float('inf')]])
def test_inverse_invalid_times(times):
    with pytest.raises(ValueError):
        inverse(lambda s: 1/s, times)


@pytest.mark.parametrize("terms", [-1, 0, 1, 4097])
def test_inverse_invalid_order(terms):
    with pytest.raises(ValueError):
        inverse(lambda s: 1/s, [1], N=terms)


def test_nonfinite_callback_raises():
    with pytest.raises(ValueError, match="Non-finite"):
        laplace(lambda t: float('nan'), 1)
    with pytest.raises(ValueError, match="Non-finite"):
        inverse(lambda s: complex('nan'), [1])


def test_information_criteria_exact_fit():
    assert utils.aic(0, 50, 2) == -math.inf
    assert utils.bic(0, 50, 2) == -math.inf


def test_zero_and_negative_log_ranges():
    t = np.array([1., 2., 3., 4.])
    assert utils.suggest_log_range(t, np.zeros(4)) == pytest.approx((.1, 40))
    assert utils.suggest_log_range(t, [-1, -.01, -.001, -.0001]) == pytest.approx((.1, 10))
