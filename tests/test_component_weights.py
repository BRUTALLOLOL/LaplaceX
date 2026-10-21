import numpy as np
import pytest
from laplacex import ExponentialComponents, decompose


def test_weights_and_legacy_significances_agree_for_signed_components():
    t=np.geomspace(.01,20,300)
    fit=decompose(t,2*np.exp(-.4*t)-.3*np.exp(-3*t),auto_deflate=False)
    expected=np.abs(fit.amplitudes)/sum(np.abs(fit.amplitudes))
    np.testing.assert_allclose(fit.weights,expected)
    np.testing.assert_array_equal(fit.weights,fit.significances)
    assert sum(fit.weights)==pytest.approx(1.)


def test_empty_components_have_empty_weights():
    assert ExponentialComponents([]).weights.shape==(0,)
