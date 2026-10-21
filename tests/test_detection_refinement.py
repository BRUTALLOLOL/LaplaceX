"""Regression observations selected from the exposed detection-development sweep."""
from pathlib import Path
import numpy as np
import pytest
from laplacex import decompose


@pytest.mark.parametrize('name', ['missed_pair','uniform_extra','log_extra','log_separated'])
def test_bic_pair_is_not_missed_or_replaced_by_cancelling_triple(name):
    # Original research seed 2026102101; grid, ratio, sigma, repeat:
    # uniform/1.1/.001/0; uniform/2/.05/2; log/2/.05/7; log/3/.05/14.
    # These are regression fixtures, not independent validation observations.
    with np.load(Path(__file__).parent/'data/component_refinement.npz') as data:
        t,y,truth,sigma=(data[name+'_'+key] for key in ('t','y','rates','sigma'))
    fit=decompose(t,y,auto_deflate=False)
    assert len(fit)==2
    assert np.all(np.isfinite(fit.amplitudes))
    if name!='missed_pair':
        np.testing.assert_allclose(sorted(fit.rates),truth,rtol=.1)
    # Compare with an independently minimized one-component Gaussian BIC.
    # A bounded golden-section solve is sufficient for the scalar rate.
    low,high=np.log(.005/t[-1]),np.log(8/t[0]);fraction=(np.sqrt(5)-1)/2
    def rss(x):
        e=np.exp(-np.exp(x)*t);r=y-e*(e@y)/(e@e)
        return r@r
    a,b=high-fraction*(high-low),low+fraction*(high-low)
    for _ in range(100):
        if rss(a)<rss(b):high,b=b,a;a=high-fraction*(high-low)
        else:low,a=a,b;b=low+fraction*(high-low)
    single=min(rss(a),rss(b));residual=y-fit.reconstruct(t);pair=residual@residual
    assert len(t)*np.log(pair/single)+2*np.log(len(t))<0

@pytest.mark.parametrize('case', ['4','6'])
def test_signed_four_mode_signal_survives_inadmissible_projected_minimum(case):
    # Signed supplement seed 2026102209, uniform grid, K=4, sigma=0,
    # repeats 4 and 6. A three-mode projected fit coalesces and is inadmissible;
    # the joint-coordinate retry must allow growth to the fourth real mode.
    with np.load(Path(__file__).parent/'data/signed_refinement.npz') as data:
        t,y,rates=(data[case+'_'+key] for key in ('t','y','rates'))
    fit=decompose(t,y,auto_deflate=False)
    assert len(fit)==4
    np.testing.assert_allclose(sorted(fit.rates),rates,rtol=2e-7)
    np.testing.assert_allclose(fit.reconstruct(t),y,atol=2e-13,rtol=1e-12)
