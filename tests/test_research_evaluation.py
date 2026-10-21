"""Independent checks on the scientific evaluation's numerical machinery."""
import numpy as np
import pytest

pytest.importorskip('scipy', reason='research baselines require the optional SciPy dependency')

from benchmarks.research.baselines import cramer_rao, matrix_pencil, projected, variable_projection


def test_projected_jacobian_against_central_differences():
    t=np.geomspace(.01,20,140)
    y=.8*np.exp(-.4*t)+.3*np.exp(-2*t)+.001*np.sin(t)
    x=np.log([.5,1.8])
    j=projected(t,y,x,True)[2]
    step=1e-5
    numerical=np.column_stack([(projected(t,y,x+np.eye(2)[k]*step)[0]
                               -projected(t,y,x-np.eye(2)[k]*step)[0])/(2*step)
                              for k in range(2)])
    np.testing.assert_allclose(j,numerical,rtol=2e-6,atol=2e-9)


@pytest.mark.parametrize('method',[variable_projection,matrix_pencil])
def test_specialist_baselines_recover_clean_pair(method):
    t=np.linspace(.01,20,180)
    y=1.2*np.exp(-.3*t)+.4*np.exp(-2*t)
    fit=method(t,y)
    np.testing.assert_allclose(fit.rates,[.3,2],rtol=1e-7)
    np.testing.assert_allclose(fit.amplitudes,[1.2,.4],rtol=1e-7)


def test_matrix_pencil_rejects_nonuniform_sampling():
    t=np.geomspace(.01,20,100)
    with pytest.raises(ValueError,match='uniform'):
        matrix_pencil(t,np.exp(-t))


def test_cramer_rao_single_mode_matches_two_parameter_information():
    t=np.linspace(.01,10,150); a=np.array([1.7]); rate=np.array([.8]); sigma=.01
    e=np.exp(-rate[0]*t); d=-a[0]*rate[0]*t*e
    information=np.array([[e@e,e@d],[d@e,d@d]])/sigma**2
    expected=np.sqrt(np.linalg.inv(information)[1,1])
    result=cramer_rao(t,a,rate,sigma)
    assert result['log_rate_sd'][0]==pytest.approx(expected,rel=1e-12)
    assert cramer_rao(t,a,rate,2*sigma)['log_rate_sd'][0]==pytest.approx(2*expected)


def test_cramer_rao_does_not_report_zero_uncertainty_for_duplicate_rates():
    result=cramer_rao(np.linspace(.01,10,100),np.array([1.,.5]),np.array([1.,1.]),.01)
    assert result['rank_deficient'] and result['log_rate_sd'] is None
