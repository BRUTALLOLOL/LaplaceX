"""Transparent specialist baselines and Gaussian local identifiability bounds.

These are benchmark implementations, not additional public LaplaceX APIs.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy import linalg, optimize


@dataclass
class Fit:
    rates: np.ndarray
    amplitudes: np.ndarray
    rss: float

    def values(self, t):
        return np.exp(-np.outer(t, self.rates)) @ self.amplitudes


def projected(t, y, log_rates, jacobian=False):
    """Exact full-rank variable-projection residual and analytic Jacobian."""
    rates = np.exp(log_rates)
    e = np.exp(-np.outer(t, rates))
    q, r = linalg.qr(e, mode='economic', check_finite=False)
    # SVD least squares also handles a numerically rank-deficient iterate.
    a, _, rank, _ = linalg.lstsq(e, y, cond=1e-13, check_finite=False)
    residual = e @ a-y
    if not jacobian:
        return residual, a
    de = -np.outer(t, rates)*e
    first = de*a
    if rank < len(rates):
        # Return a well-defined local finite-difference derivative at collisions.
        step=1e-5
        j=np.column_stack([(projected(t,y,log_rates+np.eye(len(rates))[k]*step)[0]
                           -projected(t,y,log_rates-np.eye(len(rates))[k]*step)[0])/(2*step)
                          for k in range(len(rates))])
    else:
        correction = np.diag(de.T @ residual)
        j = first-q @ (q.T @ first)-q @ linalg.solve_triangular(
            r.T, correction, lower=True, check_finite=False)
    return residual, a, j


def refine_varpro(t, y, initial):
    low, high = .005/t[-1], 8/t[0]
    start=np.log(np.clip(initial,low*1.001,high/1.001))
    last_x=None; cached=None
    def calculate(x):
        nonlocal last_x,cached
        if last_x is None or not np.array_equal(x,last_x):
            cached=projected(t,y,x,True); last_x=x.copy()
        return cached
    result=optimize.least_squares(lambda x:calculate(x)[0],start,
        jac=lambda x:calculate(x)[2],bounds=(np.log(low),np.log(high)),
        ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=300)
    residual,amps=projected(t,y,result.x)
    order=np.argsort(result.x)
    return Fit(np.exp(result.x)[order],amps[order],float(residual@residual))


def bic(fit, y):
    return len(y)*np.log(max(fit.rss,1e-28*float(y@y))/len(y))+2*len(fit.rates)*np.log(len(y))


def variable_projection(t, y, known_count=None, budget=5):
    scale=max(np.max(np.abs(y)),1e-100); normalized=y/scale
    best=None
    for k in ([known_count] if known_count is not None else range(1,budget+1)):
        candidates=[]
        for margin in (1.,4.):
            initial=np.geomspace(margin/t[-1],1/(margin*t[0]),k)
            candidates.append(refine_varpro(t,normalized,initial))
        fit=min(candidates,key=lambda f:f.rss)
        if best is None or bic(fit,normalized)<bic(best,normalized): best=fit
        if known_count is None and fit.rss<=1e-28*float(normalized@normalized): break
    return Fit(best.rates,best.amplitudes*scale,best.rss*scale*scale)


def matrix_pencil(t, y, refine=True, budget=5):
    step=np.diff(t)
    if not np.allclose(step,step[0],rtol=1e-7,atol=1e-12):
        raise ValueError('Matrix pencil requires originally uniform times; interpolation is not used')
    scale=max(np.max(np.abs(y)),1e-100); normalized=y/scale
    width=min(96,len(t)//2)
    windows=np.lib.stride_tricks.sliding_window_view(normalized,width+1)
    h0=windows[:,:-1]; h1=windows[:,1:]
    u,s,vh=linalg.svd(h0,full_matrices=False,check_finite=False)
    best=None
    for k in range(1,min(budget,len(s))+1):
        if s[k-1]<=s[0]*1e-14: break
        reduced=(u[:,:k].T@h1@vh[:k].T)/s[:k][None,:]
        poles=linalg.eigvals(reduced,check_finite=False)
        if np.any(abs(poles.imag)>1e-8) or np.any(poles.real<=0) or np.any(poles.real>=1):
            continue
        rates=np.sort(-np.log(poles.real)/step[0])
        if refine:
            fit=refine_varpro(t,normalized,rates)
        else:
            e=np.exp(-np.outer(t,rates)); amps=linalg.lstsq(e,normalized,check_finite=False)[0]
            residual=e@amps-normalized; fit=Fit(rates,amps,float(residual@residual))
        if best is None or bic(fit,normalized)<bic(best,normalized): best=fit
        if best.rss<=1e-28*float(normalized@normalized): break
    if best is None: raise ValueError('No admissible real decaying matrix-pencil model')
    return Fit(best.rates,best.amplitudes*scale,best.rss*scale*scale)


def cramer_rao(t, amplitudes, rates, sigma):
    """Known-K, known-sigma Gaussian bound, with amplitudes as nuisance parameters.

    Return log-rate standard deviations; locally these are relative rate SDs.
    This bound is not a guarantee for biased estimators or model selection.
    """
    e=np.exp(-np.outer(t,rates))
    j=np.column_stack((e,-np.outer(t,rates)*e*amplitudes))
    norms=np.linalg.norm(j,axis=0)
    if np.any(norms==0): return {'rank_deficient':True,'log_rate_sd':None}
    _,s,vh=linalg.svd(j/norms,full_matrices=False,check_finite=False)
    condition=float(s[0]/s[-1])
    if s[-1]<s[0]*1e-12:
        return {'rank_deficient':True,'scaled_jacobian_condition':condition,'log_rate_sd':None}
    covariance=((vh.T*(sigma/s)**2)@vh)/np.outer(norms,norms)
    return {'rank_deficient':False,'scaled_jacobian_condition':condition,
            'log_rate_sd':np.sqrt(np.diag(covariance)[len(rates):]).tolist()}
