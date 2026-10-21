"""
laplacex.api
============
Convenience free functions that delegate to a default
:class:`~laplacex.transform.LaplaceTransform` instance.
These are the functions exported at the top-level ``laplacex`` namespace.
"""

from __future__ import annotations

from typing import Callable, Optional, Union
import numpy as np

from laplacex._laplacex import (
    TransformConfig,
    DecomposeConfig,
    forward_batch    as _forward_batch,
    forward_discrete as _forward_discrete,
    inverse_weeks    as _inverse_weeks,
    decompose        as _decompose_cpp,
    decompose_deflate as _decompose_deflate_cpp,
)
from laplacex.components import ExponentialComponents


def laplace(
    f: Union[Callable[[float], float], np.ndarray],
    s: Union[complex, np.ndarray],
    t: Optional[np.ndarray] = None,
    *,
    cfg: Optional[TransformConfig] = None,
) -> Union[complex, np.ndarray]:
    """
    Compute the forward Laplace transform F(s) = ∫₀^∞ f(t) exp(-st) dt.

    Parameters
    ----------
    f : callable or array-like
        The function f(t).  For an array, ``t`` must be supplied.
    s : complex or array-like of complex
        One or more complex frequencies at which F(s) is evaluated.
    t : array-like, optional
        Time points when ``f`` is an array.
    cfg : TransformConfig, optional
        Quadrature configuration.  Uses library defaults when *None*.

    Returns
    -------
    complex or ndarray of complex
        F(s).

    Examples
    --------
    >>> import numpy as np
    >>> from laplacex import laplace
    >>> # F(s) for f(t) = exp(-2t) should be 1/(s+2)
    >>> s = np.array([1+0j, 2+0j, 3+1j])
    >>> F = laplace(lambda t: np.exp(-2*t), s)
    >>> np.allclose(F, 1/(s+2), atol=1e-6)
    True
    """
    _cfg = cfg or TransformConfig()
    scalar_s = np.isscalar(s) or (hasattr(s, 'shape') and s.shape == ())
    s_input = np.asarray(s, dtype=complex)
    s_arr = s_input.reshape(-1)

    if callable(f):
        result = _forward_batch(f, s_arr, _cfg)
    else:
        if t is None:
            raise ValueError("When `f` is an array, `t` must be provided.")
        result = _forward_discrete(
            np.asarray(t, dtype=float),
            np.asarray(f, dtype=float),
            s_arr, _cfg,
        )

    return complex(result[0]) if scalar_s else np.asarray(result).reshape(s_input.shape)


def inverse(
    F: Optional[Callable[[complex], complex]],
    t_eval: np.ndarray,
    *,
    sigma: float = 0.0,
    b: float = 0.0,
    N: int = 64,
    components: Optional[ExponentialComponents] = None,
    cfg: Optional[TransformConfig] = None,
) -> np.ndarray:
    """
    Compute the inverse Laplace transform f(t) from F(s).

    Parameters
    ----------
    F : callable
        The Laplace image F(s) : complex → complex.
    t_eval : array-like
        Time points at which f(t) is desired.
    sigma : float
        Laguerre shift parameter (auto-selected when 0).
    b : float
        Laguerre scale parameter (auto-selected when 0).
    N : int
        Number of Laguerre terms (default 64).
    components : ExponentialComponents, optional
        If supplied, use the component representation F(s) = Σ A_k/(s+λ_k)
        for the reconstruction instead of the Weeks method.
    cfg : TransformConfig, optional

    Returns
    -------
    ndarray of float
        f(t) at each point in ``t_eval``.

    Examples
    --------
    >>> import numpy as np
    >>> from laplacex import inverse
    >>> # Inverse of F(s)=1/(s+2) is exp(-2t)
    >>> t = np.linspace(0.1, 5, 100)
    >>> f_rec = inverse(lambda s: 1/(s+2), t)
    >>> np.allclose(f_rec, np.exp(-2*t), atol=1e-3)
    True
    """
    t_arr = np.asarray(t_eval, dtype=float)
    if not np.all(np.isfinite(t_arr)) or np.any(t_arr < 0):
        raise ValueError("Inverse evaluation times must be finite and nonnegative")
    _cfg  = cfg or TransformConfig()

    if components is not None:
        return components.reconstruct(t_arr)

    if not callable(F):
        raise TypeError("F must be callable when components are not supplied")
    return np.asarray(_inverse_weeks(F, t_arr.reshape(-1), sigma, b, N, _cfg)).reshape(
        t_arr.shape
    )


def decompose(
    t: np.ndarray,
    f: np.ndarray,
    *,
    cfg: Optional[DecomposeConfig] = None,
    auto_deflate: bool = True,
    max_components: int = 5,
    energy_tol: float = 1e-6,
) -> ExponentialComponents:
    """
    Decompose a time-series into a sum of exponential components by
    performing wavelet analysis in the logarithmic time plane.

    Persistent Mexican-Hat zero-crossing ridges initialize a joint rate and
    amplitude fit. With ``auto_deflate=True`` (default), excess residual
    energy beyond both the requested tolerance and estimated noise triggers
    a second search using rate-grid proposals.
    Accept the candidate only when it improves the reconstruction.

    Parameters
    ----------
    t : array-like
        Strictly increasing time values (all t > 0).
    f : array-like
        Signal values at ``t``.
    cfg : DecomposeConfig, optional
        Algorithm configuration.
    auto_deflate : bool
        Try a rate-grid search when residual energy exceeds the tolerance and noise bound.
    max_components : int
        Maximum components in both the wavelet and rate-grid stages (1 to 32).
    energy_tol : float
        Relative residual-energy threshold for trying the rate-grid fallback.

    Returns
    -------
    ExponentialComponents
        Detected components, sorted by significance (descending).

    Examples
    --------
    >>> import numpy as np
    >>> from laplacex import decompose
    >>> t = np.linspace(0.01, 10, 500)
    >>> f = 2.0 * np.exp(-0.5 * t) + 0.5 * np.exp(-5 * t)
    >>> ec = decompose(t, f)
    >>> len(ec) >= 1
    True
    """
    _cfg  = cfg or DecomposeConfig()
    t_arr = np.asarray(t, dtype=float)
    f_arr = np.asarray(f, dtype=float)

    if max_components <= 0 or not np.isfinite(energy_tol) or not 0 < energy_tol < 1:
        raise ValueError("Need max_components > 0 and 0 < energy_tol < 1")
    comps = _decompose_cpp(t_arr, f_arr, _cfg, max_components)
    if auto_deflate:
        f_hat = ExponentialComponents(comps).reconstruct(t_arr)
        scale = np.max(np.abs(f_arr))
        if scale == 0:
            return ExponentialComponents(comps)
        residual_energy = np.sum((f_arr/scale - f_hat/scale) ** 2)
        original_energy = np.sum((f_arr/scale) ** 2)

        if _cfg.noise_sigma >= 0:
            noise_sigma = _cfg.noise_sigma/scale
        else:
            weight = (t_arr[1:-1]-t_arr[:-2])/(t_arr[2:]-t_arr[:-2])
            normalized = f_arr/scale
            differences = (
                normalized[1:-1]-(1-weight)*normalized[:-2]-weight*normalized[2:]
            )/np.sqrt(1+(1-weight)**2+weight**2)
            noise_sigma = np.median(np.abs(differences-np.median(differences)))/0.6744897501960817
        noise_energy = (len(t_arr)+3*np.sqrt(2*len(t_arr)))*noise_sigma**2
        if residual_energy > max(energy_tol * original_energy, noise_energy):
            candidate = _decompose_deflate_cpp(
                t_arr, f_arr, _cfg,
                max_components, energy_tol,
            )
            candidate_hat = ExponentialComponents(candidate).reconstruct(t_arr)
            if np.sum((f_arr/scale - candidate_hat/scale) ** 2) < residual_energy:
                comps = candidate

    return ExponentialComponents(comps)


def reconstruct(
    components: ExponentialComponents,
    t_eval: np.ndarray,
) -> np.ndarray:
    """
    Reconstruct f(t) = Σ A_k exp(-λ_k t) [cos(ω_k t + φ_k)] from
    an :class:`ExponentialComponents` object.

    Parameters
    ----------
    components : ExponentialComponents
    t_eval : array-like

    Returns
    -------
    ndarray of float

    Examples
    --------
    >>> import numpy as np
    >>> from laplacex import decompose, reconstruct
    >>> t = np.linspace(0.01, 10, 300)
    >>> f = np.exp(-t)
    >>> ec = decompose(t, f)
    >>> f_rec = reconstruct(ec, t)
    >>> np.allclose(f_rec, f, atol=0.05)
    True
    """
    return components.reconstruct(t_eval)
