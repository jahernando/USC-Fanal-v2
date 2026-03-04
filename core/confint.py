"""
core.confint
============
Feldman-Cousins (FC) confidence-interval utilities for a Poisson process
with known background.

References
----------
Feldman & Cousins, Phys. Rev. D 57, 3873 (1998).

Functions
---------
fc_confsegment    : FC acceptance interval for a single μ (signal) value.
fc_confband       : FC acceptance band over a grid of μ values.
get_fc_confinterval : Build a closure that maps n_obs → FC confidence interval.
fca_segment       : FC segment from Monte Carlo ordering (t_μ-based).
"""

import numpy       as np
import scipy.stats as stats

import core.utils as ut


# Maximum number of standard deviations above (μ + b) used to auto-define
# the range of integer observations considered in the FC ordering.
_NSIGMA_AUTO = 5


def fc_confsegment(nu, bkg, cl=0.68, nrange=None):
    """Feldman-Cousins acceptance interval for one signal value *nu*.

    Constructs the set of observations {n} that belong to the FC confidence
    belt at confidence level *cl* for a Poisson signal *nu* on top of
    background *bkg*.

    Parameters
    ----------
    nu     : float – expected number of signal events.
    bkg    : float – expected number of background events.
    cl     : float – confidence level.  Default 0.68.
    nrange : tuple(int, int) | None
        Range of integer observations ``(n_min, n_max)`` to consider.
        If None, the range is set automatically using :data:`_NSIGMA_AUTO`.

    Returns
    -------
    tuple(int, int)
        Minimum and maximum observation counts included in the acceptance
        interval at the requested confidence level.
    """
    if nrange is None:
        nmax   = bkg + nu + _NSIGMA_AUTO * np.sqrt(bkg + nu)
        nrange = (0, int(nmax) + 1)

    ns      = np.arange(*nrange)
    # Best-fit signal given n observed (physical constraint: μ ≥ 0)
    nuhats  = np.maximum(ns - bkg, 0.)

    ps      = stats.poisson.pmf(ns, bkg + nu)
    ps_best = stats.poisson.pmf(ns, bkg + nuhats)

    # FC ordering variable: likelihood ratio
    ts   = -2. * (np.log(ps) - np.log(ps_best))
    vals = sorted(zip(ts, ps, ns))

    _, sorted_ps, sorted_ns = ut.list_transpose(vals)
    cum_ps = np.cumsum(sorted_ps)
    assert cum_ps[-1] > cl, \
        f'Observation range {nrange} is too small to reach CL={cl}'

    # Include observations until cumulative probability exceeds cl
    i = 0
    while cum_ps[i] < cl:
        i += 1

    included = sorted_ns[:i + 1]
    return int(np.min(included)), int(np.max(included))


def fc_confband(nus, bkg, cl=0.68, nrange=None):
    """FC acceptance band over an array of signal values *nus*.

    Parameters
    ----------
    nus    : array-like – grid of signal values to evaluate.
    bkg    : float – expected background.
    cl     : float – confidence level.  Default 0.68.
    nrange : tuple(int, int) | None – observation range (see :func:`fc_confsegment`).

    Returns
    -------
    n0s : numpy.ndarray(int) – lower edge of the acceptance interval.
    n1s : numpy.ndarray(int) – upper edge of the acceptance interval.
    """
    segs      = [fc_confsegment(nu, bkg, cl, nrange) for nu in nus]
    n0s, n1s  = ut.list_transpose(segs)
    return np.array(n0s, dtype=int), np.array(n1s, dtype=int)


def get_fc_confinterval(nus, bkg, cl=0.68, nrange=None):
    """Build a closure that computes the FC confidence interval for *n_obs*.

    Parameters
    ----------
    nus    : array-like – fine grid of signal (μ) values.
    bkg    : float – expected background.
    cl     : float – confidence level.  Default 0.68.
    nrange : tuple(int, int) | None

    Returns
    -------
    callable
        ``ci(n_obs)`` → ``numpy.ndarray([mu_low, mu_high])`` — the FC
        confidence interval on the signal strength for *n_obs* observed events.
        Accepts a scalar or a numpy array of observations.
    """
    nus  = np.asarray(nus)
    n0s, n1s = fc_confband(nus, bkg, cl, nrange)

    def ci(n_obs):
        """FC confidence interval for observed count *n_obs*.

        Parameters
        ----------
        n_obs : int | numpy.ndarray

        Returns
        -------
        numpy.ndarray – shape (2,) or (2, N) for array input.
        """
        if isinstance(n_obs, np.ndarray):
            results = [ci(ni) for ni in n_obs]
            return np.array(ut.list_transpose(results))

        mu_upper = np.max(nus[n0s <= n_obs])
        mu_lower = np.min(nus[n1s >= n_obs])
        return np.array((mu_lower, mu_upper))

    return ci


def fca_segment(tmus, ns, cl=0.9):
    """FC acceptance interval from a Monte Carlo sample using t_μ ordering.

    Parameters
    ----------
    tmus : array-like – FC ordering variable (e.g. t_μ values) for each trial.
    ns   : array-like – observable values (e.g. n_bb) for each trial.
    cl   : float      – confidence level.  Default 0.9.

    Returns
    -------
    numpy.ndarray([n_min, n_max])
        Lower and upper bounds of the acceptance interval.
    """
    sorted_vals = sorted(zip(tmus, ns))
    _, sorted_ns = ut.list_transpose(sorted_vals)

    # Number of trials to include to reach the desired coverage
    n_include = cl * len(tmus)
    ipos      = int(n_include)
    # Round to nearest integer
    if n_include - ipos >= 0.5:
        ipos += 1

    included = sorted_ns[:ipos]
    return np.array((np.min(included), np.max(included)))
