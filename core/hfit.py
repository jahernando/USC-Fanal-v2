"""
core.hfit
=========
Histogram fitting utilities.

Pre-defined fit functions (and their initial-guess and parameter-name helpers):

  * ``gaus``     – Gaussian: ``a * exp(-(x-b)^2 / (2*c^2))``
  * ``line``     – Linear:   ``a*x + b``
  * ``exp``      – Exponential: ``a * exp(-b*x)``
  * ``gausline`` – Gaussian + linear background
  * ``gausexp``  – Gaussian + exponential background

Each function ``f<name>`` is paired with a guess function ``g<name>`` that
returns sensible initial parameters from raw data, and a name list ``n<name>``
for plot annotations.

Main API
--------
hfit       : Fit a histogram to a function; returns parameters and errors.
hresiduals : Compute normalised residuals and chi²/ndf.
hfitres    : Fit + residuals in one call.
"""

import numpy          as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Registry of pre-defined function names
# ---------------------------------------------------------------------------

functions = ['gaus', 'line', 'exp', 'gausline', 'gausexp']

# Used by _predefined_function to resolve 'f<name>' and 'g<name>' look-ups
# from within this module.
_current_module = __import__(__name__)


def _predefined_function(fun, guess, x):
    """Resolve a string function name to the callable, initial guess and labels.

    Parameters
    ----------
    fun   : callable | str
        If a string, it must be one of :data:`functions`.
    guess : tuple | None
        Initial parameter guess.  If None *and* fun is a string, the
        corresponding ``g<name>`` helper is called with *x*.
    x     : array-like
        Data values (used only to auto-generate guesses).

    Returns
    -------
    fun    : callable
    guess  : tuple
    fnames : list[str] | None
        LaTeX parameter names, or None when *fun* was a callable.
    """
    fnames = None
    if isinstance(fun, str):
        assert fun in functions, \
            f"Unknown function '{fun}'. Choose from {functions}."
        if guess is None:
            # call the corresponding guess function, e.g. ggaus(x)
            guess = getattr(_current_module.hfit, 'g' + fun)(x)
        fnames = getattr(_current_module.hfit, 'n' + fun)
        fun    = getattr(_current_module.hfit, 'f' + fun)

    return fun, guess, fnames


# ---------------------------------------------------------------------------
# Core fitting routines
# ---------------------------------------------------------------------------

def hfit(x, bins, fun, guess=None, range=None):
    """Fit a histogram of *x* to *fun* using non-linear least squares.

    Parameters
    ----------
    x     : array-like
        Data values used to build the histogram.
    bins  : int | array-like
        Number of bins or bin edges.
    fun   : callable | str
        Function to fit.  If a string, one of :data:`functions`.
    guess : tuple, optional
        Initial parameter values.  Auto-generated for pre-defined functions.
    range : tuple(float, float), optional
        Histogram range.  Defaults to ``(min(x), max(x))``.

    Returns
    -------
    pars  : numpy.ndarray – best-fit parameters.
    epars : numpy.ndarray – parameter uncertainties (sqrt of covariance diagonal).
    """
    fun, guess, _ = _predefined_function(fun, guess, x)

    xrange     = range if range is not None else (np.min(x), np.max(x))
    yc, edges  = np.histogram(x, bins, range=xrange)
    xc         = 0.5 * (edges[1:] + edges[:-1])

    pars, pcov = optimize.curve_fit(fun, xc, yc, guess)
    return pars, np.sqrt(np.diag(pcov))


def hresiduals(x, bins, fun='gaus', pars=None, **kargs):
    """Compute normalised residuals ``(f(x) - counts) / sqrt(counts)``.

    Parameters
    ----------
    x    : array-like
    bins : int | array-like
    fun  : callable | str, optional
    pars : tuple, optional
        Parameters of *fun*.  Required when *fun* is a callable.
    **kargs :
        Extra keyword arguments forwarded to :func:`numpy.histogram`.

    Returns
    -------
    res  : numpy.ndarray – normalised residuals per bin.
    edges: numpy.ndarray – histogram bin edges.
    chi2 : float – total chi-squared.
    ndf  : int   – degrees of freedom = (non-empty bins) - (number of params).
    """
    fun, pars, _ = _predefined_function(fun, pars, x)

    ys, edges = np.histogram(x, bins, **kargs)
    xcs       = 0.5 * (edges[:-1] + edges[1:])
    # Use at least 1 to avoid division by zero in empty bins
    yerr      = np.maximum(np.sqrt(ys), 1.)
    res       = (fun(xcs, *pars) - ys) / yerr
    chi2      = float(np.sum(res * res))
    ndf       = int(np.sum(ys > 0)) - len(pars)

    return res, edges, chi2, ndf


def hfitres(x, bins, fun, guess=None, **kargs):
    """Fit + residuals combined.

    Returns
    -------
    pars  : numpy.ndarray – best-fit parameters.
    epars : numpy.ndarray – parameter uncertainties.
    chi2  : float – chi-squared.
    ndf   : int   – degrees of freedom.
    """
    pars, epars         = hfit      (x, bins, fun, guess=guess, **kargs)
    _, _, chi2, ndf     = hresiduals(x, bins, fun,  pars,       **kargs)
    return pars, epars, chi2, ndf


def str_parameters(pars, cov_pars, parnames=None, fmt='6.2f'):
    """Format fit parameters and uncertainties as a multi-line string.

    Parameters
    ----------
    pars     : array-like – parameter values.
    cov_pars : array-like – parameter uncertainties (1-sigma).
    parnames : list[str] | None
        Parameter labels.  Defaults to ``$a_0$, $a_1$, ...``.
    fmt      : str, optional
        Python format specification.  Default ``'6.2f'``.

    Returns
    -------
    str
    """
    s = ''
    for i, (par, err) in enumerate(zip(pars, cov_pars)):
        name = r'$a_' + str(i) + '$' if parnames is None else parnames[i]
        s += name + ' = '
        s += (f'{{0:{fmt}}}').format(par)  + r'$\pm$'
        s += (f'{{0:{fmt}}}').format(err)  + '\n'
    return s


# ---------------------------------------------------------------------------
# Pre-defined functions, guesses and parameter names
# ---------------------------------------------------------------------------

def fgaus(x, a, b, c):
    """Gaussian: ``a * exp(-(x-b)^2 / (2*c^2))``.

    Returns np.inf if *c* <= 0 (unphysical width), which forces the optimiser
    away from invalid regions.
    """
    if c <= 0.:
        return np.inf
    return a * np.exp(-(x - b)**2 / (2. * c**2))


def ggaus(x):
    """Initial guess for a Gaussian fit: (amplitude, mean, std)."""
    return (len(x), np.mean(x), np.std(x))


ngaus = [r'$N_\mu$', r'$\mu$', r'$\sigma$']


# ---- Linear ------------------------------------------------------------------

def fline(x, a, b):
    """Linear function: ``a*x + b``."""
    return a * x + b


def gline(x):
    """Initial guess for a linear fit from a 2-bin histogram slope."""
    ys, edges = np.histogram(x, 2)
    xc        = 0.5 * (edges[1:] + edges[:-1])
    a         = (ys[1] - ys[0]) / (xc[1] - xc[0])
    b         = ys[0] - a * xc[0]
    return a, b


nline = ['a', 'b']


# ---- Exponential -------------------------------------------------------------

def fexp(x, a, b):
    """Exponential: ``a * exp(-b*x)``."""
    return a * np.exp(-b * x)


def gexp(x):
    """Initial guess for an exponential fit from a 2-bin histogram ratio."""
    ys, edges = np.histogram(x, 2)
    xcs       = 0.5 * (edges[1:] + edges[:-1])
    dx        = xcs[1] - xcs[0]
    # Protect against log(0)
    ys        = np.maximum(ys, 1)
    b         = -(np.log(ys[1]) - np.log(ys[0])) / dx
    a         = ys[0] * np.exp(b * xcs[0])
    return a, b


nexp = [r'$N_\tau$', r'$\tau$']


# ---- Composite: Gaussian + linear background ---------------------------------

def fgausline(x, na, mu, sig, a, b):
    """Gaussian plus linear background."""
    return fgaus(x, na, mu, sig) + fline(x, a, b)


def ggausline(x):
    """Initial guess for Gaussian + linear background."""
    return list(ggaus(x)) + list(gline(x))


ngausline = ngaus + nline


# ---- Composite: Gaussian + exponential background ----------------------------

def fgausexp(x, na, mu, sig, nb, tau):
    """Gaussian plus exponential background."""
    return fgaus(x, na, mu, sig) + fexp(x, nb, tau)


def ggausexp(x):
    """Initial guess for Gaussian + exponential background."""
    return list(ggaus(x)) + list(gexp(x))


ngausexp = ngaus + nexp
