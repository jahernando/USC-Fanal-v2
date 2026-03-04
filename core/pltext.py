"""
core.pltext
===========
Extensions to :mod:`matplotlib.pyplot` for physics analysis plots.

Provides thin wrappers around ``plt.hist``, histogram fitting with overlay,
residual plots, and DataFrame inspection utilities.

Functions
---------
style          : Set a default colour-cycle style.
plt_text       : Add a text box to the current axes.
canvas         : Create a figure with N sub-plots and return a navigator.
karg           : Set a keyword-argument default without overwriting existing ones.
hist           : Histogram with optional statistics annotation and legend.
hfit           : Fit + plot a histogram to a function.
hresiduals     : Compute and plot normalised residuals.
hfitres        : Fit, plot, and display residuals in one call.
df_inspect     : Histogram all (or selected) columns of a DataFrame.
dfs_inspect    : Compare columns of multiple DataFrames.
df_corrmatrix  : Plot the absolute correlation matrix of a DataFrame.
"""

import numpy as np

import core.utils as ut
import core.hfit  as hfitm

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cycler import cycler


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

def style():
    """Apply a simple colour-cycle style to Matplotlib."""
    plt.rcParams['axes.prop_cycle'] = cycler(color='kbgrcmy')
    plt.style.context('seaborn-colorblind')


def plt_text(comment, x=0.05, y=0.7, **kargs):
    """Add a text box in axis-relative coordinates to the current axes.

    Parameters
    ----------
    comment : str – text to display.
    x, y    : float – position in axes fraction coordinates.
    **kargs : forwarded to :meth:`matplotlib.axes.Axes.text`.
    """
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    plt.gca().text(x, y, comment,
                   transform=plt.gca().transAxes, bbox=props, **kargs)


# ---------------------------------------------------------------------------
# Canvas / sub-plot navigation
# ---------------------------------------------------------------------------

def canvas(ns: int, ny: int = 2,
           height: float = 5., width: float = 6.):
    """Create a figure with *ns* sub-plots arranged in *ny* columns.

    Parameters
    ----------
    ns     : int   – total number of sub-plots.
    ny     : int   – number of columns.  Default 2.
    height : float – height of each sub-plot in inches.  Default 5.
    width  : float – width of each sub-plot in inches.  Default 6.

    Returns
    -------
    subplot : callable(iplot, dim='2d') → (nx, ny)
        Call ``subplot(i)`` (1-indexed) to activate the i-th sub-plot.
        Pass ``dim='3d'`` for three-dimensional axes.
    """
    nx = int(ns / ny + ns % ny)
    plt.figure(figsize=(width * ny, height * nx))

    def subplot(iplot, dim='2d'):
        assert iplot <= nx * ny, \
            f'iplot={iplot} exceeds the canvas size ({nx}×{ny})'
        plt.subplot(nx, ny, iplot)
        if dim == '3d':
            code = nx * 100 + ny * 10 + iplot
            plt.gcf().add_subplot(code, projection=dim)
        return nx, ny

    return subplot


# ---------------------------------------------------------------------------
# Keyword helpers
# ---------------------------------------------------------------------------

def karg(name, value, kargs):
    """Set *kargs[name] = value* only if *name* is not already in *kargs*.

    Parameters
    ----------
    name  : str
    value : any
    kargs : dict – modified in-place.

    Returns
    -------
    dict – the (possibly updated) *kargs*.
    """
    if name not in kargs:
        kargs[name] = value
    return kargs


# ---------------------------------------------------------------------------
# Histogram
# ---------------------------------------------------------------------------

def hist(x: np.ndarray, bins: int,
         stats: bool = True,
         xylabels=None,
         grid: bool = True,
         ylog: bool = False,
         **kargs):
    """Plot a histogram with optional statistics annotation.

    Parameters
    ----------
    x        : array-like – data values.
    bins     : int | array-like – number of bins or bin edges.
    stats    : bool – if True, append N / mean / std to the legend label.
    xylabels : str | tuple(str, str) | None
        If a string, sets the x-axis label.
        If a tuple ``(xlabel, ylabel)``, sets both axes labels.
    grid     : bool – show grid.  Default True.
    ylog     : bool – use log y-axis.  Default False.
    **kargs  : forwarded to :func:`matplotlib.pyplot.hist`.
        Special keys consumed here:
        * ``stats_format`` (str) – format string for mean/std.  Default ``'6.3f'``.

    Returns
    -------
    tuple – the return value of :func:`matplotlib.pyplot.hist`.
    """
    kargs.setdefault('histtype', 'step')

    if stats:
        val_range = kargs.get('range', None)
        fmt       = kargs.pop('stats_format', '6.3f')
        ss        = ut.str_stats(x, val_range=val_range, fmt=fmt)
        if 'label' in kargs:
            kargs['label'] += '\n' + ss
        else:
            kargs['label'] = ss

    c = plt.hist(x, bins, **kargs)

    if xylabels is not None:
        if isinstance(xylabels, str):
            plt.xlabel(xylabels)
        elif isinstance(xylabels, tuple):
            plt.xlabel(xylabels[0])
            plt.ylabel(xylabels[1])

    if 'label' in kargs:
        plt.legend()

    if grid:
        plt.grid(True)
    if ylog:
        plt.yscale('log')

    return c


# ---------------------------------------------------------------------------
# Histogram fit + overlay
# ---------------------------------------------------------------------------

def hfit(x, bins, fun, guess=None, range=None,
         parnames=None, fmt='6.2f', **kargs):
    """Fit a histogram of *x* to *fun* and overlay the fit curve.

    Parameters
    ----------
    x        : array-like
    bins     : int | array-like
    fun      : callable | str – function to fit (see :mod:`core.hfit`).
    guess    : tuple | None   – initial parameter values.
    range    : tuple | None   – histogram range.
    parnames : list[str] | None – parameter labels for the legend.
    fmt      : str | None
        Format string for parameter values in the legend.
        Pass None to suppress parameter annotation.
    **kargs  : forwarded to :func:`hist` and :func:`matplotlib.pyplot.plot`.

    Returns
    -------
    pars  : numpy.ndarray – best-fit parameters.
    epars : numpy.ndarray – parameter uncertainties.
    """
    fun, guess, fnames = hfitm._predefined_function(fun, guess, x)

    # Plot histogram without stats (fit will provide the annotation)
    ys, xs, _ = hist(x, bins, range=range, stats=False, **kargs)

    pars, epars = hfitm.hfit(x, bins, fun, guess, range)
    xcs         = 0.5 * (xs[1:] + xs[:-1])

    if fmt is not None:
        parnames = parnames if parnames is not None else fnames
        ss       = hfitm.str_parameters(pars, epars, parnames, fmt=fmt)
        kargs['label'] = ss if 'label' not in kargs else kargs['label'] + '\n' + ss

    plt.plot(xcs, fun(xcs, *pars), **kargs)
    if 'label' in kargs:
        plt.legend()

    return pars, epars


def hresiduals(x, bins, fun, pars, **kargs):
    """Compute and plot normalised residuals of a histogram fit.

    Parameters
    ----------
    x    : array-like
    bins : int | array-like
    fun  : callable | str
    pars : tuple – fit parameters.
    **kargs : forwarded to :func:`matplotlib.pyplot.bar`.

    Returns
    -------
    res   : numpy.ndarray – normalised residuals.
    edges : numpy.ndarray – bin edges.
    chi2  : float         – total chi-squared.
    ndf   : int           – degrees of freedom.
    """
    res, edges, chi2, ndf = hfitm.hresiduals(x, bins, fun=fun, pars=pars, **kargs)

    xcs   = 0.5 * (edges[:-1] + edges[1:])
    widths = edges[1:] - edges[:-1]

    karg('label', r'$\chi^2$/ndf {:6.3f}'.format(chi2 / ndf), kargs)
    # Remove histogram-specific keys that bar() does not accept
    for key in ('range', 'density', 'weights', 'normed'):
        kargs.pop(key, None)

    plt.bar(xcs, res, width=widths, **kargs)
    plt.legend()
    plt.grid()

    return res, edges, chi2, ndf


def hfitres(x, bins, fun, guess=None, **kargs):
    """Fit, overlay and display residuals in a split axes layout.

    Parameters
    ----------
    x     : array-like
    bins  : int | array-like
    fun   : callable | str
    guess : tuple | None
    **kargs : forwarded to :func:`hfit` and :func:`hresiduals`.

    Returns
    -------
    pars  : numpy.ndarray
    epars : numpy.ndarray
    chi2  : float
    ndf   : int
    """
    pars, epars = hfit(x, bins, fun, guess=guess, **kargs)

    ax      = plt.gca()
    divider = make_axes_locatable(ax)
    ax2     = divider.append_axes('bottom', size='20%', pad=0)
    ax.figure.add_axes(ax2)

    _, _, chi2, ndf = hresiduals(x, bins, fun, pars, **kargs)

    return pars, epars, chi2, ndf


# ---------------------------------------------------------------------------
# DataFrame inspection
# ---------------------------------------------------------------------------

def df_inspect(df, labels=None, bins=100, ranges={}, ncolumns=2):
    """Histogram the columns of a DataFrame.

    Parameters
    ----------
    df       : pandas.DataFrame
    labels   : list[str] | None – columns to plot.  Defaults to all columns.
    bins     : int – number of bins.  Default 100.
    ranges   : dict – mapping from column name to ``(x0, x1)`` histogram range.
    ncolumns : int – number of plot columns in the canvas.  Default 2.
    """
    if labels is None:
        labels = list(df.columns)

    subplot = canvas(len(labels), ncolumns)
    for i, label in enumerate(labels):
        subplot(i + 1)
        values = ut.remove_nan(df[label].values)
        xrange = ranges.get(label, None)
        hist(values, bins, range=xrange)
        plt.xlabel(label)

    plt.tight_layout()


def dfs_inspect(dfs, dfnames=None, labels=None,
                bins=100, ranges={}, ncolumns=2):
    """Overlay histograms of the same column from multiple DataFrames.

    Parameters
    ----------
    dfs      : list[pandas.DataFrame]
    dfnames  : list[str] | None – legend labels.  Defaults to ``'0', '1', ...``.
    labels   : list[str] | None – columns to plot.  Defaults to all columns.
    bins     : int
    ranges   : dict
    ncolumns : int
    """
    ndfs    = len(dfs)
    dfnames = [str(i) for i in range(ndfs)] if dfnames is None else dfnames
    labels  = list(dfs[0].columns) if labels is None else labels

    subplot = canvas(len(labels), ncolumns)
    for i, xlabel in enumerate(labels):
        subplot(i + 1)
        for name, df in zip(dfnames, dfs):
            values = ut.remove_nan(df[xlabel].values)
            xrange = ranges.get(xlabel, None)
            hist(values, bins, range=xrange, label=name, density=True)
        plt.xlabel(xlabel)

    plt.tight_layout()


def df_corrmatrix(xdf, xlabels):
    """Plot the absolute correlation matrix of selected columns.

    Parameters
    ----------
    xdf     : pandas.DataFrame
    xlabels : list[str] – columns to include.
    """
    _df  = xdf[xlabels]
    corr = _df.corr()

    fig = plt.figure(figsize=(12, 10))
    plt.matshow(np.abs(corr), fignum=fig.number, cmap='Greys')
    plt.xticks(range(_df.shape[1]), _df.columns, fontsize=14, rotation=45)
    plt.yticks(range(_df.shape[1]), _df.columns, fontsize=14)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
