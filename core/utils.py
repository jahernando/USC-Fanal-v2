"""
core.utils
==========
General-purpose utilities for array manipulation, statistical summaries,
range-based selections and efficiency calculations over NumPy arrays and
Pandas DataFrames.

Functions
---------
list_transpose        : Transpose a list of lists.
list_to_df            : Convert a list of lists into a pandas DataFrame.
remove_nan            : Remove NaN entries from an array.
in_range              : Boolean mask for values inside an interval.
centers               : Mid-points between consecutive partition edges.
stats                 : Basic statistics (mean, std, counts) with optional range.
str_stats             : Human-readable statistics string.
efficiency            : Binomial efficiency and its uncertainty.
selection             : Boolean mask for DataFrame rows satisfying variable ranges.
selection_efficiency  : Efficiency of a DataFrame selection.
selection_sample      : Sub-DataFrame of rows that pass a selection.
"""

import numpy  as np
import pandas as pd

from functools import reduce


# ---------------------------------------------------------------------------
# List utilities
# ---------------------------------------------------------------------------

def list_transpose(ll):
    """Transpose a list-of-lists (m lists of length n → n lists of length m).

    Parameters
    ----------
    ll : list[list]
        Input rectangular list of lists; all inner lists must have the same
        length.

    Returns
    -------
    list[list]
        Transposed list of lists.
    """
    m  = len(ll[0])
    lt = [[row[i] for row in ll] for i in range(m)]
    return lt


def list_to_df(ll, names):
    """Convert a list of columns into a :class:`pandas.DataFrame`.

    Parameters
    ----------
    ll    : list[array-like]
        One element per column; each element is an iterable of column values.
    names : list[str]
        Column names.  Must have the same length as *ll*.

    Returns
    -------
    pandas.DataFrame
    """
    assert len(ll) == len(names), \
        f'required same number of lists ({len(ll)}) and names ({len(names)})'
    return pd.DataFrame(dict(zip(names, ll)))


# ---------------------------------------------------------------------------
# Array utilities
# ---------------------------------------------------------------------------

def remove_nan(vals: np.ndarray) -> np.ndarray:
    """Return *vals* with NaN entries removed.

    Parameters
    ----------
    vals : numpy.ndarray

    Returns
    -------
    numpy.ndarray
    """
    return vals[~np.isnan(vals)]


def in_range(vals: np.ndarray,
             val_range=None,
             upper_limit_in: bool = False) -> np.ndarray:
    """Boolean mask selecting elements of *vals* that lie inside *val_range*.

    Parameters
    ----------
    vals         : numpy.ndarray
    val_range    : None | scalar | tuple(float, float)
        * ``None``          → all elements selected (True everywhere).
        * ``(x0, x1)``      → ``x0 <= vals < x1`` (or ``<= x1`` if
          *upper_limit_in* is True).
        * scalar *v*        → ``vals == v`` (exact equality).
    upper_limit_in : bool, optional
        If True the upper bound is inclusive (``vals <= x1``).  Default False.

    Returns
    -------
    numpy.ndarray of bool
    """
    if val_range is None:
        return np.ones(len(vals), dtype=bool)

    if isinstance(val_range, (list, tuple)):
        x0, x1 = val_range[0], val_range[1]
        sel_lo = vals >= x0
        sel_hi = vals <= x1 if upper_limit_in else vals < x1
        return sel_lo & sel_hi

    return vals == val_range


def centers(xs: np.ndarray) -> np.ndarray:
    """Return the mid-points between consecutive elements of *xs*.

    Parameters
    ----------
    xs : numpy.ndarray
        Partition edges (length n).

    Returns
    -------
    numpy.ndarray
        Array of length n-1 with the centres of each interval.
    """
    return 0.5 * (xs[1:] + xs[:-1])


# ---------------------------------------------------------------------------
# Statistical summaries
# ---------------------------------------------------------------------------

def stats(vals: np.ndarray, val_range=None):
    """Compute basic statistics of *vals*, optionally restricted to *val_range*.

    Parameters
    ----------
    vals      : array-like
    val_range : None | tuple(float, float), optional
        If given, only elements inside the range are included.

    Returns
    -------
    evts  : int   – number of selected elements.
    mean  : float – mean of selected elements.
    std   : float – standard deviation of selected elements.
    oevts : int   – number of elements outside the range.
    """
    vals  = np.array(vals)
    vals  = remove_nan(vals)
    sel   = in_range(vals, val_range)
    vv    = vals[sel]
    mean  = np.mean(vv)
    std   = np.std(vv)
    evts  = len(vv)
    oevts = len(vals) - evts
    return evts, mean, std, oevts


def str_stats(vals, val_range=None, fmt='6.2f'):
    """Return a multi-line string summarising statistics of *vals*.

    Parameters
    ----------
    vals      : array-like
    val_range : None | tuple(float, float), optional
    fmt       : str, optional
        Python format specification for mean and std.  Default ``'6.2f'``.

    Returns
    -------
    str
    """
    evts, mean, std, _ = stats(vals, val_range)
    s  = f'entries {evts}\n'
    s += (f'mean {{0:{fmt}}}').format(mean) + '\n'
    s += (f'std  {{0:{fmt}}}').format(std)
    return s


# ---------------------------------------------------------------------------
# Efficiency
# ---------------------------------------------------------------------------

def efficiency(sel, n=None):
    """Binomial efficiency and its uncertainty from a boolean selection mask.

    Parameters
    ----------
    sel : array-like of bool
        Boolean array where ``True`` marks selected elements.
    n   : int or None, optional
        Denominator.  If None, ``len(sel)`` is used.

    Returns
    -------
    eff  : float – efficiency = sum(sel) / n.
    ueff : float – binomial uncertainty = sqrt(eff * (1 - eff) / n).
    """
    n    = n if n is not None else len(sel)
    eff  = np.sum(sel) / n
    ueff = np.sqrt(eff * (1.0 - eff) / n)
    return eff, ueff


# ---------------------------------------------------------------------------
# DataFrame selections
# ---------------------------------------------------------------------------

def selection(df, varname, varrange, oper=np.logical_and):
    """Boolean mask selecting rows of *df* where column(s) fall inside range(s).

    Parameters
    ----------
    df       : pandas.DataFrame
    varname  : str | list[str]
        Column name, or list of column names for a multi-variable selection.
    varrange : tuple | list[tuple]
        Range ``(x0, x1)`` for a single variable, or a list of ranges for
        multiple variables.  Each range is passed to :func:`in_range`.
    oper     : callable, optional
        Logical reduction applied when multiple variables are given.
        Default :func:`numpy.logical_and`.

    Returns
    -------
    numpy.ndarray of bool
        Same length as *df*.
    """
    _isiter = lambda x: isinstance(x, (list, tuple))

    if _isiter(varname):
        assert len(varname) == len(varrange), \
            'required same length of variables and ranges'
        sels = [selection(df, ivar, irange)
                for ivar, irange in zip(varname, varrange)]
        return reduce(oper, sels)

    return in_range(df[varname].values, varrange)


def selection_efficiency(df, varname, varrange):
    """Efficiency of a selection applied to *df*.

    Parameters
    ----------
    df       : pandas.DataFrame
    varname  : str | list[str]
    varrange : tuple | list[tuple]

    Returns
    -------
    eff  : float – efficiency.
    ueff : float – binomial uncertainty.
    """
    return efficiency(selection(df, varname, varrange))


def selection_sample(df, varname, varrange):
    """Sub-DataFrame of rows in *df* that pass the selection.

    Parameters
    ----------
    df       : pandas.DataFrame
    varname  : str | list[str]
    varrange : tuple | list[tuple]

    Returns
    -------
    pandas.DataFrame
    """
    return df[selection(df, varname, varrange)]
