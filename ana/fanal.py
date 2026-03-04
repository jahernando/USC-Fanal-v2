"""
ana.fanal
=========
Core analysis module for the FANAL double-beta neutrino-less decay search
(NEXT experiment, ¹³⁶Xe).

Workflow
--------
1. Define selection cuts and energy ranges.
2. Build an extended likelihood (ELL) object from MC templates.
3. Fit real or pseudo-data to extract the number of signal (bb0ν) events.
4. Compute test statistics (t_μ, q_μ, q₀) for hypothesis testing.
5. Run many toy Monte Carlo experiments and summarise results in a DataFrame.

Two fit strategies are supported:
  * ``ell``    – signal region only.
  * ``simell`` – simultaneous signal + control region.

Physical constants / defaults
------------------------------
erange  : (2.4, 2.7) MeV   – energy analysis window.
eroi    : (2.43, 2.48) MeV – signal Region of Interest.
eblob2  : 0.4 MeV          – minimum energy of the second blob.
NA      : Avogadro number.
abundance : 0.9             – ¹³⁶Xe isotopic abundance.
Qbb     : 2458 keV          – Q-value of ¹³⁶Xe 0νββ.
W       : 135.9 g/mol       – atomic weight of ¹³⁶Xe.
"""

import numpy  as np
import pandas as pd

from collections import namedtuple

import scipy.constants as constants
import scipy.stats     as stats

import core.utils as ut
import core.efit  as efit


# ---------------------------------------------------------------------------
# Sample labels and physical constants
# ---------------------------------------------------------------------------

ssamples  = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']

erange    = (2.4, 2.7)   # MeV – full analysis energy window
eroi      = (2.43, 2.48) # MeV – signal Region of Interest
eblob2    = 0.4          # MeV – min blob-2 energy

NA        = constants.Avogadro
abundance = 0.9          # ¹³⁶Xe isotopic fraction
Qbb       = 2458         # keV
W         = 135.9        # g/mol


# ---------------------------------------------------------------------------
# Physics helpers
# ---------------------------------------------------------------------------

def half_life(nbb, exposure, eff, abundance=abundance, W=W):
    """Compute the ¹³⁶Xe 0νββ half-life from the fitted signal count.

    Parameters
    ----------
    nbb       : float – fitted number of signal events in the RoI.
    exposure  : float – detector exposure in kg·y.
    eff       : float – total signal detection efficiency (fraction).
    abundance : float – isotopic abundance of ¹³⁶Xe.  Default 0.9.
    W         : float – atomic weight of ¹³⁶Xe in g/mol.  Default 135.9.

    Returns
    -------
    float – half-life in years.
    """
    return 1e3 * eff * abundance * (exposure / nbb) * (NA / W) * np.log(2.)


# ---------------------------------------------------------------------------
# Event selections
# ---------------------------------------------------------------------------

def selection_blind(df, eroi=eroi, eblob2=eblob2):
    """Return a boolean mask for the *blind* sample (events outside the RoI).

    Events are blinded if they fall in the energy RoI **and** pass the blob-2
    cut.  The returned mask selects events that are *not* blinded, restricted
    to the main energy analysis window.

    Parameters
    ----------
    df     : pandas.DataFrame
    eroi   : tuple(float, float) – energy RoI in MeV.
    eblob2 : float               – minimum blob-2 energy in MeV.

    Returns
    -------
    numpy.ndarray of bool
    """
    in_roi        = (df.track0_E >= eroi[0]) & (df.track0_E < eroi[1])
    blob2_pass    = df.blob2_E > eblob2
    blinded       = in_roi | blob2_pass
    in_ewindow    = (df.E >= erange[0]) & (df.E < erange[1])
    return (~blinded) & in_ewindow


# ---------------------------------------------------------------------------
# Monte Carlo experiment generation
# ---------------------------------------------------------------------------

def generate_mc_experiment(mcs, nevts):
    """Sample a pseudo-data set from MC templates with Poisson fluctuations.

    Parameters
    ----------
    mcs   : list[pandas.DataFrame] – MC samples (one per component).
    nevts : array-like             – expected number of events per component.

    Returns
    -------
    pandas.DataFrame – concatenated pseudo-data set.
    """
    ns   = [int(stats.poisson.rvs(ni, size=1)[0]) for ni in nevts]
    xmcs = [mc.sample(n=ni) for mc, ni in zip(mcs, ns)]
    return pd.concat(xmcs, ignore_index=True)


# ---------------------------------------------------------------------------
# Extended likelihood construction
# ---------------------------------------------------------------------------

def get_ell(mcs, refnames, refranges,
            varname='E', varrange=erange, bins=100):
    """Build an :class:`~core.efit.ExtComPDF` from MC template histograms.

    Parameters
    ----------
    mcs       : list[pandas.DataFrame] – MC samples.
    refnames  : str | list[str]        – variables used to select the PDF region.
    refranges : tuple | list[tuple]    – ranges for the PDF-region selection.
    varname   : str                    – variable to fit.  Default ``'E'``.
    varrange  : tuple(float, float)    – histogram range.  Default :data:`erange`.
    bins      : int                    – number of histogram bins.

    Returns
    -------
    ExtComPDF
    """
    # Select MC events in the reference region to build the PDF templates
    refmcs = [ut.selection_sample(mc, refnames, refranges) for mc in mcs]

    histos = [np.histogram(mc[varname], bins, range=varrange) for mc in refmcs]
    pdfs   = [stats.rv_histogram(h) for h in histos]
    return efit.ExtComPDF(pdfs)


# ---------------------------------------------------------------------------
# Signal-region fit
# ---------------------------------------------------------------------------

def prepare_fit_ell(mcs, nevts, varnames, varranges,
                    refnames=None, refranges=None,
                    varname='E', varrange=erange, bins=100):
    """Prepare a closure that fits data using a signal-region ELL.

    The returned function accepts a DataFrame and returns the fit result.

    Parameters
    ----------
    mcs       : list[pandas.DataFrame]
    nevts     : array-like – total expected events per MC component.
    varnames  : str | list[str]  – selection variables.
    varranges : tuple | list     – selection ranges.
    refnames  : str | list | None  – variables for PDF construction.
        Defaults to *varnames*.
    refranges : tuple | list | None – ranges for PDF construction.
        Defaults to *varranges*.
    varname   : str   – fit variable.  Default ``'E'``.
    varrange  : tuple – fit histogram range.  Default :data:`erange`.
    bins      : int   – histogram bins.  Default 100.

    Returns
    -------
    callable – ``fit(data) → (result, values, ell, nevts_exp)``
    """
    refnames  = varnames  if refnames  is None else refnames
    refranges = varranges if refranges is None else refranges

    # Expected counts in the signal selection for each component
    effs      = np.array([ut.selection_efficiency(mc, varnames, varranges)[0]
                          for mc in mcs])
    nevts_exp = effs * np.array(nevts)

    ell = get_ell(mcs, refnames, refranges, varname, varrange, bins)

    def _fit(data):
        data_sel  = ut.selection_sample(data, varnames, varranges)
        values    = data_sel[varname].values
        result    = ell.best_estimate(values, *nevts_exp)
        return result, values, ell, nevts_exp

    return _fit


# ---------------------------------------------------------------------------
# Simultaneous signal + control region fit
# ---------------------------------------------------------------------------

def prepare_fit_simell(mcs, nevts, varnames, varranges,
                       refnames, refranges, connames, conranges,
                       varname='E', varrange=erange, bins=100):
    """Prepare a simultaneous signal + control region ELL fit.

    Parameters
    ----------
    mcs       : list[pandas.DataFrame]
    nevts     : array-like – total expected events per component.
    varnames  : list – signal selection variables.
    varranges : list – signal selection ranges.
    refnames  : list – variables for signal PDF construction.
    refranges : list – ranges for signal PDF construction.
    connames  : list – variables for control PDF construction.
    conranges : list – ranges for control PDF construction.
    varname   : str   – fit variable.  Default ``'E'``.
    varrange  : tuple – fit range.  Default :data:`erange`.
    bins      : int   – histogram bins.

    Returns
    -------
    callable – ``fit(data) → (result, values, ell, nevts_exp)``
    """
    effs_sig    = np.array([ut.selection_efficiency(mc, varnames, varranges)[0]
                             for mc in mcs])
    effs_con    = np.array([ut.selection_efficiency(mc, connames, conranges)[0]
                             for mc in mcs])
    nevts_exp   = effs_sig * np.array(nevts)
    # Ratio of control to signal expected counts (for the simultaneous fit)
    ratios      = effs_con / effs_sig

    ell_signal  = get_ell(mcs, refnames, refranges, varname, varrange, bins)
    ell_control = get_ell(mcs, connames,  conranges, varname, varrange, bins)
    ell         = efit.SimulExtComPDF(ell_signal, ell_control, ratios)

    def _fit(data):
        sig_data = ut.selection_sample(data, varnames, varranges)
        con_data = ut.selection_sample(data, connames, conranges)
        values   = (sig_data[varname].values, con_data[varname].values)
        result   = ell.best_estimate(values, *nevts_exp)
        return result, values, ell, nevts_exp

    return _fit


# ---------------------------------------------------------------------------
# Log-likelihood scan and test statistics
# ---------------------------------------------------------------------------

def tmu_scan(values, pars, ell, sizes=2., nbins=50):
    """Profile log-likelihood scan over all parameters.

    Parameters
    ----------
    values : array-like – data (or tuple of arrays for simultaneous fit).
    pars   : array-like – best-fit parameters.
    ell    : ComPDF     – object with a ``loglike`` method.
    sizes  : float | array-like
        Width of the scan range in units of sqrt(n) for each parameter.
        A single float applies to all parameters.
    nbins  : int – number of scan points.

    Returns
    -------
    nis  : list[numpy.ndarray] – scan values per parameter.
    tmus : list[numpy.ndarray] – -2ΔlogL values per parameter.
    """
    npars = len(pars)
    sizes = npars * (sizes,) if isinstance(sizes, float) else sizes

    def _scan_range(n, size):
        """Linear range around the best-fit value n, width = size * sqrt(n)."""
        n0 = max(0., n - size * np.sqrt(n))
        n1 = n + size * np.sqrt(n) if n > 1 else 20.
        return np.linspace(n0, n1, nbins)

    result = []
    for i, (ni, si) in enumerate(zip(pars, sizes)):
        scan_pts = _scan_range(ni, si)
        tmu_vals = efit.llike_scan(values, ell, pars, scan_pts, i)
        result.append((scan_pts, tmu_vals))

    return ut.list_transpose(result)


def tmu_values(values, par_est, ell, par_exp):
    """Compute hypothesis-test statistics for a fitted experiment.

    Parameters
    ----------
    values  : array-like – data used in the fit.
    par_est : array-like – best-fit (estimated) parameters.
    ell     : ComPDF     – PDF object.
    par_exp : array-like – expected (nominal) parameter values.

    Returns
    -------
    tmun : float – t_μ comparing full vectors (ipos=-1).
    tmu  : float – t_μ for the signal parameter (ipos=0).
    qmu  : float – one-sided statistic q_μ (= t_μ if signal < expected, else 0).
    q0   : float – q₀ = t_μ(μ=0), tests the background-only hypothesis.
    """
    tmun = efit.tmu(values, ell, par_est, par_exp,     ipos=-1)
    tmu  = efit.tmu(values, ell, par_est, par_exp[0],  ipos=0)
    q0   = efit.tmu(values, ell, par_est, 0.,           ipos=0)
    # q_μ is zero when the estimated signal exceeds the hypothesised value
    qmu  = tmu if par_est[0] < par_exp[0] else 0.

    return tmun, tmu, qmu, q0


# ---------------------------------------------------------------------------
# Toy Monte Carlo experiments
# ---------------------------------------------------------------------------

# Named tuple to collect the results of one toy experiment
ExpResult = namedtuple(
    'ExpResult',
    ('nbb', 'nBi', 'nTl',       # fitted event counts
     'nbb0', 'nBi0', 'nTl0',    # expected event counts
     'tmun', 'tmu', 'qmu', 'q0') # test statistics
)


def prepare_experiment_ell(mcs, nevts, *args, **kargs):
    """Closure that generates and analyses one signal-region toy experiment.

    Parameters
    ----------
    mcs     : list[pandas.DataFrame] – MC templates.
    nevts   : array-like             – expected counts per component.
    *args   : positional args for :func:`prepare_fit_ell`.
    **kargs : keyword args for :func:`prepare_fit_ell`.

    Returns
    -------
    callable – ``experiment() → (success, mc_data, ExpResult | None)``
    """
    fit = prepare_fit_ell(mcs, nevts, *args, **kargs)
    return _prepare_experiment(mcs, nevts, fit)


def prepare_experiment_simell(mcs, nevts, *args, **kargs):
    """Closure that generates and analyses one simultaneous toy experiment.

    Parameters
    ----------
    mcs     : list[pandas.DataFrame]
    nevts   : array-like
    *args   : positional args for :func:`prepare_fit_simell`.
    **kargs : keyword args for :func:`prepare_fit_simell`.

    Returns
    -------
    callable – ``experiment() → (success, mc_data, ExpResult | None)``
    """
    fit = prepare_fit_simell(mcs, nevts, *args, **kargs)
    return _prepare_experiment(mcs, nevts, fit)


def _prepare_experiment(mcs, nevts, fit):
    """Internal factory: wraps a *fit* function in an experiment loop."""

    def _experiment():
        mc_data = generate_mc_experiment(mcs, nevts)
        result, values, ell, nevts_exp = fit(mc_data)

        if not result.success:
            return False, mc_data, None

        nevts_est = result.x
        tmu_vals  = tmu_values(values, nevts_est, ell, nevts_exp)
        eresult   = ExpResult(*nevts_est, *nevts_exp, *tmu_vals)

        return True, mc_data, eresult

    return _experiment


def run(experiment, size=1):
    """Run *size* toy experiments and return results in a DataFrame.

    Failed experiments (where the fit did not converge) are silently dropped.

    Parameters
    ----------
    experiment : callable – zero-argument function returning
                 ``(success, data, ExpResult)``.
    size       : int – number of toys.  Default 1.

    Returns
    -------
    pandas.DataFrame
        One row per successful experiment, columns = :attr:`ExpResult._fields`.
    """
    eresults = []
    for _ in range(size):
        success, _, eresult = experiment()
        if success:
            eresults.append(eresult)

    eresults = ut.list_transpose(eresults)
    return ut.list_to_df(eresults, ExpResult._fields)
