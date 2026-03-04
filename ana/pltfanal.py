"""
ana.pltfanal
============
Plotting utilities for the FANAL double-beta neutrino-less decay analysis.

Functions
---------
plot_fit_ell       : Overlay fitted PDF on energy histogram with residuals.
plot_fit_simell    : Call :func:`plot_fit_ell` for signal and control regions.
plot_tmu_scan      : Profile log-likelihood scan plots.
plot_nevts         : Distributions of fitted event counts over toy experiments.
plot_gaus_domain   : Check that test-statistic distributions follow chi-squared.
plot_exps_fc_confint : Feldman-Cousins confidence band vs. true signal count.
plot_exps_z0       : Discovery significance Z₀ vs. true signal count.
plot_contributions : Data histogram with per-component MC overlays.
"""

import numpy  as np
import scipy.stats as stats

import core.utils   as ut
import core.confint as confint
import core.pltext  as pltext

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

ssamples = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']
nbins    = 100
erange   = (2.4, 2.7)  # MeV


# ---------------------------------------------------------------------------
# ELL fit visualisation
# ---------------------------------------------------------------------------

def plot_fit_ell(x, par, ell,
                 bins=nbins,
                 parnames=ssamples,
                 plot_residuals=True,
                 title=''):
    """Plot data histogram with the best-fit PDF and per-component curves.

    Parameters
    ----------
    x              : array-like – energy values (data or pseudo-data).
    par            : array-like – best-fit event counts per component.
    ell            : ExtComPDF  – fitted extended likelihood object.
    bins           : int        – number of histogram bins.  Default 100.
    parnames       : list[str]  – component labels for the legend.
    plot_residuals : bool       – if True, attach a residual panel.  Default True.
    title          : str        – figure title.
    """
    subplot = pltext.canvas(1, 1, 8, 10)
    subplot(1)

    counts, edges = np.histogram(x, bins)
    centers  = 0.5 * (edges[1:] + edges[:-1])
    bin_width = centers[1] - centers[0]
    ecounts  = np.sqrt(counts)
    has_data = ecounts > 0

    plt.errorbar(centers[has_data], counts[has_data], yerr=ecounts[has_data],
                 marker='o', ls='', label='data')

    # Total fitted PDF
    norm_total = np.sum(par) * bin_width
    plt.plot(centers, norm_total * ell.pdf(centers, *par), label='ELL fit')

    # Individual component PDFs
    for ni, ipdf, name in zip(par, ell.pdfs, parnames):
        norm = ni * bin_width
        plt.plot(centers, norm * ipdf.pdf(centers),
                 label=f'{name} : {ni:6.2f}')

    plt.legend()
    plt.grid()
    plt.title(title)

    if plot_residuals:
        ax      = plt.gca()
        divider = make_axes_locatable(ax)
        ax2     = divider.append_axes('bottom', size='20%', pad=0)
        ax.figure.add_axes(ax2)
        # Wrap ell.pdf to match the signature expected by hresiduals
        fun = lambda x, *p: norm_total * ell.pdf(x, *p)
        pltext.hresiduals(x, bins, fun, par)


def plot_fit_simell(values, pars, simell, **kargs):
    """Plot signal and control region fits side by side.

    Parameters
    ----------
    values : tuple(array-like, array-like) – (signal data, control data).
    pars   : array-like – best-fit signal-region event counts.
    simell : SimulExtComPDF – simultaneous ELL object.
    **kargs : forwarded to :func:`plot_fit_ell`.
    """
    signal, control = values[0], values[1]
    plot_fit_ell(signal,  pars,
                 simell.ell_signal,
                 title='signal data',  **kargs)
    plot_fit_ell(control, pars * simell.ratios,
                 simell.ell_control,
                 title='control data', **kargs)


# ---------------------------------------------------------------------------
# Log-likelihood scan
# ---------------------------------------------------------------------------

def plot_tmu_scan(nis, tmus, cls=(0.68, 0.9), titles=ssamples):
    """Plot the profile log-likelihood scan for each parameter.

    Parameters
    ----------
    nis    : list[array-like] – scan points per parameter.
    tmus   : list[array-like] – -2ΔlogL values per parameter.
    cls    : tuple – confidence levels to mark with horizontal lines.
    titles : list[str] – subplot titles (one per parameter).
    """
    npars   = len(nis)
    subplot = pltext.canvas(npars, npars)

    for i, (ni, tmu) in enumerate(zip(nis, tmus)):
        subplot(i + 1)
        plt.plot(ni, tmu)
        for cl in cls:
            t0 = stats.chi2.ppf(cl, 1)
            plt.axhline(t0, ls='--', label=f'CL {100*cl:.0f}%')
        plt.xlabel('number of events')
        plt.ylabel(r'$\Delta -2 \mathrm{log}\,\mathcal{L}$')
        plt.grid()
        plt.title(titles[i])

    plt.tight_layout()


# ---------------------------------------------------------------------------
# Toy-experiment distributions
# ---------------------------------------------------------------------------

def plot_nevts(nevts, nbins=50, labels=ssamples):
    """Histogram the fitted event count distributions from toy experiments.

    Parameters
    ----------
    nevts  : list[array-like] – fitted counts from many experiments, one per component.
    nbins  : int – histogram bins.  Default 50.
    labels : list[str] – component labels.
    """
    npars   = len(nevts)
    subplot = pltext.canvas(npars, npars)

    for i, (nv, label) in enumerate(zip(nevts, labels)):
        subplot(i + 1)
        pltext.hist(nv, nbins, density=True)
        plt.xlabel('number of events', fontsize=12)
        plt.title(label)

    plt.tight_layout()


def plot_gaus_domain(tmun, tmu, q0, nbins=50, df=3, title=''):
    """Verify that test-statistic distributions follow the expected chi-squared.

    Parameters
    ----------
    tmun  : array-like – t_μ values from full-vector comparison.
    tmu   : array-like – t_μ values for signal parameter only.
    q0    : array-like – q₀ = t_μ(μ=0) values.
    nbins : int
    df    : int – degrees of freedom for the chi-squared overlay on tmun.
    title : str
    """
    subplot = pltext.canvas(3, 3)

    subplot(1)
    _, xs, _ = pltext.hist(tmun[tmun >= 0], nbins, density=True)
    xcs = 0.5 * (xs[1:] + xs[:-1])
    plt.plot(xcs, stats.chi2(df).pdf(xcs), label=fr'$\chi^2$({df})')
    plt.xlabel(r'$t_\mu(x, n)$', fontsize=12)
    plt.title(title)
    plt.legend()

    subplot(2)
    _, xs, _ = pltext.hist(tmu[tmu >= 0], nbins, density=True)
    xcs = 0.5 * (xs[1:] + xs[:-1])
    plt.plot(xcs, stats.chi2(1).pdf(xcs), label=r'$\chi^2$(1)')
    plt.xlabel(r'$t_\mu(x, 1)$', fontsize=12)
    plt.title(title)
    plt.legend()

    subplot(3)
    pltext.hist(np.sqrt(q0[q0 > 0]), nbins)
    plt.xlabel(r'$Z_0(x)$', fontsize=14)

    plt.tight_layout()


# ---------------------------------------------------------------------------
# Feldman-Cousins confidence intervals
# ---------------------------------------------------------------------------

def plot_exps_fc_confint(dfs, cls=(0.68, 0.9)):
    """Plot FC confidence band on signal count vs. true signal count.

    Parameters
    ----------
    dfs : list[pandas.DataFrame] – one DataFrame per true signal hypothesis.
        Each must have columns ``nbb0`` (true signal), ``nbb`` (fitted signal),
        and ``tmu``.
    cls : tuple – confidence levels to display.
    """
    n0s   = np.array([np.mean  (df.nbb0[df.tmu >= 0]) for df in dfs])
    nns   = np.array([np.median(df.nbb [df.tmu >= 0]) for df in dfs])

    subplot = pltext.canvas(2, 2, 6, 8)
    subplot(1)
    plt.plot(nns, n0s)
    for cl in cls:
        ci = [confint.fca_segment(df.tmu[df.tmu >= 0],
                                  df.nbb[df.tmu >= 0], cl) for df in dfs]
        ci = ut.list_transpose(ci)
        plt.fill_betweenx(n0s, *ci, alpha=0.5, color='y',
                          label=f'FC CI {100*cl:.0f} % CL')
    plt.grid()
    plt.legend()
    plt.xlabel(r'$n_{\beta\beta}$',       fontsize=14)
    plt.ylabel(r'$n_{\beta\beta}$ true', fontsize=14)
    plt.tight_layout()


def plot_exps_z0(dfs):
    """Plot discovery significance Z₀ vs. true signal count (and half-life).

    Parameters
    ----------
    dfs : list[pandas.DataFrame] – same format as for :func:`plot_exps_fc_confint`.
        Must also contain columns ``tau0`` (true half-life) and ``q0``.
    """
    n0s   = np.array([np.mean  (df.nbb0[df.tmu >= 0]) for df in dfs])
    tau0s = np.array([np.mean  (df.tau0[df.tmu >= 0]) for df in dfs])
    q0s   = np.array([np.median(df.q0  [df.q0  >= 0]) for df in dfs])
    z0s   = np.sqrt(q0s)

    subplot = pltext.canvas(2, 2, 6, 8)

    for iplot, (yvals, ylabel) in enumerate(
            [(n0s,   r'$n_{\beta\beta}$ true'),
             (tau0s, r'$T_{\beta\beta}$')], start=1):

        subplot(iplot)
        plt.plot(z0s, yvals)

        for cl in (0.68, 0.9):
            ci = [confint.fca_segment(df.tmu[df.q0 >= 0],
                                      df.q0 [df.q0 >= 0], cl) for df in dfs]
            ci = [np.sqrt(c) for c in ut.list_transpose(ci)]
            plt.fill_betweenx(yvals, *ci, alpha=0.5, color='y',
                              label=f'FC CI {100*cl:.0f} % CL')

        plt.axvline(3, ls='--', color='k', label=r'$3\sigma$')
        plt.axvline(5, ls='--', color='r', label=r'$5\sigma$')
        plt.grid(which='both')
        plt.legend()
        plt.ylabel(ylabel, fontsize=14)
        plt.xlabel(r'$Z_0$', fontsize=14)
        if iplot == 2:
            plt.yscale('log')

    plt.tight_layout()


# ---------------------------------------------------------------------------
# Data vs. MC comparison
# ---------------------------------------------------------------------------

def plot_contributions(data, mcs, ns,
                        varname='E',
                        varrange=erange,
                        nbins=80,
                        ssamples=ssamples):
    """Histogram data and overlay scaled MC component histograms.

    Parameters
    ----------
    data     : pandas.DataFrame
    mcs      : list[pandas.DataFrame] – one per component.
    ns       : array-like – expected event counts used to scale each MC.
    varname  : str   – variable to plot.  Default ``'E'``.
    varrange : tuple – histogram range.  Default :data:`erange`.
    nbins    : int   – number of bins.
    ssamples : list[str] – component labels.
    """
    counts, edges = np.histogram(data[varname], nbins, range=varrange)
    centers  = 0.5 * (edges[1:] + edges[:-1])
    ecounts  = np.sqrt(counts)
    has_data = ecounts > 0

    plt.errorbar(centers[has_data], counts[has_data], yerr=ecounts[has_data],
                 marker='o', ls='', label='data')
    plt.ylabel('counts')

    totals = np.zeros(len(counts))
    for n, mc, label in zip(ns, mcs, ssamples):
        mc_counts, _ = np.histogram(mc[varname], edges)
        scaled        = n * mc_counts / np.sum(mc_counts)
        totals       += scaled
        plt.plot(centers, scaled, label=label)

    plt.plot(centers, totals, label='total')
    plt.grid()
    plt.title(varname)
    plt.legend()
