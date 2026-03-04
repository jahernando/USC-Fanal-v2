"""
core.efit
=========
Extended (and simultaneous) maximum-likelihood fitting with composite PDFs.

The module provides:

* :func:`best_estimate`          – generic parameter optimisation via -2 log L.
* :func:`llike_scan`             – profile-likelihood scan over one parameter.
* :func:`tmu`                    – test statistic t_μ for a hypothesis.
* :func:`tmu_conf_int`           – confidence interval from a t_μ scan.
* :class:`ComPDF`                – weighted mixture of PDFs.
* :class:`ExtComPDF`             – extended version (Poisson normalisation term).
* :class:`ConstrainedExtComPDF`  – extended fit with Gaussian-constrained nuisance params.
* :class:`SimulExtComPDF`        – simultaneous signal + control-region fit.

JAH, 25/1/2021
"""

import numpy          as np
import scipy.stats    as stats
import scipy.optimize as optimize


# ---------------------------------------------------------------------------
# Profile-likelihood utilities
# ---------------------------------------------------------------------------

def llike_scan(x, ell, muhat, mus, ipos=0):
    """Scan the profile log-likelihood over parameter *ipos* fixing other params.

    Parameters
    ----------
    x     : array-like – observed data.
    ell   : object     – must expose ``loglike(x, *pars)`` and
                         ``best_estimate(x, *pars, mask=...)``.
    muhat : array-like – best-fit parameter vector.
    mus   : array-like – values of parameter *ipos* to scan.
    ipos  : int        – index of the parameter to scan.  Default 0.

    Returns
    -------
    numpy.ndarray
        ``-2 * (ll(mu) - ll(muhat))`` for each value in *mus*.
    """
    assert isinstance(ipos, int) and 0 <= ipos < len(muhat), \
        f'ipos={ipos} is invalid; must be an integer in [0, {len(muhat) - 1}]'

    ll_best = ell.loglike(x, *muhat)

    # Mask to fix all parameters except index ipos during conditional fit
    mask       = [True] * len(muhat)
    mask[ipos] = False

    def _conditional_ll(ni):
        """Log-likelihood with parameter ipos fixed to ni."""
        nni       = np.copy(muhat)
        nni[ipos] = ni
        res = ell.best_estimate(x, *nni, mask=mask)
        return ell.loglike(x, *res.x)

    lls = np.array([_conditional_ll(ni) for ni in mus])
    return -2. * (lls - ll_best)


def tmu(x, ell, muhat, mu, ipos=0):
    """Test statistic t_μ = -2(log L(μ) - log L(μ̂)).

    Parameters
    ----------
    x     : array-like – data.
    ell   : object     – PDF with ``loglike`` method.
    muhat : array-like – best-fit (unconditional) parameters.
    mu    : float | array-like
        Hypothesised parameter value(s).  When *ipos* == -1 the full
        parameter vector is compared.
    ipos  : int
        Index of the parameter of interest.  Use -1 to compare full vectors.

    Returns
    -------
    float – value of t_μ.
    """
    if ipos == -1:
        ll_best = ell.loglike(x, *muhat)
        ll_val  = ell.loglike(x, *mu)
        return -2. * (ll_val - ll_best)

    return llike_scan(x, ell, muhat, (mu,), ipos)[0]


def tmu_conf_int(mus, tmus, pvalue=0.68, df=1):
    """Confidence interval from a t_μ scan using chi-squared critical value.

    Parameters
    ----------
    mus    : numpy.ndarray – scanned parameter values.
    tmus   : numpy.ndarray – corresponding t_μ values.
    pvalue : float         – confidence level.  Default 0.68 (≈ 1σ).
    df     : int           – degrees of freedom for the chi-squared CDF.

    Returns
    -------
    tuple(float, float) – (lower, upper) confidence bound.
    """
    t_crit = stats.chi2.ppf(pvalue, df)
    inside = tmus <= t_crit
    return float(np.min(mus[inside])), float(np.max(mus[inside]))


# ---------------------------------------------------------------------------
# Generic best-estimate optimiser
# ---------------------------------------------------------------------------

def best_estimate(ll, x, *par0, mask=None, bounds=None):
    """Minimise -2 log L to find best-fit parameters.

    Parameters
    ----------
    ll     : callable(x, *pars) → float – log-likelihood function.
    x      : array-like – data.
    *par0  : float – initial parameter values.
    mask   : list[bool] | None
        If given, only parameters where mask[i] is True are varied; the rest
        are fixed at their *par0* values.
    bounds : list[tuple] | None
        Parameter bounds ``[(lo, hi), ...]``.  Defaults to ``(0, +∞)`` for
        all parameters (non-negative counts).

    Returns
    -------
    scipy.optimize.OptimizeResult
        ``result.x`` contains the full best-fit parameter vector (including
        fixed parameters when a *mask* is supplied).
    """
    npar   = len(par0)
    fun    = lambda par: -2. * ll(x, *par)
    bounds = bounds if bounds is not None else npar * ((0., np.inf),)

    if mask is None:
        return optimize.minimize(fun, par0, bounds=bounds)

    # --- Conditional optimisation with some parameters fixed ---
    mask   = np.array(mask, dtype=bool)
    assert len(mask) == npar, \
        'mask length must equal the number of parameters'

    par_arr  = np.array(par0, dtype=float)
    free_p0  = par_arr[mask]
    free_bnd = [b for b, m in zip(bounds, mask) if m]

    def masked_fun(free_pars):
        par_arr[mask] = free_pars
        return fun(par_arr)

    res = optimize.minimize(masked_fun, free_p0, bounds=free_bnd)

    # Reconstruct full parameter vector
    full_best          = np.copy(par_arr)
    full_best[mask]    = res.x
    res.x              = full_best
    return res


# ---------------------------------------------------------------------------
# Composite PDF
# ---------------------------------------------------------------------------

class ComPDF:
    """Weighted mixture of probability density functions.

    Each component PDF enters with a weight proportional to its expected
    number of events.  Internally, weights are always normalised to sum to 1.

    Parameters
    ----------
    pdfs : list
        PDF objects, each exposing ``.pdf(x)``, ``.logpdf(x)``, and
        ``.rvs(size)`` methods (compatible with :mod:`scipy.stats`).
    *ns  : float
        Expected number of events (weights) for each PDF.  If omitted,
        equal weights are used.
    """

    def __init__(self, pdfs, *ns):
        self.pdfs = pdfs
        self.ns   = np.array(ns, dtype=float) if ns else np.ones(len(pdfs))
        assert len(self.ns) == len(self.pdfs), \
            'number of PDFs and number of weights must match'

    def weights(self, *ns):
        """Normalise *ns* to fractions and return the total.

        Parameters
        ----------
        *ns : float
            Number of events per PDF.  If empty, ``self.ns`` is used.

        Returns
        -------
        ntot : float – sum of *ns*.
        fis  : numpy.ndarray – fractional weight of each component.
        """
        ns   = np.array(ns) if ns else self.ns
        assert len(ns) == len(self.pdfs), \
            'number of weights must equal number of PDFs'
        ntot = np.sum(ns)
        fis  = ns / ntot
        return ntot, fis

    def pdf(self, x, *ns):
        """Evaluate the mixture PDF at *x*.

        Parameters
        ----------
        x   : array-like
        *ns : float
            Component weights.  Defaults to ``self.ns``.

        Returns
        -------
        numpy.ndarray | float
        """
        ntot, fis = self.weights(*ns)
        if np.any(fis < 0.):
            return np.inf
        p = np.sum([fi * ipdf.pdf(x) for fi, ipdf in zip(fis, self.pdfs)], axis=0)
        return p

    def loglike(self, x, *ns):
        """Log-likelihood  Σ log p(xᵢ | ns).

        Parameters
        ----------
        x   : array-like – data.
        *ns : float – component weights (default ``self.ns``).

        Returns
        -------
        float
        """
        return float(np.sum(np.log(self.pdf(x, *ns))))

    def best_estimate(self, x, *ns, mask=None):
        """Minimise -2 log L to find the best-fit *ns*.

        Parameters
        ----------
        x    : array-like – data.
        *ns  : float – initial parameter values.
        mask : list[bool] | None – parameters to optimise (True = free).

        Returns
        -------
        scipy.optimize.OptimizeResult
        """
        return best_estimate(self.loglike, x, *ns, mask=mask)

    def rvs(self, *ns, size=1):
        """Generate *size* random variates from the mixture.

        Parameters
        ----------
        *ns  : float – component weights (default ``self.ns``).
        size : int   – number of samples.

        Returns
        -------
        numpy.ndarray of shape (*size*,)
        """
        ntot, fis = self.weights(*ns)
        # Cumulative fractions used as category boundaries
        cum_fis   = np.cumsum(fis)
        u         = np.random.uniform(size=size)
        # Assign each uniform draw to a component
        components = np.digitize(u, cum_fis)
        x = np.array([self.pdfs[i].rvs() for i in components])
        return x


# ---------------------------------------------------------------------------
# Extended (Poisson-normalised) composite PDF
# ---------------------------------------------------------------------------

class ExtComPDF(ComPDF):
    """Extended maximum-likelihood composite PDF.

    Adds a Poisson term ``log P(n_obs | Σ nᵢ)`` to the log-likelihood so that
    the total number of events is also a fitted quantity.

    Parameters
    ----------
    pdfs : list – see :class:`ComPDF`.
    *ns  : float – expected number of events per component.
    """

    def __init__(self, pdfs, *ns):
        ComPDF.__init__(self, pdfs, *ns)

    def loglike(self, x, *ns):
        """Extended log-likelihood = Σ log p(xᵢ) + log Poisson(n_obs | Σnᵢ).

        Parameters
        ----------
        x   : array-like – data.
        *ns : float – expected counts per component.

        Returns
        -------
        float
        """
        ll_shape = ComPDF.loglike(self, x, *ns)
        n_obs    = len(x)
        mu_tot, _ = self.weights(*ns)
        ll_norm  = stats.poisson.logpmf(n_obs, mu_tot)
        return ll_shape + ll_norm

    def rvs(self, *ns, size=1):
        """Generate *size* experiments, each with a Poisson-fluctuated number of events.

        Parameters
        ----------
        *ns  : float – expected counts per component.
        size : int   – number of independent experiments.

        Returns
        -------
        numpy.ndarray (if size == 1) or list of numpy.ndarray
        """
        mu_tot, _ = self.weights(*ns)

        def _single_experiment():
            n = int(stats.poisson.rvs(mu=mu_tot, size=1))
            return ComPDF.rvs(self, *ns, size=n)

        xs = [_single_experiment() for _ in range(size)]
        return xs[0] if size == 1 else xs


# ---------------------------------------------------------------------------
# Constrained extended composite PDF
# ---------------------------------------------------------------------------

class ConstrainedExtComPDF(ExtComPDF):
    """Extended fit with Gaussian constraints on nuisance parameters.

    Adds penalty terms ``log N(nᵢ | n₀ᵢ, σᵢ)`` for each component with a
    non-zero uncertainty *σᵢ*.

    Parameters
    ----------
    pdfs : list – component PDFs.
    ns   : array-like – nominal expected event counts.
    uns  : array-like – 1-sigma uncertainties on each *ns*.
    """

    def __init__(self, pdfs, ns, uns):
        ComPDF.__init__(self, pdfs, *ns)
        self.ns  = np.array(ns,  dtype=float)
        self.uns = np.array(uns, dtype=float)

    def loglike(self, x, *ns):
        """Extended log-likelihood with Gaussian penalty terms."""
        ll = ExtComPDF.loglike(self, x, *ns)
        # Add Gaussian constraints only for components with non-zero uncertainty
        penalty = sum(
            stats.norm.logpdf(n, n0, u)
            for n, n0, u in zip(ns, self.ns, self.uns)
            if u != 0.
        )
        return ll + penalty


# ---------------------------------------------------------------------------
# Simultaneous signal + control region fit
# ---------------------------------------------------------------------------

class SimulExtComPDF:
    """Simultaneous extended likelihood over a signal and a control region.

    The number of events in the control region is linked to the signal region
    via fixed *ratios*: ``n_control = n_signal * ratios``.

    Parameters
    ----------
    ell_signal  : ComPDF – composite PDF for the signal region.
    ell_control : ComPDF – composite PDF for the control region.
    ratios      : array-like
        Scale factors converting signal-region counts to control-region counts
        for each component.
    """

    def __init__(self, ell_signal, ell_control, ratios):
        self.ell_signal  = ell_signal
        self.ell_control = ell_control
        self.ratios      = np.array(ratios, dtype=float)

    def ns_control(self, ns_signal):
        """Return expected counts in the control region.

        Parameters
        ----------
        ns_signal : array-like

        Returns
        -------
        numpy.ndarray
        """
        return np.array(ns_signal) * self.ratios

    def loglike(self, data, *ns_signal):
        """Combined log-likelihood of signal and control regions.

        Parameters
        ----------
        data      : tuple(array-like, array-like)
            ``data[0]`` = signal-region data, ``data[1]`` = control-region data.
        *ns_signal : float – expected signal-region counts per component.

        Returns
        -------
        float
        """
        signal, control = data[0], data[1]
        ll_sig = self.ell_signal .loglike(signal,  *ns_signal)
        ll_con = self.ell_control.loglike(control, *self.ns_control(ns_signal))
        return ll_sig + ll_con

    def best_estimate(self, data, *ns_signal, mask=None):
        """Minimise -2 log L to find the best-fit signal-region counts.

        Parameters
        ----------
        data       : tuple(array-like, array-like) – signal and control data.
        *ns_signal : float – initial parameter values.
        mask       : list[bool] | None

        Returns
        -------
        scipy.optimize.OptimizeResult
        """
        return best_estimate(self.loglike, data, *ns_signal, mask=mask)
