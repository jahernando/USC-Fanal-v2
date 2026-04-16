"""
crib_fc
=======
Helper functions for the Feldman-Cousins confidence interval notebook.

Each function encapsulates one conceptual block: plotting, table printing,
or computation used in `crib_feldman_cousins.ipynb`.
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import core.confint as confint


# ---------------------------------------------------------------------------
# 1. NOMAD problem
# ---------------------------------------------------------------------------

def plot_nomad_poisson(bkg=3.0, nmax=15):
    """Plot the Poisson distribution for background-only hypothesis
    and highlight n=0."""
    p_zero = stats.poisson.pmf(0, bkg)
    print('P(n=0 | mu=0, b={:.0f}) = {:.4f}  ({:.1f}%)'.format(bkg, p_zero, 100 * p_zero))

    ns = np.arange(0, nmax)
    plt.figure(figsize=(7, 4))
    plt.bar(ns, stats.poisson.pmf(ns, bkg), color='steelblue', alpha=0.7,
            label=r'Poisson($\mu=0, b={:.0f}$)'.format(bkg))
    plt.axvline(0, color='red', ls='--', lw=2, label='n=0 observed')
    plt.xlabel('n (observed events)')
    plt.ylabel('P(n)')
    plt.title('NOMAD problem: we expect b={:.0f} background, observe n=0'.format(bkg))
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()


# ---------------------------------------------------------------------------
# 2. Classical central interval
# ---------------------------------------------------------------------------

def classical_central_interval(mu, bkg, cl=0.90, nmax=30):
    """Classical central acceptance interval for signal mu with background bkg.

    Splits probability (1-cl)/2 in each tail.
    Returns (n_low, n_high).
    """
    alpha = 1 - cl
    ns = np.arange(0, nmax + 1)
    pmf = stats.poisson.pmf(ns, mu + bkg)
    cdf = np.cumsum(pmf)

    n_low = 0
    while n_low < nmax and cdf[n_low] < alpha / 2:
        n_low += 1

    n_high = 0
    while n_high < nmax and cdf[n_high] < 1 - alpha / 2:
        n_high += 1

    return n_low, n_high


def classical_central_band(mus, bkg, cl=0.90):
    """Build the classical central acceptance band over a grid of mu values.

    Returns (n_lows, n_highs) as numpy arrays.
    """
    n_lows, n_highs = [], []
    for mu in mus:
        nl, nh = classical_central_interval(mu, bkg, cl)
        n_lows.append(nl)
        n_highs.append(nh)
    return np.array(n_lows), np.array(n_highs)


def plot_classical_band(bkg=3.0, cl=0.90, n_obs=0, mu_max=15, nmu=300):
    """Plot the classical central acceptance band and check if n_obs falls inside."""
    mus = np.linspace(0, mu_max, nmu)
    n_lows, n_highs = classical_central_band(mus, bkg, cl)

    mask = (n_lows <= n_obs) & (n_highs >= n_obs)

    plt.figure(figsize=(8, 6))
    plt.fill_between(mus, n_lows, n_highs, alpha=0.3, color='steelblue',
                     label='Central acceptance band {:.0f}% CL'.format(100 * cl))
    plt.plot(mus, n_lows, 'b-', lw=1)
    plt.plot(mus, n_highs, 'b-', lw=1)
    plt.axhline(n_obs, color='red', ls='--', lw=2,
                label=r'$n_{{obs}} = {:d}$'.format(n_obs))

    if np.any(mask):
        mu_low = np.min(mus[mask])
        mu_high = np.max(mus[mask])
        plt.axvspan(mu_low, mu_high, alpha=0.3, color='red',
                    label=r'CI: [{:.1f}, {:.1f}]'.format(mu_low, mu_high))
        print('Classical central interval for n_obs={:d}: [{:.2f}, {:.2f}]'.format(
            n_obs, mu_low, mu_high))
    else:
        plt.text(5, 0.5, 'EMPTY INTERVAL\nfor $n_{obs}=0$', fontsize=14, color='red',
                 ha='center', fontweight='bold',
                 bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='red'))
        print('EMPTY INTERVAL! No mu >= 0 whose central acceptance band '
              'includes n={:d}'.format(n_obs))

    plt.xlabel(r'$\mu$ (true signal)')
    plt.ylabel(r'$n$ (observed events)')
    plt.title(r'Classical central acceptance band ({:.0f}% CL, $b={:.0f}$)'.format(
        100 * cl, bkg))
    plt.legend(loc='upper left')
    plt.grid(alpha=0.3)
    plt.ylim(-0.5, 20)
    plt.tight_layout()

    return mus, n_lows, n_highs


# ---------------------------------------------------------------------------
# 3. FC ordering example
# ---------------------------------------------------------------------------

def print_fc_ordering_table(mu=1.0, bkg=3.0, nmax=15):
    """Print the FC ordering table: n, mu_best, P(n|mu+b), P(n|best), R, t."""
    ns = np.arange(0, nmax)
    mu_best = np.maximum(ns - bkg, 0.)
    p_mu = stats.poisson.pmf(ns, mu + bkg)
    p_best = stats.poisson.pmf(ns, mu_best + bkg)
    R = p_mu / p_best
    t = -2 * (np.log(p_mu + 1e-300) - np.log(p_best + 1e-300))

    print('  n  | mu_best |  P(n|mu+b) | P(n|best)  |    R     |    t')
    print('-----|---------|------------|------------|----------|--------')
    for i, n in enumerate(ns):
        print(' {:2d}  |  {:4.1f}   | {:.6f}  | {:.6f}  | {:.4f}  | {:6.3f}'.format(
            n, mu_best[i], p_mu[i], p_best[i], R[i], t[i]))

    return ns, mu_best, p_mu, p_best, R, t


def print_fc_inclusion_order(mu=1.0, bkg=3.0, cl=0.90, nmax=15):
    """Print the step-by-step FC inclusion order and return the acceptance interval."""
    ns = np.arange(0, nmax)
    mu_best = np.maximum(ns - bkg, 0.)
    p_mu = stats.poisson.pmf(ns, mu + bkg)
    p_best = stats.poisson.pmf(ns, mu_best + bkg)
    t = -2 * (np.log(p_mu + 1e-300) - np.log(p_best + 1e-300))

    order = np.argsort(t)
    cum_prob = 0.
    included = []

    print('FC inclusion order (mu={:.1f}, b={:.1f}, CL={:.0f}%):'.format(mu, bkg, 100 * cl))
    print('  step |  n  |    t   |  P(n|mu+b)  |  cumulative P | included')
    print('-------|-----|--------|-------------|---------------|----------')
    for step, idx in enumerate(order):
        cum_prob += p_mu[idx]
        if len(included) == 0 or (cum_prob - p_mu[idx]) < cl:
            included.append(ns[idx])
            mark = '<--'
        else:
            included.append(ns[idx])
            mark = '<-- (reaches CL)'

        print('   {:2d}  | {:2d}  | {:5.3f} |  {:.6f}   |   {:.6f}     | {:s}'.format(
            step, ns[idx], t[idx], p_mu[idx], cum_prob, mark))

        if cum_prob >= cl:
            break

    n_min, n_max = min(included), max(included)
    print()
    print('FC acceptance interval: [{:d}, {:d}]'.format(n_min, n_max))
    print('Cumulative probability: {:.4f} >= {:.2f}'.format(cum_prob, cl))

    return n_min, n_max


# ---------------------------------------------------------------------------
# 4. FC confidence belt
# ---------------------------------------------------------------------------

def plot_fc_belt(bkg=3.0, cl=0.90, mu_max=15., nmu=300):
    """Build and plot the FC confidence belt. Print CI for several n_obs values.

    Returns (mus, n0s, n1s, fc_ci) where fc_ci is the CI closure.
    """
    mus = np.linspace(0., mu_max, nmu)
    n0s, n1s = confint.fc_confband(mus, bkg, cl=cl)
    fc_ci = confint.get_fc_confinterval(mus, bkg, cl=cl)

    plt.figure(figsize=(8, 6))
    plt.fill_between(mus, n0s, n1s, alpha=0.3, color='green',
                     label='FC band {:.0f}% CL'.format(100 * cl))
    plt.plot(mus, n0s, 'g-', lw=1)
    plt.plot(mus, n1s, 'g-', lw=1)

    plt.axhline(0, color='red', ls='--', lw=2, label='$n_{obs} = 0$')
    ci_0 = fc_ci(0)
    plt.axvspan(ci_0[0], ci_0[1], alpha=0.2, color='red')

    print('FC {:.0f}% CI for n_obs=0, b={:.0f}: mu in [{:.2f}, {:.2f}]'.format(
        100 * cl, bkg, *ci_0))
    print('  -> Upper limit: mu < {:.2f}'.format(ci_0[1]))
    for n_test in [0, 1, 2, 3, 5, 8]:
        ci = fc_ci(n_test)
        print('  n_obs = {:d}: mu in [{:.2f}, {:.2f}]'.format(n_test, *ci))

    plt.xlabel(r'$\mu$ (true signal)')
    plt.ylabel(r'$n$ (observed events)')
    plt.title(r'Feldman-Cousins confidence belt ({:.0f}% CL, $b={:.0f}$)'.format(
        100 * cl, bkg))
    plt.legend(loc='upper left')
    plt.grid(alpha=0.3)
    plt.ylim(-0.5, 20)
    plt.tight_layout()

    return mus, n0s, n1s, fc_ci


# ---------------------------------------------------------------------------
# 5. Comparison: classical central vs FC
# ---------------------------------------------------------------------------

def plot_comparison_classical_fc(bkg=3.0, cl=0.90, mu_max=15., nmu=300):
    """Side-by-side comparison of classical central and FC acceptance bands."""
    mus = np.linspace(0., mu_max, nmu)
    n_lows_c, n_highs_c = classical_central_band(mus, bkg, cl)
    n0s_fc, n1s_fc = confint.fc_confband(mus, bkg, cl=cl)
    fc_ci = confint.get_fc_confinterval(mus, bkg, cl=cl)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left: classical central
    ax = axes[0]
    ax.fill_between(mus, n_lows_c, n_highs_c, alpha=0.3, color='steelblue',
                    label='Classical central')
    ax.plot(mus, n_lows_c, 'b-', lw=1)
    ax.plot(mus, n_highs_c, 'b-', lw=1)
    ax.axhline(0, color='red', ls='--', lw=2, label='$n_{obs} = 0$')
    ax.set_xlabel(r'$\mu$'); ax.set_ylabel(r'$n$')
    ax.set_title('Classical central ({:.0f}% CL)'.format(100 * cl))
    ax.set_ylim(-0.5, 20); ax.legend(loc='upper left'); ax.grid(alpha=0.3)
    ax.text(5, 1, 'EMPTY\nINTERVAL', fontsize=12, color='red',
            ha='center', fontweight='bold')

    # Right: FC
    ax = axes[1]
    ax.fill_between(mus, n0s_fc, n1s_fc, alpha=0.3, color='green',
                    label='Feldman-Cousins')
    ax.plot(mus, n0s_fc, 'g-', lw=1)
    ax.plot(mus, n1s_fc, 'g-', lw=1)
    ax.axhline(0, color='red', ls='--', lw=2, label='$n_{obs} = 0$')
    ci_0 = fc_ci(0)
    ax.axvspan(ci_0[0], ci_0[1], alpha=0.2, color='red',
               label=r'CI: $\mu \in [{:.1f}, {:.1f}]$'.format(*ci_0))
    ax.set_xlabel(r'$\mu$'); ax.set_ylabel(r'$n$')
    ax.set_title('Feldman-Cousins ({:.0f}% CL)'.format(100 * cl))
    ax.set_ylim(-0.5, 20); ax.legend(loc='upper left'); ax.grid(alpha=0.3)

    plt.suptitle(r'NOMAD case: $b = {:.0f}$, $n_{{obs}} = 0$'.format(bkg),
                 fontsize=14, y=1.02)
    plt.tight_layout()


# ---------------------------------------------------------------------------
# 6. Upper limit ↔ two-sided transition
# ---------------------------------------------------------------------------

def plot_fc_transition(bkg=3.0, cl=0.90, mu_max=20., nmu=400, n_obs_max=16):
    """Plot the automatic UL → two-sided transition of FC intervals."""
    mus = np.linspace(0., mu_max, nmu)
    fc_ci = confint.get_fc_confinterval(mus, bkg, cl=cl)

    ns_demo = np.arange(0, n_obs_max)
    cis = [fc_ci(n) for n in ns_demo]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: horizontal bars
    ax = axes[0]
    for i, n in enumerate(ns_demo):
        ci = cis[i]
        color = 'darkorange' if ci[0] == 0 else 'seagreen'
        tipo = 'UL' if ci[0] == 0 else 'two-sided'
        ax.barh(n, ci[1] - ci[0], left=ci[0], height=0.6, color=color, alpha=0.7,
                edgecolor='k', lw=0.5)
        ax.text(ci[1] + 0.3, n, tipo, va='center', fontsize=8, color=color)

    ax.axhline(bkg, color='gray', ls=':', lw=1.5,
               label='$b = {:.0f}$'.format(bkg))
    ax.set_ylabel(r'$n_{obs}$')
    ax.set_xlabel(r'$\mu$ (signal)')
    ax.set_title('FC {:.0f}% CI for each $n_{{obs}}$'.format(100 * cl))
    ax.legend()
    ax.grid(alpha=0.3)

    # Right: interval bounds
    ax2 = axes[1]
    mu_lows = [ci[0] for ci in cis]
    mu_highs = [ci[1] for ci in cis]
    ax2.plot(ns_demo, mu_lows, 'o-', color='seagreen', label=r'$\mu_{low}$')
    ax2.plot(ns_demo, mu_highs, 's-', color='darkorange', label=r'$\mu_{high}$')
    ax2.axvline(bkg, color='gray', ls=':', lw=1.5,
                label='$b = {:.0f}$'.format(bkg))
    ax2.set_xlabel(r'$n_{obs}$')
    ax2.set_ylabel(r'$\mu$')
    ax2.set_title(r'CI bounds as a function of $n_{obs}$')
    ax2.legend()
    ax2.grid(alpha=0.3)

    plt.tight_layout()

    print('Transition from upper limit to two-sided interval:')
    for n in ns_demo:
        ci = fc_ci(n)
        tipo = 'Upper limit' if ci[0] == 0 else 'Two-sided'
        print('  n={:2d}: mu in [{:5.2f}, {:5.2f}]  ({:s})'.format(
            n, ci[0], ci[1], tipo))


# ---------------------------------------------------------------------------
# 7. Effect of background on the FC belt
# ---------------------------------------------------------------------------

def plot_fc_vs_background(bkgs=(0., 3., 10.), cl=0.90, mu_max=15., nmu=300):
    """Compare FC bands for different background levels."""
    fig, axes = plt.subplots(1, len(bkgs), figsize=(16, 5))
    colors = ['royalblue', 'seagreen', 'darkorange']
    mus = np.linspace(0., mu_max, nmu)

    for ax, b, col in zip(axes, bkgs, colors):
        fc = confint.get_fc_confinterval(mus, b, cl=cl)
        n0, n1 = confint.fc_confband(mus, b, cl=cl)

        ax.fill_between(mus, n0, n1, alpha=0.3, color=col)
        ax.plot(mus, n0, '-', color=col, lw=1)
        ax.plot(mus, n1, '-', color=col, lw=1)

        ax.axhline(b, color='gray', ls=':', label='$n = b = {:.0f}$'.format(b))

        ci = fc(0)
        ax.axhline(0, color='red', ls='--', lw=1.5)
        ax.set_title('$b = {:.0f}$,  FC {:.0f}% CI($n=0$) = [{:.1f}, {:.1f}]'.format(
            b, 100 * cl, *ci), fontsize=10)
        ax.set_xlabel(r'$\mu$')
        ax.set_ylabel(r'$n$')
        ax.set_ylim(-0.5, 25)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

    plt.suptitle('Feldman-Cousins bands at {:.0f}% CL for different backgrounds'.format(
        100 * cl), fontsize=13)
    plt.tight_layout()


# ---------------------------------------------------------------------------
# 8. Comparison of confidence levels
# ---------------------------------------------------------------------------

def plot_fc_vs_cl(bkg=3.0, cls=(0.95, 0.90, 0.68), mu_max=15., nmu=300):
    """Compare FC intervals at different confidence levels."""
    mus = np.linspace(0., mu_max, nmu)
    ns_range = np.arange(0, 20)
    alphas = [0.2, 0.35, 0.5]
    colors = ['gold', 'orange', 'darkorange']

    plt.figure(figsize=(8, 6))
    for cl, alpha, col in zip(cls, alphas, colors):
        fc = confint.get_fc_confinterval(mus, bkg, cl=cl)
        ys = fc(ns_range)
        plt.fill_between(ns_range, ys[0], ys[1], alpha=alpha, color=col,
                         label='FC {:.0f}% CL'.format(100 * cl))

    plt.axvline(bkg, color='gray', ls=':', label='$b = {:.0f}$'.format(bkg))
    plt.xlabel(r'$n_{obs}$ (observed events)')
    plt.ylabel(r'$\mu$ (signal)')
    plt.title(r'FC intervals for $b = {:.0f}$ at different CLs'.format(bkg))
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()

    print('\nSummary for n_obs = 0, b = {:.0f}:'.format(bkg))
    for cl in cls:
        fc = confint.get_fc_confinterval(mus, bkg, cl=cl)
        ci = fc(0)
        print('  {:2.0f}% CL: mu in [{:.2f}, {:.2f}]  (UL = {:.2f})'.format(
            100 * cl, ci[0], ci[1], ci[1]))
