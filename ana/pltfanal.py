#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 09:22:52 2022

@author: hernando
"""

import numpy  as np
#import pandas as pd

#import scipy.constants as constants
import scipy.stats     as stats
#import scipy.optimize  as optimize

import core.utils as ut
#import core.efit  as efit
import core.confint as confint

#import ana.fanal  as fn


import core.pltext as pltext
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


ssamples = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']
nbins    = 100
erange   = (2.4, 2.7)


def plot_fit_ell(x,
                 par, 
                 ell, 
                 bins = nbins,
                 parnames = ssamples,
                 plot_residuals = True,
                 title = ''):
    """ plot the data x, and superimposed the pdf with parameters (par)
    inputs:
        x    : np.array(float), data to plot
        par  : tuple(float), parameters of the pdf
        pdf  : function, pdf(x, *par), the pdf values of the distribution along x
        pdfs : tuple(pdfs), the mc pdfs for each sample
        parnames: tuple(str), list of the parameters (and samples) for the legend
        title   : str, title of the plot, default = ''
    """

    subplot = pltext.canvas(1, 1, 8, 10)
    subplot(1)

    counts, edges = np.histogram(x, bins);
    centers = 0.5 * (edges[1:] + edges[:-1])
    ecounts = np.sqrt(counts)
    sel     = ecounts > 0
    nn      = np.sum(par)
    factor  = nn * (centers[1] - centers[0])
    plt.errorbar(centers[sel], counts[sel], yerr = ecounts[sel], 
                 marker = 'o', ls = '', label = 'data')

    label  = 'ELL fit \n'
    plt.plot(centers, factor * ell.pdf(centers, *par), label = label)

    i = 0
    for ni, ipdf in zip(par, ell.pdfs):
        factor = ni * (centers[1] - centers[0])
        label  = ' {:s} : {:6.2f} \n'.format(parnames[i], ni)
        plt.plot(centers, factor * ipdf.pdf(centers), label = label)
        i += 1
    plt.legend(); plt.grid();
    plt.title(title)

    if (not plot_residuals): return
    
    ax = plt.gca()
    divider = make_axes_locatable(ax)
    ax2 = divider.append_axes("bottom", size = '20%', pad = 0)
    ax.figure.add_axes(ax2)
    fun = lambda x, *p : factor * ell.pdf(x, *p)
    pltext.hresiduals(x, bins, fun, par)

    return

def plot_fit_simell(values,
                    pars, 
                    simell, 
                    **kargs):
    
    signal, control = values[0], values[1]
     
    plot_fit_ell(signal , pars, simell.ell_signal, title = 'signal data ', **kargs)
    plot_fit_ell(control, pars * simell.ratios,
                 simell.ell_control, title = 'control data', **kargs)
    return


def plot_tmu_scan(nis, tmus, cls = (0.68, 0.9), titles = ssamples):
    
    npars = len(nis)
    subplot = pltext.canvas(npars, npars)
    for i in range(npars):
        ni, tmu = nis[i], tmus[i]
        subplot(i + 1)
        plt.xlabel('number of events'); plt.ylabel(r'$\Delta -2 \mathrm{log} \mathcal{L}$')
        plt.plot(ni, tmu);
        for cl in cls:
            t0 = stats.chi2.ppf(cl, 1)
            plt.plot(ni, t0 * np.ones(len(ni)), '-', label = 'CL {:4.0f}'.format(100*cl));
        plt.grid(); plt.title(titles[i])
    plt.tight_layout()
    
    

#  Plot Experiments

def plot_nevts(nevts, nbins = 50, labels = ssamples):
    
    npars = len(nevts)
    
    subplot = pltext.canvas(npars, npars)
    
    for i in range(npars):
        subplot( i + 1)
        pltext.hist(nevts[i], nbins, density = True)
        plt.xlabel(r'number of events', fontsize = 12)
        plt.title(labels[i])

    plt.tight_layout()
    
    
def plot_gaus_domain(tmun, tmu, q0, nbins = 50, df = 3, title = ''):
    
    subplot = pltext.canvas(3, 3)
    
    subplot(1)
    _, xs, _ = pltext.hist(tmun[tmun >= 0], nbins, density = True);
    xcs = 0.5*(xs[1:] + xs[:-1])
    plt.plot(xcs, stats.chi2(df).pdf(xcs))
    plt.xlabel(r'$t_\mu(x, n)$', fontsize = 12)
    plt.title(title)

    #left, bottom, width, height = [0.2, 0.2, 0.5, 0.5]
    #ax2 = plt.gcf().add_axes([left, bottom, width, height])
    #pltext.hist(stats.chi2(df).cdf(tmun), nbins);
    #ax2.set_xlabel('p value')

    subplot(2)
    _, xs, _ = pltext.hist(tmu[tmu >= 0], nbins, density = True)
    xcs = 0.5*(xs[1:] + xs[:-1])
    plt.plot(xcs, stats.chi2(1).pdf(xcs))
    plt.xlabel(r'$t_\mu(x, 1)$', fontsize = 12)
    plt.title(title)

    
    #left, bottom, width, height = [0.5 + 0.26, 0.4, 0.2, 0.3]
    #ax3 = plt.gcf().add_axes([left, bottom, width, height])
    #pltext.hist(stats.chi2(1).cdf(tmu), nbins);
    #plt.xlabel('p value')


    subplot(3)
    pltext.hist(np.sqrt(q0[q0>0]), nbins)
    plt.xlabel(r'$Z_0(x)$', fontsize = 14)

    plt.tight_layout()


def plot_exps_fc_confint(dfs, cls = (0.68, 0.9)):

    n0s   = np.array([np.mean  (df.nbb0[df.tmu >= 0]) for df in dfs]) 
    tau0s = np.array([np.mean  (df.tau0[df.tmu >= 0]) for df in dfs])
    nns   = np.array([np.median(df.nbb [df.tmu >= 0]) for df in dfs])

    subplot = pltext.canvas(2, 2, 6, 8)
    subplot(1)
    plt.plot(nns, n0s);
    for cl in cls:
        ci = [confint.fca_segment(df.tmu[df.tmu >= 0], df.nbb[df.tmu >= 0], cl) for df in dfs]
        ci = ut.list_transpose(ci)
        plt.fill_betweenx(n0s, *ci, alpha = 0.5, color = 'y', label = 'FC CI {:2.0f} % CL'.format(100*cl))
    plt.grid(); plt.legend();
    plt.xlabel(r'$n_{\beta\beta}$ ', fontsize = 14); 
    plt.ylabel(r'$n_{\beta\beta}$ true ', fontsize = 14);

    subplot(2)
    plt.plot(nns, tau0s);
    for cl in cls:
        ci = [confint.fca_segment(df.tmu[df.tmu >= 0], df.nbb[df.tmu >= 0], cl) for df in dfs]
        ci = ut.list_transpose(ci)
        plt.fill_betweenx(tau0s, *ci, alpha = 0.5, color = 'y', label = 'FC CI {:2.0f} % CL'.format(100*cl))
    plt.grid(which = 'both'); plt.legend();
    plt.xlabel(r'$n_{\beta\beta}$ ', fontsize = 14); 
    plt.ylabel(r'$T_{\beta\beta}$ ', fontsize = 14);
    #plt.xscale('log'); 
    plt.yscale('log'); 
    plt.tight_layout()
    
    
def plot_exps_z0(dfs):

    n0s   = np.array([np.mean  (df.nbb0[df.tmu >= 0]) for df in dfs])
    tau0s = np.array([np.mean  (df.tau0[df.tmu >= 0]) for df in dfs])
    q0s   = np.array([np.median(df.q0  [df.q0  >= 0]) for df in dfs])
    z0s   = np.sqrt(q0s)

    subplot = pltext.canvas(2, 2, 6, 8)
    subplot(1)
    plt.plot(z0s, n0s);
    for cl in (0.68, 0.9):
        ci = [confint.fca_segment(df.tmu[df.q0 >= 0], df.q0[df.q0 >= 0], cl) for df in dfs]
        ci = ut.list_transpose(ci)
        ci = [np.sqrt(c) for c in ci]
        plt.fill_betweenx(n0s, *ci, alpha = 0.5, color = 'y', label = 'FC CI {:2.0f} % CL'.format(100*cl))
    plt.plot(3 * np.ones(len(n0s)), n0s)
    plt.plot(5 * np.ones(len(n0s)), n0s)
    plt.grid(which = 'both'); plt.legend();
    plt.ylabel(r'$n_{\beta\beta}$ true', fontsize = 14); plt.xlabel(r'$Z_0$', fontsize = 14);

    subplot(2)
    plt.plot(z0s, tau0s);
    for cl in (0.68, 0.9):
        ci = [confint.fca_segment(df.tmu[df.q0 >= 0], df.q0[df.q0 >= 0], cl) for df in dfs]
        ci = ut.list_transpose(ci)
        ci = [np.sqrt(c) for c in ci]
        plt.fill_betweenx(tau0s, *ci, alpha = 0.5,  color = 'y', label = 'FC CI {:2.0f} % CL'.format(100*cl))
    plt.plot(3 * np.ones(len(n0s)), tau0s)
    plt.plot(5 * np.ones(len(n0s)), tau0s)
    plt.grid(which = 'both'); plt.legend();
    plt.ylabel(r'$T_{\beta\beta}$', fontsize = 14); plt.xlabel(r'$Z_0$', fontsize = 14);
    plt.yscale('log');
    
    plt.tight_layout()
    
    
def _plot_exps_ul(dfs, cls = (0.9, 0.95)):
    
    n0s   = np.array([np.mean  (df.nbb0[df.qmu >=  0]) for df in dfs]) 
    tau0s = np.array([np.mean  (df.tau0[df.tmu >= 0]) for df in dfs])
    nns   = np.array([np.median(df.nbb[df.qmu  >= 0]) for df in dfs])

    subplot = pltext.canvas(2, 2, 6, 8)
    subplot(1)
    plt.plot(nns, n0s);
    for cl in cls:
        ci = [confint.fca_segment(df.qmu[df.qmu >= 0], df.nbb[df.qmu >= 0], cl) for df in dfs]
        ci = ut.list_transpose(ci)
        plt.fill_betweenx(n0s, *ci, alpha = 0.2, label = 'Upper Limit {:2.0f} % CL'.format(100*cl))
        plt.plot(ci[0], n0s, alpha = 1., ls = '--', label = 'Upper Limit {:2.0f} % CL'.format(100*cl))
    plt.grid(which = 'both'); plt.legend();
    plt.ylabel(r'$n_{\beta\beta}$ true', fontsize = 14); plt.xlabel(r'$n_{\beta\beta}$', fontsize = 14);

    
    subplot(2)
    plt.plot(nns, tau0s);
    for cl in cls:
        ci   = [confint.fca_segment(df.qmu[df.qmu >= 0], df.nbb[df.qmu >= 0], cl) for df in dfs]
        ci   = ut.list_transpose(ci)
        ctau = [tau(np.array(c)) for c in ci]
        plt.fill_between(n0s, *ctau, alpha = 0.2, label = 'FC CI {:2.0f} % CL'.format(100*cl))
        plt.plot(n0s, ctau[0], alpha = 1., ls = '--', label = 'Upper Limit {:2.0f} % CL'.format(100*cl))
    plt.grid(which = 'both'); plt.legend();
    plt.xlabel(r'$n_{\beta\beta}$', fontsize = 14); plt.ylabel(r'$T_{\beta\beta}$', fontsize = 14);
    plt.yscale('log');
    
    plt.tight_layout()
    

#---
    
def plot_contributions(data,
                       mcs,
                       ns,
                       varname  = 'E',
                       varrange = erange,
                       nbins    = 80,
                       ssamples = ssamples):
    
    # plot data
    counts, bins = np.histogram(data[varname], nbins, range = varrange)
    cbins = 0.5 * (bins[1:] + bins[:-1])
    esel  = counts > 0
    ecounts = np.sqrt(counts)
    plt.errorbar(cbins[esel], counts[esel], yerr = ecounts[esel], marker = 'o', ls = '', label = 'data');
    plt.ylabel('counts')
        
    #ax = plt.gca().twinx()
    i = 0
    utots = np.zeros(len(counts))
    for n, mc in zip(ns, mcs):
        ucounts, _   = np.histogram(mc[varname], bins)
        ucounts      = n * ucounts/np.sum(ucounts)
        utots        += ucounts
        plt.plot(cbins, ucounts, label = ssamples[i]);
        i += 1
    plt.plot(cbins, utots, label = 'total');
    plt.grid(); plt.title(varname), plt.legend();
    return


#------

def plt_data_fanal(data, sel, mcs, nbkgs, nbins = 20):
    ranges = {'E' : eroi, 'blob2_E': eblob2_range}
    def _plot(varname):
        labels = [r"$^{214}$Bi", r"$^{108}$Tl"]
        counts, bins = np.histogram(data[sel][varname], nbins, range = ranges[varname])
        cbins = 0.5 * (bins[1:] + bins[:-1])
        esel  = counts > 0
        ecounts = np.sqrt(counts)
        plt.errorbar(cbins[esel], counts[esel], yerr = ecounts[esel], marker = 'o', ls = '', label = 'data');
        i = 0
        for n, mc in zip(nbkgs, mcs):
            ucounts, _   = np.histogram(mc[varname], bins)
            ucounts      = n * ucounts/np.sum(ucounts)
            plt.plot(cbins, ucounts, label = labels[i]);
            i+=1
        plt.grid(); plt.title(varname), plt.legend(); 
        plt.xlabel('Energy (MeV)'); plt.ylabel('counts')
        
    subplot = pltext.canvas(2)
    subplot(1); _plot('E')
    subplot(2); _plot('blob2_E')
