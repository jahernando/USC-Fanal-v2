#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:28:02 2022

@author: hernando
"""

import numpy as np
import core.utils as ut
import scipy.stats    as stats

nsigma = 5


def fc_confsegment(nu, bkg, cl = 0.68, nrange = None):
    """
    
    Compute Feldman-Cousing segment for nu, signal, and bkg, events
    
    Parameters
    ----------
    nu  : float, number of signal events
    bkg : float, number of bkg events
    cl  : float, confidence level.The default is 0.68.
    nrange : (int, int), range of expected number of events, default None (auto-defined)

    Returns
    -------
    int : (int, int), range of possibe number of events at CL

    """
    

    if (nrange == None):    
        nmax   = bkg + nu + nsigma * np.sqrt(bkg + nu)
        nrange = (0, nmax + 1)
     
    ns      = np.arange(*nrange)
    nuhats  = ns - bkg
    nuhats[nuhats <= 0] = 0
    ps     = stats.poisson.pmf(ns, bkg + nu)
    psbest = stats.poisson.pmf(ns, bkg + nuhats)
    ts = -2 * (np.log(ps) - np.log(psbest))
    vals = sorted(zip(ts, ps, ns))
    _, ops, ons = ut.list_transpose(vals)
    cops = np.cumsum(ops)
    assert (cops[-1] > cl), 'not enough range to compute CL'
    i = 0
    while (cops[i] < cl): i += 1
    int = np.min(ons[:i+1]), np.max(ons[:i+1])
    return int
    

def fc_confband(nus, bkg, cl = 0.68, nrange = None):
    """
    
    Parameters
    ----------
    nus    : np.array(float), array with the scan on number of signal values
    bkg    : float, number of bkg events
    cl     : float, confidence level. The default is 0.68.
    nrange : (int, int), range of expected number of events, default None (auto-defined)


    Returns
    -------
    n0s    : np.array(int), lower number of events of the CL band 
    n1s    : np.array(int), upper number of events of the CL band

    """
    
    vals   = [fc_confsegment(nu, bkg, cl, nrange) for nu in  nus] 
    n0s, n1s =  ut.list_transpose(vals)
    return np.array(n0s, int), np.array(n1s, int)


def get_fc_confinterval(nus, bkg, cl = 0.68, nrange = None):
    """
    
    return a function to comppute the FC confidence intervals for a given observation

    Parameters
    ----------
    nus    : np.array(float), list of possible mu values
    bkg    : float, value of the background
    cl     : float, confidence value. The default is 0.68.
    nrange : (float, float) or None, range of the possible observation values.
    The default is None.

    Returns
    -------
    ci     : function that computes the FC CI for a given CL

    """
    
    n0s, n1s = fc_confband(nus, bkg, cl, nrange)
    
    def _ci(nobs):
        """        
        return cover interval at cl for number of observed events, *nobs*

        """
        if (isinstance(nobs, np.ndarray)):
            ys = [_ci(ni) for ni in nobs]
            ys = ut.list_transpose(ys)
            return np.array(ys)
        nu1 = np.max(nus[n0s <= nobs])
        nu0 = np.min(nus[n1s >= nobs])
        return np.array((nu0, nu1))
    
    return _ci
        
def fca_segment(tmus, ns, cl = 0.9):
    """
    
    Return the FC segment of a list of observations (ns) with they FC ordering variable (tmus)
    at a givel CL (cl). 

    Parameters
    ----------
    tmus : np.array(float), values of the FC ordering varialbe of the observations
    ns   : np.array(float), values of the observations
    cl   : float, confidence level. The default is 0.9.

    Returns
    -------
    cint : (float, float), tuple with the minimum and maximum of the value sof the FC segment

    """
    vals = zip(tmus, ns)
    vals = sorted(vals)
    _, ons = ut.list_transpose(vals)
    xpos = cl * len(tmus)
    ipos = int(xpos)
    ipos = ipos if xpos - ipos < 0.5 else ipos + 1
    cint = np.array( (np.min(ons[:ipos]), np.max(ons[:ipos])) )
    return cint