#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 11:28:02 2022

@author: hernando
"""

def fc_confint_segment(nu, bkg, nrange, cl = 0.68):
    """
    

    Parameters
    ----------
    nu  : float, number of signal events
    bkg : float, number of bkg events
    ns  : np.array(int), arry of possible number of events (integers)
    cl  : float, confidence level.The default is 0.68.

    Returns
    -------
    int : (int, int), range of possibe number of events at CL

    """
    
    ns      = np.arange(*nrange)
    nuhats  = ns - bkg
    nuhats[nuhats <= 0] = 0
    ps     = stats.poisson.pmf(ns, bkg + nu)
    psbest = stats.poisson.pmf(ns, bkg + nuhats)
    ts = -2 * (np.log(ps) - np.log(psbest))
    vals = sorted(zip(ts, ps, ns))
    _, ops, ons = ut.list_transpose(vals)
    cops = np.cumsum(ops)
    i = 0
    while (cops[i] < cl): i += 1
    int = np.min(ons[:i+1]), np.max(ons[:i+1])
    return int
    

def fc_confint_band(nus, bkg, nrange, cl = 0.68):
    """
    
    Parameters
    ----------
    nus    : np.array(float), array with the scan on number of signal values
    bkg    : float, number of bkg events
    nrange : (int, int), range of expected number of events
    cl     : float, confidence level. The default is 0.68.

    Returns
    -------
    n0     : np.array(int), lower number of events of the CL band 
    n1     : np.array(int), upper number of events of the CL band

    """
    
    vals   = [fc_confint_segment(nu, bkg, nrnage, cl) for nu in  nus] 
    n0, n1 =  ut.list_transpose(vals)
    return n0, n1
