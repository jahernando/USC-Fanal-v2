import numpy  as np
import pandas as pd

from collections import namedtuple

import scipy.constants as constants
import scipy.stats     as stats
#import scipy.optimize  as optimize

import core.utils as ut
import core.efit  as efit

#import core.pltext as pltext
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable


ssamples = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']

erange     = (2.400, 2.650)
eroi       = (2.440, 2.475)
keys       = ['E', 'num_tracks', 'blob2_E', 'RoI']
varnames   = ['E', 'num_tracks', 'blob2_E', 'E']
varranges  = [erange, (1., 1.1), (0.4, np.inf), eroi]

#label   = 'track0_E'
erange   = (2.4, 2.7)
eroi     = (2.43,  2.48)
eblob2   = 0.4
#bins    = 100

NA       = constants.Avogadro
abundace = 0.9
Qbb      = 2458 # keV Qbb value
W        = 135.9


def half_life(nbb, exposure, eff, abundance = abundace, W = W):
    """  Compute the half-life time
    inputs:
        nbb       : float, number of signal events in RoI
        exposure  : float, (kg y)
        eff       : total signal efficiency, (in fraction)
        abundance : float, isotope fraction (0.9)
        W         : float, Atomic weight, (135.9 g/mol for 136Xe)
    returns:
        tau       : float, half-life (y)
    """
    tau  = 1e3 * eff * abundance * (exposure / nbb) * (NA / W) * np.log(2.)
    return tau


def selection_blind(df, eroi = eroi, eblob2 = eblob2):
    """
    
    returns the blind selection

    Parameters
    ----------
    df     : DataFrame, data
    eroi   : tuple(float, float), Energy Track0 blind range
    eblob2 : tuple(float, float), Energy range of the 2nd blob

    Returns
    -------
    sel    : np.array(bool), selection of blind events

    """
    sel0 = (df.track0_E >= eroi[0]) &  (df.track0_E < eroi[1])
    sel1 = (df.blob2_E  > eblob2)
    sel  = np.logical_or(sel0, sel1)
    return ~sel


def generate_mc_experiment(mcs, nevts):
    """ generate a MC experiment with the mcs samples, mcs, and the events in each sample, nevents
    inputs  :
        mcs     : tuple(DF), DFs with the mc samples
        nevents : tuple(int), number of events to select in each sample
    returns :
        mc      : DF, DF with the mix of number event of events, nevents, of the initial samples, mcs
    """
    def _ni(ni):
        return stats.poisson.rvs(ni, size = 1)[0]
    nns  = [_ni(ni) for ni in nevts]
    xmcs = [mc.sample(n = ni) for mc, ni in zip(mcs, nns)] # sample ni events of mc-sample mci
    mc   = pd.concat(xmcs) # concatenate all background to generate a mc sample that mimics the data-blind sample
    return mc



def get_ell(mcs, 
            refnames,
            refranges,
            varname  = 'E',
            varrange = erange,
            bins     = 100):
    """
    Construct a Composite PDF object
    

    Parameters
    ----------
    mcs       : tuple(DataFrames), data frames of the mc samples
    refnames  : tuple(str), list of the variables of the selection to create the pdfs
    refranges : tuple((float, float)), list of the ranges of the selection to create the pdfs
    varname   : str, name of the variable of the pdf, The default is 'E'.
    varrange  : (float, float), range of the pdf variable
    bins      : int, number of bins of the histograms to create the pdfs

    Returns
    -------
    ell      : object, Extended Maximum LL object to do a fit to combined pdfs
    """
    
    refmcs = [ut.selection_sample(mc, refnames, refranges) for mc in mcs]

    # generate the PDFs using the blind mc samples
    histos   = [np.histogram(mc[varname], bins, range = varrange) for mc    in refmcs]
    pdfs     = [stats.rv_histogram(histo)                         for histo in histos]
    ell      = efit.ExtComPDF(pdfs)
    
    return ell


def prepare_fit_ell(mcs, 
                    nevts,
                    varnames,
                    varranges,
                    refnames = [],
                    refranges = [],
                    varname  = 'E',
                    varrange = erange,
                    bins     = 100):
    """
    
    Parameters
    ----------
    mcs       : tuple(DataFrames), data frames of the mc samples
    varnames  : tuple(str), list of the variables of the selection
    varranges : tuple((float, float)), list of the ranges of the selection
    refnames  : tuple(str), list of the variables of the selection to create the pdfs.
        If empty, the same as varnames
    refranges : tuple((float, float)), list of the ranges of the selection to create the pdfs.
        If empty, the same as varranges
    varname   : str, name of the variable of the pdf, The default is 'E'.
    varrange  : (float, float), range of the pdf variable
    bins      : int, number of bins of the histograms to create the pdfs

    Returns
    -------
    fit       : functio(DataFrame) to fit data to the Composite PDF

    """
    
    refnames  = varnames  if refnames  == [] else refnames
    refranges = varranges if refranges == [] else refranges
        
    # expected number of events for each mc sample
    effs            = [ut.selection_efficiency(mc, varnames, varranges)[0] for mc in mcs]
    nevts_exp       = effs * np.array(nevts)
   
    
    ell = get_ell(mcs, refnames, refranges)
        
    def _fit(data):
            
        # select the data
        datana  = ut.selection_sample(data, varnames, varranges)

        # fit the energy values of the data 
        values  = datana[varname].values
        result  = ell.best_estimate(values, *nevts_exp)
            
        return result, values, ell, nevts_exp

    return _fit
            
    
def prepare_fit_simell(mcs, 
                       nevts,
                       varnames,
                       varranges,
                       refnames,
                       refranges,
                       connames,
                       conranges,
                       varname  = 'E',
                       varrange = erange,
                       bins     = 100):
    """
    
    Parameters
    ----------
    mcs       : tuple(DataFrames), data frames of the mc samples
    varnames  : tuple(str), list of the variables of the selection
    varranges : tuple((float, float)), list of the ranges of the selection
    refnames  : tuple(str), list of the variables of the selection to create the signal pdfs.
    refranges : tuple((float, float)), list of the ranges of the selection to create the signal pdfs.
    connames  : tuple(str), list of the variables of the selection to create the control pdfs.
    conranges : tuple((float, float)), list of the ranges of the selection to create the control pdfs.
    varname   : str, name of the variable of the pdf, The default is 'E'.
    varrange  : (float, float), range of the pdf variable
    bins      : int, number of bins of the histograms to create the pdfs

    Returns
    -------
    fit       : functio(DataFrame) to fit data to the Composite PDF

    """
    
    
    
    effs_signal     = np.array([ut.selection_efficiency(mc, varnames, varranges)[0] for mc in mcs])
    effs_control    = np.array([ut.selection_efficiency(mc, connames, conranges)[0] for mc in mcs])
    nevts_exp       = effs_signal * nevts
    factor_control  = effs_control/effs_signal

    # generate the ELL instace to fit the energy distribution to the energy distribution of the three mc samples
    ell_signal      = get_ell(mcs, refnames, refranges)
    ell_control     = get_ell(mcs, connames, conranges)
    ell             = efit.SimulExtComPDF(ell_signal, ell_control, factor_control)

        
    def _fit(data):
        
        data_signal  = ut.selection_sample(data, varnames, varranges)
        data_control = ut.selection_sample(data, connames, conranges)

        values       = (data_signal[varname].values, data_control[varname].values)
        result       = ell.best_estimate(values, *nevts_exp)
        
        return result, values, ell, nevts_exp

    return _fit
            


def tmu_scan(values, pars, ell, sizes = 2., nbins = 50):
    """
    
    So a -2 loglike scan in the parameters

    Parameters
    ----------
    values : data, 
    pars   : np.array, parameters
    ell    : a ComPDF object, musht have a loglike method
    sizes  : sizes of the range of the parameters to scan
    nbins  : int, number of points in the scan

    Returns
    -------
    
    nis    : list of the scan points of each parameter
    tmus   : -2loglike values of the scan

    """

    
    npars  = len(pars)
    sizes = npars * (sizes,) if isinstance(sizes, float) else sizes
    
    def n_scan(n, size = 2.):
        n0 = max(0., n - size * np.sqrt(n))
        n1 = n + size * np.sqrt(n) if n > 1 else 20
        return np.linspace(n0, n1, nbins)
    
    tmus = []
    for i in range(npars):
        nis   = n_scan(pars[i], sizes[i])
        itmus  = efit.llike_scan(values, ell, pars, nis, i)
        tmus.append((nis,itmus))
        
    return ut.list_transpose(tmus)


def tmu_values(values, par_est, ell, par_exp):
    """
    
    Compute generic tmu-values for hypothesis testing

    Parameters
    ----------
    values  : data
    par_est : np.array(float), estimated parameters
    ell     : PDF object, it must have a liglike method, loglike(data, *pars)
    par_exp : np.arrayt(float), expected parameters
        DESCRIPTION.

    Returns
    -------
    tmun  : tmu(mu, muhat)
    tmu   : tmu(mu, nuhat(x), muhat), only for the first parameter
    qmu   : qmu(mu, nuhat(x), muhat), test the alternative hypothesis
    q0    : qmu(0, nuhat(x), muhat), test the null hypothesis

    """
        
    tmun = efit.tmu(values, ell, par_est, par_exp, -1) 
    tmu  = efit.tmu(values, ell, par_est, par_exp[0], 0) 
    q0   = efit.tmu(values, ell, par_est, 0., 0)        
    qmu  = tmu if par_est[0] < par_exp[0] else 0
    
    return tmun, tmu, qmu, q0
    

#---- MC experiments


ExpResult = namedtuple('ExpResult', 
                       ('nbb', 'nBi', 'nTl', 'nbb0', 'nBi0', 'nTl0', 
                        'tmun', 'tmu', 'qmu', 'q0'))


def prepare_experiment_ell(mcs, nevts, *args, **kargs): 
    """
    Return a function that generates and analyzes data of a random mc experiment with nevts.
    Fit the data into the signal region
    
    Parameters
    ----------
    mcs   : tuple(DataFrames), list of MC DataFrame samples
    nevts : np.array(float), number of events in each sample
    *args : arguments of prepapre_fit_ell
    **kargs : key arguments of prepare_fit_simell

    Returns
    -------
    exp   : function, that generates and analyzes data of a random mc experiment

    """

    fit   = prepare_fit_ell(mcs, nevts, *args, **kargs)
    return _prepare_experiment(mcs, nevts, fit)


def prepare_experiment_simell(mcs, nevts, *args, **kargs): 
    """
    Return a function that generates and analyzes data of a random mc experiment with nevts.
    Fit data to Simulaneous signal and control samples
    
    Parameters
    ----------
    mcs   : tuple(DataFrames), list of MC DataFrame samples
    nevts : np.array(float), number of events in each sample
    *args : arguments of prepapre_fit_simell
    **kargs : key arguments of prepare_fit_simell

    Returns
    -------
    exp   : function, that generates and analyzes data of a random mc experiment

    """


    fit   = prepare_fit_simell(mcs, nevts, *args, **kargs)
    return _prepare_experiment(mcs, nevts, fit)

    
def _prepare_experiment(mcs, nevts, fit):
    
    def _experiment():
        
        mcdata = generate_mc_experiment(mcs, nevts)
        
        result, values, ell, nevts_exp = fit(mcdata)
        
        if (not result.success): 
            return result.success, mcdata, None

        nevts_est = result.x    
        tmuvals  = tmu_values(values, nevts_est, ell, nevts_exp)
        eresult  = ExpResult(*nevts_est, *nevts_exp, *tmuvals)
        
        return result.success, mcdata, eresult
    
    return _experiment

def run(experiment, size = 1):
    """
    
    run the experiment function ntimes (size), returns a dataFrame with rhe results-

    Parameters
    ----------
    experiment : function
    size       : int, number of random experiments to generate and analyze. The default is 1.

    Returns
    -------
    df         : DataFrame, with the resutls of the analysis of the size experiments

    """
    
    eresults = []
    for i in range(size):
        success, _, eresult = experiment()
        if (success): eresults.append(eresult)
            
    eresults = ut.list_transpose(eresults)
    df       = ut.list_to_df(eresults, ExpResult._fields)
    
    return df      





#--------------------------

# def ana_samples(data, 
#                 mcs, 
#                 varnames   = varnames[:-1],
#                 varranges = varranges[:-1],
#                 mc_level   = -1,
#                 verbose    = True):
#     """
    
#     Return the sample to analysis after a selection in the variables *varnames* 
#     in ranges *varranges*

#     Parameters
#     ----------
#     data      : DataFrame, data
#     mcs       : tuple(DataFrame), mc data frames for the different samples
#     varnames  : tuple(float), list of the variables names of the selection
#     varranges : tuple( (float, float)), list of the variable ranges of the selection
#     mc_level  : int, reduce the selection for mc samples to get more stats for the pdfs
#                 a given level, default = -1
#     verbose   : bool, print into

#     Returns
#     -------
#     anadata   : DataFrame, the selected data events
#     anamcs    : tuple(DataFrame), the selected mc events (to get the pdfs)
#     effs      : tuple(float), selection efficiencies in the different samples

#     """
    
#     effs    = [efficiencies(mc, varnames, varranges)[0][-1]  for mc in mcs]
#     anamcs  = [mc[selection(mc, varnames[:mc_level], varranges[:mc_level])] for mc in mcs]
#     anadata = data[selection(data, varnames, varranges)]
    
#     if (verbose):
#         print('selection variables  :', varnames)
#         print('selection ranges     :', varranges)
#         print('selection mc samples :', varnames[:mc_level])
#         print('data size            :', len(anadata))
#         print('mc sizes             :', [len(mc) for mc in anamcs])
#         print('efficiencies         : {:6.2f}, {:1.2e}, {:1.2e}'.format(* effs))
    
#     return anadata, anamcs, effs



# def ana_experiment(data,
#                    mcs,
#                    nevts,
#                    unevts    = None,
#                    varnames  = varnames[:-1],
#                    varranges = varranges[:-1],
#                    mc_level  = -1,
#                    varname   = 'E',
#                    varrange  = erange,
#                    bins      = 150,
#                    verbose   = True): 
#     """

#     perform the analysis

#     Parameters
#     ----------
#     data   : DataFrame, data
#     mcs    : tuple(DataFrame), mc data (bb, Bi, Tl)
#     nevts  : tuple(float), number of expected total events of each sample
#     unevts : tuple(float), optional, uncertainties of the expected total events in each sample
#         if provided, the fit is constrained
#     varnames  : tuple(strig), list of the variables of the selection
#     varranges : tuple( (float, float)), list with the selection variable ranges
#     mc_level  : int, optional, the mc samples for the pdfs use a 'downscaled' selection 
#     varname   : float, variable of the fit. The default is 'E'.
#     varrange  : (float, float), range of the fit. The default is erange.
#     bins      : int, optional, number of bins. The default is 150.
#     verbose   : bool, optional, print out and plots. The default is True.

#     Returns
#     -------
#     result    : object, result of the fit with the LL, and the estimated parameters
#     enes      : np.array(float), data values entered into the fit
#     ell       : object, extended LL fit object, access to the loglike
#     pdfs      : tuple(object), pdf objects associated to each sample

#     """
    
#     anadata, anamcs, effs = ana_samples(data, mcs, varnames, varranges, 
#                                         mc_level = mc_level, verbose = verbose)
#     nns  = [eff * ni  for eff, ni  in zip(effs, nevts)]
#     unns = None if unevts == None else [eff * uni for eff, uni in zip(effs, unevts)]
#     result, enes, ell  = fit_ell(anadata, anamcs, nns, unns,
#                                  varname = varname, varrange = varrange, 
#                                  bins = bins)
    
#     if (verbose):
#         print('Initial       Events : {:6.2f}, {:6.2f}, {:6.2f}'.format(* nns))
#         if (unns != None):
#             print('Uncertainties Events : {:6.2f}, {:6.2f}, {:6.2f}'.format(*unns))
#         print('Fit success          : ', result.success)
#         print('Estimated     Events : {:6.2f}, {:6.2f}, {:6.2f}'.format(*result.x))
#         plot_fit_ell(enes, result.x, ell)
    
#     return result, enes, ell, np.array(effs)
    