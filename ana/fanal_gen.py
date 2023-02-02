import numpy  as np
import pandas as pd

import scipy.constants as constants

import ana.fanal as fana

#import scipy.stats     as stats

#import core.utils as ut
#import core.efit  as efit


#import core.pltext as pltext
#import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable


label   = 'track0_E'
erange  = (2.4, 2.7)
eroi    = (2.43,  2.48)
eblob2  = 0.4
#bins    = 100

tau      = 1e26 # s
sigma0   = 5.3  # keV
exposure = 500  # kg year
bkgindex = 1e-4 # counts / (keV kg y )
abundace = 0.9
Qbb      = 2458 # keV Qbb value
W        = 135.9
fbi      = 0.25


# Load data
#--------

dirpath  = '/Users/hernando/work/investigacion/NEXT/data/MC/FANAL/'
samples  = ['bb0nu', 'bi214', 'tl208']
lsamples = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']
mcs_numbers = {'bb0nu' : 0, 'bi214' : 1, 'tl208': 2}


def load_df(sample, dirpath = dirpath):
    """ Returns the original DF of a given sample
    """
    labels = ('num_tracks', 'track0_E', 'E', 'track0_length')
    df = pd.read_hdf(dirpath + sample  + '.h5', 'events')
    df['mc'] = mcs_numbers[sample]
    df = df.rename(columns = {"smE": "E"})
    ntot = len(df)
    df.dropna(subset = labels, inplace=True)
    acc = len(df)/float(ntot)
    return df, acc

def load_dfs(dirpath = dirpath, samples = samples):
    """
    Returns the original MC DF

    Parameters
    ----------
    dirpath : (string), path to the input data
    samples : (typle of string), names of the samples: bb0nu, Bi, Tl

    Returns
    -------
    dfs    : (tuple of DF) the DF for the different samples
    accs   : (tuple of floats) acceptance (events non-null in the DF respect the total)

    """
    _dfs    = [load_df(name, dirpath) for name in samples] 
    dfs  = [x[0] for x in _dfs] 
    accs = [x[1] for x in _dfs]
    return dfs, accs
    
#   Energy operations
#----------------------------


def energy_effect(df, sigma = sigma0, efactor = 1.):    
    
    labels = ['E', 'track0_E', 'blob1_E', 'blob2_E', 'track1_E']
    factor = np.sqrt(sigma**2 - sigma0**2)/Qbb if sigma > sigma0 else 0.
    df1 = df.copy()
    for label in labels:
        df1[label] += np.random.normal(0., factor * df[label].values)
        df1[label] *= efactor
    return df1



#  Number of Event operations
#----------------------------------

def nevents_bb0nu(exposure, tau = tau, W = W, abundance = abundace):
    """ number of bb0nu
    inputs:
        Mt  : float, exposure in kg/y
        tau : float, bb0nu half-life y
        W   : float, atomic number gr/mol
        epsilon: fliat, efficiency factor
    return:
        number of bb0nu events (float)
    """
    NA  = constants.Avogadro
    nbb = 1e3 * abundance * (exposure / tau) * (NA / W) * np.log(2.)
    return nbb


def nevents_bkg(exposure, roi, bkgindex, acbi, actl, fbi = fbi):
    """
    Returns the number total Bi and Tl, and the number of events in RoI
    
    Parameters
    ----------
    exposure : float, exposure, kg year
    roi      : float, size of energy RoI (keV)
    bkgindex : float, counts/ (keV kg yeat)
    acbi     : float, acceptance fraction (total acceptance of Bi events in RoI)
    actl     : float, acceptance fraciton (total acceptance of Tl events in RoI)
    fbi      : float, fraction, fraction of Bi bkg events in RoI of the total

    Returns
    -------
    nbi      : float, total number of Bi events in exposure
    ntl      : float, total number of Tl events in exposure
    nbkg_roi : float, number of bkg events in Energy RoI
    """
    
    
    nbkg_roi = bkgindex * roi * exposure
    acc      = (fbi * acbi + (1-fbi) * actl)
    nbi      = nbkg_roi * fbi /acc
    ntl      = nbkg_roi * (1-fbi) /acc
    return nbi, ntl, nbkg_roi


def selection_analysis(xdf,  
                       xroi   = eroi, 
                       eblob2 = eblob2):
    sel = (xdf.num_tracks == 1) & (xdf.E >= xroi[0]) & \
        (xdf.E < xroi[1]) & (xdf.blob2_E > eblob2)
    return sel


def selection_blind(df, eroi = eroi, eblob2 = eblob2):
    sel0 = (df.track0_E >= eroi[0]) &  (df.track0_E < eroi[1])
    sel1 = (df.blob2_E  > eblob2)
    sel  = np.logical_or(sel0, sel1)
    return sel


def test_experiment(df,
                    accs, 
                    exposure = exposure,
                    eroi     = eroi,
                    eblob2   = eblob2):
    
    sel      = selection_analysis(df, eroi, eblob2)
    nns      = [np.sum(df[sel].mc == i) for i in range(3)]
    tau      = fana.half_life(nns[0], exposure, accs[0])
    mms      = [n / acc for n, acc in zip(nns, accs)]
    roi      = 1e3*(eroi[1] - eroi[0])
    bkgindex = np.sum(nns[1:])/(roi*exposure)
    print(' --- Test --- ' )
    print('true  events in RoI    : {:3d}, {:3d}, {:3d}'.format(*nns))
    print('total events           : {:6.3f}, {:1.2e}, {:1.2e}'.format(*mms))
    print('tau             (s)    : {:1.2e}'.format(tau))
    print('bkg index c/(keV kg y) : {:1.2e}'.format(bkgindex))
    return nns, mms,tau, bkgindex



def experiment(eroi     = eroi,  # energy range
               eblob2   = eblob2,
               tau      = tau,   # s
               exposure = exposure,    # kg y
               bkgindex = bkgindex, # counts / (kg y keV)
               fbi      = fbi,
               sigma    = sigma0,
               efactor  = 1.,
               collname = ''):
    """
    
    Generate the data for the open-exercise of fanal: a double-beta 
    neutrino-less search

    Parameters
    ----------
    eroi     : (float, float), range of the Energy RoI in MeV, default: eroi
    eblob2   : (float), cut on the energy of the 2nd blob in MeV, default: eblob2
    tau      : (float, half-life time in s, detault: tau
    exposure : (float), exposure in kg y, default: exposure
    bkgindex : (float), bkg index in RoI, counts/ (keV kg y), default: exposure
    fbi      : (float), fraction of Bi respect Tl, default: fbi
    sigma    : (float), energy resolution in KeV, default: sigma0
    efactor  : (float), energy scale factor, default: efactor
    collname : (string), name of the collaboration, if given a name, the 
        dataframe are written into an output file, default: ''.
    Returns
    -------
    dfmcs    : (tuple) of Data Frames for the MC: bb0mu, Bi, ÇTl
    dfcal    : DF, with the Tl calibration data
    dfdat    : DF, data
    dfblind  : DF, data to do the blind analysis (with the signal region removed)
    dfroi    : DF, data in the signal region
    """
    

    print(' --- Generation : {:s} ---'.format(collname))
    print('tau       : {:1.2e} s'.format(tau))
    print('exposure  : {:6.3f} kg y'.format(exposure))
    print('bkg index : {:1.2e}'.format(bkgindex))
    print('roi range : ({:6.3f}, {:6.3f}) MeV'. format(*eroi))
    print('sigma     : {:6.3f} keV'.format(sigma))
    print('f Bi      : {:6.3f}'.format(fbi))
    print('efactor   : {:6.3f}'.format(efactor))
    print('eblob2    : {:6.3f} MeV'.format(eblob2))
    factor = np.sqrt(sigma**2 - sigma0**2)/Qbb if sigma > sigma0 else 0.
    print('efactor   : {:1.3e}'.format(factor))
    print('  --- Events --- ')
    # load data
    dfs, accs = load_dfs()
    sels      = [selection_analysis(df, eroi, eblob2) for df in dfs]
    effs      = [float(np.sum(sel)/len(df)) for sel, df in zip(sels, dfs)]
    
    # total acceptance
    taccs     = [acc * eff for acc, eff in zip(accs, effs)]
    print('acceptance        : {:6.3f}, {:1.2e}, {:1.2e}'.format(*accs))
    print('efficiencies      : {:6.3f}, {:1.2e}, {:1.2e}'.format(*effs))
    print('total acceptance  : {:6.3e}, {:1.2e}, {:1.2e}'.format(*taccs))
    
    # number of signal events
    nbb       = nevents_bb0nu(exposure, tau = tau)
    nbb_roi   = nbb * taccs[0] # number of signal events in RoI
    
    # number of bkg events
    roi                = 1e3 * (eroi[1] - eroi[0]) # size of RoI in keV
    print('Roi (keV)         : {:6.3f} keV'.format(roi))
    accbi, acctl       = taccs[1], taccs[2] # total Bi and Tl acceptances
    nbi, ntl, nbkg_roi = nevents_bkg(exposure, roi, bkgindex, accbi, acctl, fbi)
    nbi_roi, ntl_roi   = nbi * taccs[1], ntl * taccs[2]
    
    # random int number of total events - Poisson
    mms_  = [n * acc for n, acc in zip((nbb, nbi, ntl), accs)]
    mms   = [np.random.poisson(n, 1)[0] for n in mms_]

    print('number of events  : {:6.2f}, {:6.2f}, {:6.2f} '.format(nbb, nbi, ntl))
    print('number in acc     : {:6.2f}, {:6.2f}, {:6.2f} '.format(*mms_))
    print('number of bkr Roi : {:6.2f}'.format(nbkg_roi))
    print('number of Roi     : {:6.2f}, {:6.2f}, {:6.2f} '.format(nbb_roi, 
                                                                  nbi_roi, 
                                                                  ntl_roi))
    print('number int in acc : {:6.2f}, {:6.2f}, {:6.2f} '.format(*mms))

    # smeare the Energy of MC
    dfmcs = dfs
    dfmcs = [energy_effect(df, sigma, 1.) for df in dfs]
    
    # generate the data
    xdfs  = [df.sample(n = n) for df, n in zip(dfmcs, mms)]
    xdf   = pd.concat((xdfs), ignore_index = True)
    dfdat = xdf.iloc[np.random.permutation(xdf.index)].reset_index(drop = True)
    
    dfcal  = dfmcs[2].sample(n = int(len(dfmcs[2])/2))
    
    
    sel     = selection_blind(dfdat, eroi = eroi, eblob2 = eblob2)

    # escale the energy (conversion to pes)    
    dfdat   = energy_effect(dfdat, 0., efactor)
    dfcal   = energy_effect(dfcal, 0., efactor)
    dfblind = dfdat [~sel]
    dfroi   = dfdat [sel]
    
    print('total, roi, blind : {:6.2f}, {:6.2f}, {:6.2f} '.format(len(dfdat),
                                                                  np.sum(sel),
                                                                  np.sum(~sel)))
    
    test_experiment(dfdat, taccs, exposure, eroi, eblob2)
    
    if (collname != ''):
        write_experiment(collname, dfmcs, dfcal, dfdat, dfblind, dfroi)
    
    return dfmcs, dfcal, dfdat, dfblind, dfroi



def write_experiment(collname, 
                     mcs,
                     cal,
                     dat,
                     blind,
                     dfroi):
    
    ofile = 'fanal_test_' + collname + '.h5'

    labels = ['mc', 'mcE', 'E',
              'num_tracks', 'num_voxels', 
              'track0_E', 'track0_voxels', 'track0_length',
              'blob1_E', 'blob2_E', 
              'track1_E', 'track1_voxels', 'track1_length']
    
    mcs[0][labels[1:]].to_hdf(ofile, key = 'mc/bb0nu', mode = 'a')
    mcs[1][labels[1:]].to_hdf(ofile, key = 'mc/bi214', mode = 'a')
    mcs[2][labels[1:]].to_hdf(ofile, key = 'mc/tl208', mode = 'a')
    dat   [labels]    .to_hdf(ofile, key = 'mc/dummy', mode = 'a')

    cal[labels[2:]]   .to_hdf(ofile, key = 'data/tl208', mode = 'a')

    blind [labels[2:]].to_hdf(ofile, key = 'data/blind', mode = 'a')
    dfroi [labels[2:]].to_hdf(ofile, key = 'data/roi'  , mode = 'a')
    return


#----------------

# ssamples = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']

# erange     = (2.400, 2.650)
# eroi       = (2.440, 2.475)
# keys       = ['E', 'num_tracks', 'blob2_E', 'RoI']
# varnames   = ['E', 'num_tracks', 'blob2_E', 'E']
# varranges  = [erange, (1., 1.1), (0.4, np.inf), eroi]

# blindvar   = 'track0_E'
# blindrange = (2.420, 2.520)


# def half_life(nbb, exposure, eff, acc = 0.8, W = 136):
#     """  Compute the half-life time
#     inputs:
#         nbb     : float, number of events in RoI
#         exposure: float, (kg y)
#         eff     : signal efficiency, (in fraction)
#         acc     : float, isotope fraction (0.9)
#         W       : float, Atomic weight, (136 g/mol for 136Xe)
#     returns:
#         tau     : float, half-life (y)
#     """
#     NA   = constants.Avogadro
#     tau  = 1e3 * eff * acc * (exposure / nbb) * (NA / W) * np.log(2.)
#     return tau


# def efficiencies(df, names = varnames, ranges = varranges):
#     """ returns the efficiencies and its uncertatines for a serie of selections
#     inputs:
#         df    : DF
#         names : tuple(str), names of the variables to select a given range
#         ranges: tuple( (float, float), ), ranges of the variables
#     returns:
#         eff   : tuple(float), efficiency of each cut in the serie
#         ueff  : tuple(float), uncertainty of the efficiency in each cut of the serie
#     """

#     sels = ut.selections(df, names, ranges)
#     effs = [ut.efficiency(sel) for sel in sels]
#     eff  = [x[0] for x in effs]
#     ueff = [x[1] for x in effs]
#     return eff, ueff


# def blind_mc_samples(mcs, blindvar = blindvar, blindrange = blindrange):
#     """ return MC sample with a blind region in the variable *blindvar* and range *blindrange*
#     inputs:
#         mcs        : tuple(DF), DFs with the MC samples
#         blindvar   : str, name of the variable used to blind
#         blindrange : tuple(float, float), range of the variable to blind
#     returns:
#         mcs        : tuple(DF), blind DFs
#     """
#     bmcs = []
#     for mc in mcs:
#         sel = ~ut.selection(mc, blindvar, blindrange)
#         bmcs.append(mc[sel])
#     return bmcs


# def generate_mc_experiment(mcs, nevents):
#     """ generate a MC experiment with the mcs samples, mcs, and the events in each sample, nevents
#     inputs  :
#         mcs     : tuple(DF), DFs with the mc samples
#         nevents : tuple(int), number of events to select in each sample
#     returns :
#         mc      : DF, DF with the mix of number event of events, nevents, of the initial samples, mcs
#     """
#     xmcs = [mc.sample(n = int(ni)) for mc, ni in zip(mcs, nevents)] # sample ni events of mc-sample mci
#     mc    = pd.concat(xmcs) # concatenate all background to generate a mc sample that mimics the data-blind sample
#     return mc


# def fit_ell(data, mcs, ns, varname = 'E', bins = 100, varrange = erange):
#     """ Fit the variable 'varname' ('E') of the DF of the data (data) to a combined PDF
#     of the MC samples (mcs), using a histogram with number of bins (bins) in the range of
#     the variable (varrange)
#     inputs:
#         data    : DF, DF of the data
#         mcs     : tuple(DF), list of mcs samples DFs
#         varname : str, name of the variable to fit, to generate the combined PDF
#         varrange: (float, float), range of the histogram of the PDFs
#     returns:
#         results : A FitResult object with the fit results (see optimize.minimize in Scipy)
#         x       : the data used in the fit,
#                   that is the values of the data varialbe (varname) in range (varrange)
#         ell     : ExtComPDF object, An object with the functionality of a Extended Combined PDF
#                   (see efit module)
#     """

#     # generate pdfs from the MC
#     pdfs = [stats.rv_histogram(np.histogram(mc[varname], bins, range = varrange))
#             for mc in mcs]

#     # generate ELL object with the PDFs
#     ell  = efit.ExtComPDF(pdfs, *ns)

#     # get the data variable
#     sel = ut.selection(data, varname, varrange)
#     x = data[varname][sel]

#     # fit
#     result = ell.best_estimate(x, *ns)
#     nsbest = result.x

#     return result, x, ell


# def plot_fit_ell(x, par, pdf, bins = 100, parnames = ssamples):
#     """ plot the data x, and superimposed the pdf with parameters (par)
#     inputs:
#         x   : np.array(float), data to plot
#         par : tuple(float), parameters of the pdf
#         pdf : function, pdf(x, *par), returns the pdf values of the distribution along x
#         parnames: tuple(str), list of the parameters names to write them in the plot legend
#     """

#     subplot = pltext.canvas(1, 1, 8, 10)
#     subplot(1)

#     counts, edges = np.histogram(x, bins);
#     centers = 0.5 * (edges[1:] + edges[:-1])
#     ecounts = np.sqrt(counts)
#     sel     = ecounts > 0
#     nn     = np.sum(par)
#     factor = nn * (centers[1] - centers[0])
#     plt.errorbar(centers[sel], counts[sel], yerr = ecounts[sel], marker = 'o', ls = '', label = 'data')

#     label  = 'ELL fit \n'
#     for si, ni in zip(parnames, par):
#         label += ' {:s} : {:6.2f} \n'.format(si, ni)
#     plt.plot(centers, factor * pdf(centers, *par), label = label)
#     plt.legend(); plt.grid();

#     ax = plt.gca()
#     divider = make_axes_locatable(ax)
#     ax2 = divider.append_axes("bottom", size = '20%', pad = 0)
#     ax.figure.add_axes(ax2)
#     fun = lambda x, *p : factor * pdf(x, *p)
#     pltext.hresiduals(x, bins, fun, par)

#     return



# def ana_experiment(data, mcs, nevts, level_data = 2):

#     #level_mc = 1 if selection is None else level_mc
#     #print(level_mc)
#     #vnames, vranges = (varnames[:-1], varranges[:-1]) if selection is None else selection

#     vnames   = ['E', 'num_tracks', 'blob2_E']
#     vranges  = [(2.400, 2.650), (1., 1.1), (0.4, np.inf)]
#     level_mc = min(1, level_data) # MC level selction E in range, 1 track

#     sels = [ut.selections(mc, vnames, vranges) for mc  in mcs]
#     effs = [ut.efficiency(sel[level_data])     for sel in sels]
#     tmcs = [mc[sel[level_mc]]                  for mc, sel in zip(mcs, sels)]
#     tdat = data[ut.selections(data, vnames, vranges)[level_data]]

#     # expected number of events after selection
#     ns_exp       = [ni * eff[0] for ni, eff in zip(nevts, effs)]
#     res, x, ell  = fit_ell(tdat, tmcs, ns_exp)

#     return (res, x, ell, effs)






# #----------------
# #
# # def test_pdf(dat, mcs, varname, bins, varrange, ssamples = ssamples):
# #     ene  = dat[varname].values
# #     ene  = ene[ut.in_range(ene, erange)]
# #     print('events', len(ene))
# #
# #     subplot = pltext.canvas(2)
# #
# #     subplot(1)
# #     pltext.hist(ene, bins, label = 'blind data');
# #     plt.xlabel(varname);
# #
# #     subplot(2)
# #     for i, mc in enumerate(mcs):
# #         pltext.hist(mc[varname], bins, range = erange, density = True, label = ssamples[i])
# #     plt.xlabel(varname);
