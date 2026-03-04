"""
ana.fanal_gen
=============
Data-generation module for the FANAL exercise.

Simulates a realistic NEXT-like double-beta neutrino-less decay experiment:

1. Load MC templates (bb0ν, ²¹⁴Bi, ²⁰⁸Tl) from HDF5 files.
2. Apply energy smearing (energy resolution effect).
3. Compute the expected number of signal and background events for a given
   exposure, background index, and half-life.
4. Sample a pseudo-data set and split it into blind / RoI sub-samples.
5. Optionally write the generated data sets to an HDF5 output file.

Physical defaults
-----------------
tau      : 1e26 s            – signal half-life.
sigma0   : 5.3 keV           – baseline energy resolution at Qbb.
exposure : 500 kg·y          – detector exposure.
bkgindex : 1e-4 c/(keV kg y) – background index in the RoI.
abundance: 0.9               – ¹³⁶Xe isotopic fraction.
Qbb      : 2458 keV          – ¹³⁶Xe Q-value.
W        : 135.9 g/mol       – atomic weight.
fbi      : 0.25              – fraction of Bi background vs. Tl.
"""

import numpy  as np
import pandas as pd

import scipy.constants as constants

import ana.fanal as fana


# ---------------------------------------------------------------------------
# Physical / experimental constants
# ---------------------------------------------------------------------------

erange   = (2.4,  2.7)    # MeV – full energy analysis window
eroi     = (2.43, 2.48)   # MeV – signal Region of Interest
eblob2   = 0.4            # MeV – minimum second-blob energy

tau      = 1e26            # s   – nominal bb0ν half-life
sigma0   = 5.3             # keV – nominal energy resolution at Qbb
exposure = 500             # kg·y
bkgindex = 1e-4            # counts / (keV kg y)
abundance = 0.9            # ¹³⁶Xe isotopic fraction
Qbb      = 2458            # keV
W        = 135.9           # g/mol
fbi      = 0.25            # Bi fraction of total background

samples  = ['bb0nu', 'bi214', 'tl208']
lsamples = [r'$\beta\beta0\nu$', r'$^{214}$Bi', r'$^{208}$Tl']

# Map sample name → integer label stored in the 'mc' column
_MC_LABEL = {'bb0nu': 0, 'bi214': 1, 'tl208': 2}

# Columns that must be non-NaN in the raw MC files
_REQUIRED_LABELS = ('num_tracks', 'track0_E', 'E', 'track0_length')


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_df(sample, dirpath):
    """Load one MC sample from an HDF5 file.

    Parameters
    ----------
    sample  : str – one of ``'bb0nu'``, ``'bi214'``, ``'tl208'``.
    dirpath : str – directory that contains the ``<sample>.h5`` files.

    Returns
    -------
    df  : pandas.DataFrame – cleaned MC sample.
    acc : float – acceptance = fraction of events with all required columns.
    """
    df   = pd.read_hdf(dirpath + sample + '.h5', 'events')
    df['mc'] = _MC_LABEL[sample]
    df   = df.rename(columns={'smE': 'E'})
    ntot = len(df)
    df.dropna(subset=_REQUIRED_LABELS, inplace=True)
    acc  = len(df) / float(ntot)
    return df, acc


def load_dfs(dirpath, samples=samples):
    """Load all MC samples and return them with their acceptances.

    Parameters
    ----------
    dirpath : str – directory containing the HDF5 files.
    samples : list[str] – sample names.  Default: all three samples.

    Returns
    -------
    dfs  : list[pandas.DataFrame]
    accs : list[float] – per-sample acceptance fractions.
    """
    loaded = [load_df(name, dirpath) for name in samples]
    dfs    = [x[0] for x in loaded]
    accs   = [x[1] for x in loaded]
    return dfs, accs


# ---------------------------------------------------------------------------
# Energy smearing
# ---------------------------------------------------------------------------

def energy_effect(df, sigma=sigma0, efactor=1.):
    """Apply additional energy smearing and/or a global energy scale.

    Adds Gaussian smearing with width
    ``σ_add = sqrt(σ² - σ₀²) / Qbb * E``
    (in quadrature on top of the existing baseline resolution *σ₀*), then
    multiplies all energy columns by *efactor*.

    Parameters
    ----------
    df      : pandas.DataFrame – MC events.
    sigma   : float – target resolution in keV.  If ≤ σ₀, no smearing is added.
    efactor : float – multiplicative energy scale factor.  Default 1.

    Returns
    -------
    pandas.DataFrame – copy of *df* with modified energy columns.
    """
    energy_cols = ['E', 'track0_E', 'blob1_E', 'blob2_E', 'track1_E']
    # Additional smearing width (relative to Qbb)
    extra = np.sqrt(max(0., sigma**2 - sigma0**2)) / Qbb

    df1 = df.copy()
    for col in energy_cols:
        if col not in df1.columns:
            continue
        df1[col] += np.random.normal(0., extra * df[col].values)
        df1[col] *= efactor

    return df1


# ---------------------------------------------------------------------------
# Event-count calculations
# ---------------------------------------------------------------------------

def nevents_bb0nu(exposure, tau=tau, W=W, abundance=abundance):
    """Expected number of bb0ν signal events for a given exposure and half-life.

    Parameters
    ----------
    exposure  : float – exposure in kg·y.
    tau       : float – half-life in s.  Default :data:`tau`.
    W         : float – atomic weight in g/mol.  Default :data:`W`.
    abundance : float – isotopic fraction.  Default :data:`abundance`.

    Returns
    -------
    float – expected signal count.
    """
    NA = constants.Avogadro
    return 1e3 * abundance * (exposure / tau) * (NA / W) * np.log(2.)


def nevents_bkg(exposure, roi, bkgindex, acbi, actl, fbi=fbi):
    """Expected total Bi and Tl background counts and total RoI background.

    Parameters
    ----------
    exposure : float – exposure in kg·y.
    roi      : float – RoI width in keV.
    bkgindex : float – background index in counts / (keV kg y).
    acbi     : float – total Bi acceptance × efficiency.
    actl     : float – total Tl acceptance × efficiency.
    fbi      : float – fraction of RoI background that is Bi.

    Returns
    -------
    nbi      : float – total Bi events in the exposure.
    ntl      : float – total Tl events in the exposure.
    nbkg_roi : float – total background events in the energy RoI.
    """
    nbkg_roi = bkgindex * roi * exposure
    acc_mix  = fbi * acbi + (1. - fbi) * actl
    nbi      = nbkg_roi * fbi       / acc_mix
    ntl      = nbkg_roi * (1. - fbi) / acc_mix
    return nbi, ntl, nbkg_roi


# ---------------------------------------------------------------------------
# Event selection
# ---------------------------------------------------------------------------

def selection_analysis(xdf, xroi=eroi, eblob2=eblob2):
    """Boolean selection mask for the signal region.

    Parameters
    ----------
    xdf    : pandas.DataFrame
    xroi   : tuple(float, float) – energy RoI in MeV.
    eblob2 : float               – minimum blob-2 energy in MeV.

    Returns
    -------
    numpy.ndarray of bool
    """
    return ((xdf.num_tracks == 1) &
            (xdf.E >= xroi[0]) & (xdf.E < xroi[1]) &
            (xdf.blob2_E > eblob2))


def selection_blind(df, eroi=eroi, eblob2=eblob2):
    """Boolean mask for events that should be *blinded* (removed from the sample).

    Parameters
    ----------
    df     : pandas.DataFrame
    eroi   : tuple(float, float)
    eblob2 : float

    Returns
    -------
    numpy.ndarray of bool – True for events to blind.
    """
    in_roi     = (df.track0_E >= eroi[0]) & (df.track0_E < eroi[1])
    blob2_pass = df.blob2_E > eblob2
    return in_roi | blob2_pass


# ---------------------------------------------------------------------------
# Diagnostic
# ---------------------------------------------------------------------------

def test_experiment(df, accs, exposure=exposure, eroi=eroi, eblob2=eblob2):
    """Print a quick summary of the number of events in the signal region.

    Parameters
    ----------
    df       : pandas.DataFrame – generated data.
    accs     : list[float] – total acceptance per component.
    exposure : float
    eroi     : tuple
    eblob2   : float

    Returns
    -------
    nns      : list[int]   – RoI counts per component.
    mms      : list[float] – total counts per component (extrapolated).
    tau      : float       – estimated half-life.
    bkgindex : float       – estimated background index.
    """
    sel      = selection_analysis(df, eroi, eblob2)
    nns      = [int(np.sum(df[sel].mc == i)) for i in range(3)]
    tau_est  = fana.half_life(nns[0], exposure, accs[0])
    mms      = [n / acc for n, acc in zip(nns, accs)]
    roi_keV  = 1e3 * (eroi[1] - eroi[0])
    bkg_idx  = np.sum(nns[1:]) / (roi_keV * exposure)

    print('--- Test ---')
    print('true  events in RoI    : {:3d}, {:3d}, {:3d}'.format(*nns))
    print('total events           : {:6.3f}, {:1.2e}, {:1.2e}'.format(*mms))
    print('tau             (s)    : {:1.2e}'.format(tau_est))
    print('bkg index c/(keV kg y) : {:1.2e}'.format(bkg_idx))

    return nns, mms, tau_est, bkg_idx


# ---------------------------------------------------------------------------
# Full experiment generator
# ---------------------------------------------------------------------------

def experiment(eroi     = eroi,
               eblob2   = eblob2,
               tau      = tau,
               exposure = exposure,
               bkgindex = bkgindex,
               fbi      = fbi,
               sigma    = sigma0,
               efactor  = 1.,
               dirpath  = None,
               collname = ''):
    """Generate a complete pseudo-experiment for the FANAL exercise.

    Loads MC templates, smears energy, samples events according to the
    expected counts, and splits data into blind / RoI sub-samples.

    Parameters
    ----------
    eroi     : tuple(float, float) – energy RoI in MeV.
    eblob2   : float               – minimum blob-2 energy in MeV.
    tau      : float               – bb0ν half-life in s.
    exposure : float               – exposure in kg·y.
    bkgindex : float               – background index in c/(keV kg y).
    fbi      : float               – Bi fraction of total background.
    sigma    : float               – energy resolution in keV.
    efactor  : float               – energy scale factor.
    dirpath  : str                 – path to the HDF5 MC files.
        **Required**: there is no default path; must be provided explicitly.
    collname : str                 – collaboration name tag.  If non-empty,
        the generated data are written to ``fanal_test_<collname>.h5``.

    Returns
    -------
    dfmcs   : list[pandas.DataFrame] – smeared MC samples (bb0ν, Bi, Tl).
    dfcal   : pandas.DataFrame       – Tl calibration sample (half of Tl MC).
    dfdat   : pandas.DataFrame       – full pseudo-data.
    dfblind : pandas.DataFrame       – blind sample (outside RoI).
    dfroi   : pandas.DataFrame       – RoI sample (signal region).
    """
    if dirpath is None:
        raise ValueError(
            'dirpath must be provided (e.g. dirpath="/path/to/MC/files/")')

    print(f' --- Generation : {collname} ---')
    print(f'tau       : {tau:1.2e} s')
    print(f'exposure  : {exposure:6.3f} kg y')
    print(f'bkg index : {bkgindex:1.2e}')
    print(f'roi range : ({eroi[0]:6.3f}, {eroi[1]:6.3f}) MeV')
    print(f'sigma     : {sigma:6.3f} keV')
    print(f'f Bi      : {fbi:6.3f}')
    print(f'efactor   : {efactor:6.3f}')
    print(f'eblob2    : {eblob2:6.3f} MeV')

    # ---- Load templates ----
    dfs, accs = load_dfs(dirpath)
    sels  = [selection_analysis(df, eroi, eblob2) for df in dfs]
    effs  = [float(np.sum(s) / len(df)) for s, df in zip(sels, dfs)]
    taccs = [a * e for a, e in zip(accs, effs)]  # total acceptance

    print('acceptance        : {:6.3f}, {:1.2e}, {:1.2e}'.format(*accs))
    print('efficiencies      : {:6.3f}, {:1.2e}, {:1.2e}'.format(*effs))
    print('total acceptance  : {:6.3e}, {:1.2e}, {:1.2e}'.format(*taccs))

    # ---- Expected event counts ----
    nbb     = nevents_bb0nu(exposure, tau=tau)
    nbb_roi = nbb * taccs[0]

    roi_keV               = 1e3 * (eroi[1] - eroi[0])
    nbi, ntl, nbkg_roi    = nevents_bkg(exposure, roi_keV, bkgindex,
                                         taccs[1], taccs[2], fbi)
    nbi_roi = nbi * taccs[1]
    ntl_roi = ntl * taccs[2]

    # Poisson-fluctuated total counts (after acceptance, before efficiency)
    expected_total = [nbb * accs[0], nbi * accs[1], ntl * accs[2]]
    mms = [int(np.random.poisson(n)) for n in expected_total]

    print('number of events  : {:6.2f}, {:6.2f}, {:6.2f}'.format(nbb, nbi, ntl))
    print('number in acc     : {:6.2f}, {:6.2f}, {:6.2f}'.format(*expected_total))
    print('number of bkg RoI : {:6.2f}'.format(nbkg_roi))
    print('number of RoI     : {:6.2f}, {:6.2f}, {:6.2f}'.format(nbb_roi,
                                                                   nbi_roi,
                                                                   ntl_roi))
    print('sampled int in acc: {:6d}, {:6d}, {:6d}'.format(*mms))

    # ---- Smear energy ----
    dfmcs = [energy_effect(df, sigma, 1.) for df in dfs]

    # ---- Sample pseudo-data ----
    sampled = [df.sample(n=n) for df, n in zip(dfmcs, mms)]
    xdf     = pd.concat(sampled, ignore_index=True)
    # Shuffle rows
    dfdat   = xdf.sample(frac=1.).reset_index(drop=True)

    # ---- Calibration sample (half of Tl MC) ----
    dfcal = dfmcs[2].sample(n=int(len(dfmcs[2]) / 2))

    # ---- Blind / RoI split ----
    blind_mask = selection_blind(dfdat, eroi=eroi, eblob2=eblob2)

    # Apply energy scale (no extra smearing, sigma=0 means extra=0)
    dfdat = energy_effect(dfdat, 0., efactor)
    dfcal = energy_effect(dfcal, 0., efactor)

    dfblind = dfdat[ blind_mask]
    dfroi   = dfdat[~blind_mask]

    print('total, blind, roi : {:6d}, {:6d}, {:6d}'.format(
        len(dfdat), int(np.sum(blind_mask)), int(np.sum(~blind_mask))))

    test_experiment(dfdat, taccs, exposure, eroi, eblob2)

    if collname:
        write_experiment(collname, dfmcs, dfcal, dfdat, dfblind, dfroi)

    return dfmcs, dfcal, dfdat, dfblind, dfroi


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_experiment(collname, mcs, cal, dat, blind, dfroi):
    """Write all generated DataFrames to a single HDF5 file.

    Parameters
    ----------
    collname : str – collaboration tag used in the output file name.
    mcs      : list[pandas.DataFrame] – MC samples (bb0ν, Bi, Tl).
    cal      : pandas.DataFrame       – calibration sample.
    dat      : pandas.DataFrame       – full data.
    blind    : pandas.DataFrame       – blind sample.
    dfroi    : pandas.DataFrame       – RoI sample.
    """
    ofile = 'fanal_test_' + collname + '.h5'

    mc_cols = ['mc', 'mcE', 'E',
               'num_tracks', 'num_voxels',
               'track0_E', 'track0_voxels', 'track0_length',
               'blob1_E', 'blob2_E',
               'track1_E', 'track1_voxels', 'track1_length']
    data_cols = mc_cols[2:]  # exclude 'mc' and 'mcE' for public data

    mcs[0][mc_cols[1:]].to_hdf(ofile, key='mc/bb0nu',   mode='a')
    mcs[1][mc_cols[1:]].to_hdf(ofile, key='mc/bi214',   mode='a')
    mcs[2][mc_cols[1:]].to_hdf(ofile, key='mc/tl208',   mode='a')
    dat   [mc_cols    ].to_hdf(ofile, key='mc/dummy',   mode='a')

    cal  [data_cols   ].to_hdf(ofile, key='data/tl208', mode='a')
    blind[data_cols   ].to_hdf(ofile, key='data/blind', mode='a')
    dfroi[data_cols   ].to_hdf(ofile, key='data/roi',   mode='a')
