#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 11:17:59 2022

@author: hernando
"""

from collections import namedtuple

Measurement = namedtuple('Measurement', ('value', 'uncertainty'))

# energies
sigma = 8    # keV, energy resolution at Qbb
Qbb   = 2458 # keV, energy at Qbb value 
EBi   = 2448 # keV, energy of the Bi photo-peak
ETl   = 2614 # keV, energy of the Tl phoot-peak

# exposure
exposure = 500 # kg y

# Bkg Events from fit to the blind sample
nevts_Bi_blind   = Measurement(117.13, 14.6)
nevts_Tl_blind   = Measurement(791.87, 30.1)

# efficiencies
acc_bb        = 0.794
acc_Bi        = 2.36e-01
acc_Tl        = 6.33e-01

eff_Bi_blind  = 0.793
eff_Tl_blind  = 0.709

eff_bb_RoI    = 0.51885
eff_Bi_RoI    = 1.55e-2 
eff_Tl_RoI    = 1.94e-4

# Estimated Bkg Event (total and in RoI)
nevts_Bi     = Measurement( 147.8, 28.9)
nevts_Tl     = Measurement(1116.3, 21.1)

nevts_Bi_RoI = Measurement(2.29, 0.232)
nevts_Tl_RoI = Measurement(0.22, 0.003)

# bkg index
bkgindex     = 1e-4 # counts/ (keV kg y) 
