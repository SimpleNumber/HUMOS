#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This file contains all parameters used for simulation
'''

peptide_collection_size = 100 #number of modeled peptides
TIC = 1e9 #total ion current
low_mass = 350 #low mass cutoff
high_mass = 1500 #high mass cutoff
nScans = 2 #number of BoxCar scans
nBoxes = 12 #number of boxes per scan
box_overlap = 1 #overlap between boxes
ms2length = 43 #length of a single MS/MS scan in milliseconds (DOI: 10.1021/pr500985w)
LC_time = 60 #LC method length in minutes
resolutions_list = [7500*(2**x) for x in range(0,6)] #list of resolutions
agc_list = [1e5, 3e5, 5e5, 1e6, 3e6, 1e7] #list of AGC targets
transients = dict([(resolutions_list[x], 16*(2**x) + 6) for x in range(0,6)]) # transients from PMC4256516
ion_models = ['equal', 'lognormal', 'lognormal-major'] #list of ion models


#Tooltips for headers
dyn_range_descript = 'Shows actual dynamic range of analysing mixture (Mixture) and dynamic range of acquired spectrum (Spectral).'
pep_distr_descript = 'Generates three different abundance distribution of the peptide mixture.\nEquimolar produce a mass spectrum using equal amount of each peptide.\nRegular set all peptide abundances to imitate those experimentally observed in a proteomics experiment.\nRegular with majors is identical to regular, with the addition, that the quantity of a few peptides is greatly enhanced, i.e. in plasma samples.'
resolution_descript = 'Allows changing the mass spectral resolution'
AGC_discript = 'AGC target is the maximum total number of ions that can be collected prior to ion detection in the Orbitrap.'
IT_descript = 'The maximum injection time is the time allowed to be spent to accumulate ions and potentially reach the AGC target value.'
acquisition_discript = 'Ways to acquire ions.\n In usual MS1 all ions are accumulated simultaneously.\nIn BoxCar acquisition ions are selected in narrow “packages” according to their m/z.'
topN_discript = 'Controls the number of fragmentation events (MS2 scans) that is scheduled for each parent ion scan.'