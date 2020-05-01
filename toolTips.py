#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 16:39:46 2020

@author: julia
"""

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

#Tooltips for headers
tooltip_style = {'background-color':'white',
                 'box-shadow':  '2px 2px 5px 3px #b3b3b3',
                 'width': '350px',
                 'marginRight': '10px',
                 'marginleft': '10px',
                 'padding-left': '5px',
                 'padding-right': '5px'
                }

def logo_tooltip():
    return dbc.Tooltip(html.Img(src='/assets/trap.gif'), delay={'show':2000}, target="logo",
                        placement='left-start')

def text_tooltip(text, target):
    return dbc.Tooltip(dcc.Markdown(text, style=tooltip_style),
                       target=target, delay={'show':1000},  placement='bottom-start')

#AGC info table and Dynamic range graph
table_descript = '''Table with actual acquisition parameters.'''
dynRange_descript = '''
**Peptide** - dynamic range of peptide mixture

**MS1**, **BoxCar scans** - dynamic range of regular individual spectra

**Spectrum** - combined dynamic range of all spectra

Number to the left corresponds to covered orders of magnitude

**Visible peptides** percent of peptides in the mixture detected with current settings
'''

#Block1 distribution and Acquisition method
pep_distr_descript = '''
**Equimolar** - fictional peptide mixture having equal amount of each peptide.

**Regular** - peptide mixture similar to common proteomics samples, for
example, cell digest.
 
**Regular with majors** - peptide mixture similar to samples wherre small
number of peptides corresponding to the most of the sample, for example, plasma
samples.
'''

ionCurrent_discript = 'Ion current for acquiring spectra.'

acquisition_discript = '''
**Usual MS1** all ions are accumulated simultaneously.

**BoxCar** ions are selected in narrow __packages__ according to their _m/z_.'''

#Block2 MS1 parameters
resolution_descript = 'Allows changing the mass spectral resolution of MS1 spectra.'

AGC_discript = '''**Automatic Gate Control (AGC)** target is the maximum total number of 
               ions that can be collected prior to ion detection in the Orbitrap.'''
               
IT_descript = '''**Maximum injection time** is the time allowed to be spent to 
accumulate ions and potentially reach the corresponding AGC target value.'''


#Block3 MS2 parameters
resolutionMS2_descript = 'Allows changing the mass spectral resolution of MS2 spectra.'

topN_discript = '''Controls the number of fragmentation events (MS2 scans) that 
is scheduled for each parent ion scan.'''

parallel_descript = '''Acquire MS1 and MS2 spectra in parallel.'''