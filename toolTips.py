#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tooltip data
"""

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

#Tooltips for headers
tooltip_style = {'background-color':'white',
                 'box-shadow':  '2px 2px 5px 3px #b3b3b3',
                 'width': '350px',
                 'margin-right': '10px',
                 'margin-left': '10px',
                 'padding-left': '5px',
                 'padding-right': '5px'
                }

def logo_tooltip():
    return dbc.Tooltip(html.Img(src='/assets/trap.gif'), delay={'show':2000}, target="logo",
                        placement='left-start')

def text_tooltip(text, target):
    return dbc.Tooltip(dcc.Markdown(text, style=tooltip_style),
                       target=target, delay={'show':1000},  placement='top-start')

#AGC info table and Dynamic range graph
info_table = 'Table with actual acquisition parameters.'

dynamic_range = '''
**Peptide** - dynamic range of peptide mixture

**MS1**, **BoxCar scans** - dynamic range of individual spectra

**Spectrum** - combined dynamic range of all spectra

Number to the left corresponds to covered orders of magnitude
'''

observed_peptides = '''
**Detected peptides** percent of peptides in the mixture detected with current
settings
'''

#Block1 distribution and Acquisition method
peptide_distribution = '''
**Equimolar** - fictional peptide mixture having equal amount of each peptide.

**Regular** - peptide mixture similar to common proteomics samples, for
example, cell digest.
 
**Regular with majors** - peptide mixture similar to samples wherre small
number of peptides corresponding to the most of the sample, for example, plasma
samples.
'''

ion_current = 'Ion current for acquiring spectra.'

acquisition = '''
**Usual MS1** all ions are accumulated simultaneously.

**BoxCar** ions are selected in narrow __packages__ according to their _m/z_.
'''

#Block2 MS1 parameters
resolution = 'Allows changing the mass spectral resolution of MS1 spectra.'

AGC = '''
**Automatic Gain Control (AGC)** target is the maximum total number of 
ions that can be collected prior to ion detection in the Orbitrap.
'''
               
MaxIT = '''
**Maximum injection time** is the time allowed to be spent to 
accumulate ions and potentially reach the corresponding AGC target value.
'''

#Block3 MS2 parameters
resolutionMS2 = '''
Allows changing the mass spectral resolution of MS2 spectra.

**IT** - MS2 scan acquired in Ion Trap.'''

topN = '''
Controls the number of fragmentation events (MS2 scans) that is scheduled for
each parent ion scan.'''

parallel = '''Acquire MS1 and MS2 spectra in parallel.'''

#Cycle time
cycle_time = '''
Visual representation of duty cycle. Each circle correspond to a specific mass
spectrometer part.

Thick parts show when a part is busy, colors correspond to spectra types.

The text in the middle is the length of the duty cycle.
'''