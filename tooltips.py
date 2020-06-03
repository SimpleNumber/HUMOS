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
info_table = 'The table with the actual acquisition parameters.'

dynamic_range = '''
**Peptide** - the dynamic range of the peptide mixture

**MS1**, **BoxCar scans** - the dynamic range of individual spectra

**Spectrum** - the combined dynamic range of all spectra

Numbers to the right correspond to the covered orders of magnitude
'''

observed_peptides = '''
**Detected peptides** - percent of peptides in the mixture that are detected
with the current settings
'''

#Block1 distribution and Acquisition method
peptide_distribution = '''
**Equimolar** - the fictional peptide mixture having equal amount of each
peptide.

**Regular** - the peptide mixture observed for common proteomics samples, for
example, a cell digest.
 
**Regular with majors** - the peptide mixture observed for samples where a
small number of peptides make the most of the sample, for example, plasma
samples.
'''

ion_current = 'The ion current used for spectra aquisition.'

acquisition = '''
**Usual MS1** - all ions are accumulated simultaneously.

**BoxCar** - ions are accumulated in narrow __packages__ according to their
 _m/z_.
'''

#Block2 MS1 parameters
resolution = 'Mass spectral resolution for MS1 spectra.'

AGC = '''
**Automatic Gain Control (AGC)** target is the maximum total number of 
ions that can be collected before the ion detection in the Orbitrap.
'''
               
MaxIT = '''
**Maximum injection time** is the time allowed to be spent to accumulate ions
and (potentially) reach the corresponding AGC target value.
'''

#Block3 MS2 parameters
resolutionMS2 = '''
Mass spectral resolution for MS2 spectra.

**IT** - MS2 spectra are acquired in an ion trap.'''

topN = '''
The number of fragmentation events (MS2 spectra) that are scheduled for
each parent ion scan (MS1 spectrum).
'''

parallel = '''Allow using different instrument parts in parallel (if possible)'''

#Cycle time
cycle_time = '''
The visual representation of the duty cycle. Each circle corresponds to a
specific mass spectrometer part.

Thick parts show when a part is busy, colors correspond to the spectra types.

The text in the middle is the length of the duty cycle.
'''