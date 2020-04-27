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
                 'width': '400px'}
def logo_tooltip():
    return dbc.Tooltip( html.Img (src='/assets/trap.gif'), delay={'show':2000}, target="logo")

def text_tooltip(text, target):
    return dbc.Tooltip( dcc.Markdown(text, style=tooltip_style),
                       target=target)
 #AGC info table and Dynamic range graph
table_descript = '''Table with actual acquisition parameters.'''

dynRange_descript = '''Shows actual dynamic range of analysing mixture (**Mixture**) 
                    and dynamic range of acquired spectrum (**Spectral**).'''
dynRange_bar1_descript = '''Shows actual dynamic range of analysing mixture (**Mixture**) 
                         and dynamic range of acquired spectrum (**Spectral**).'''
observed_peptides_descript = '''Shows persent of peptides in the mixture visiable in the experiment.'''    


#Block1 distribution and Acquisition method
pep_distr_descript = '''Generates three different abundance distribution of the peptide mixture.'''
pep_distr_rb_descript = '''**Equimolar** produce a mass spectrum using equal amount of each peptide.
                        **Regular** set all peptide abundances to imitate those experimentally observed in a proteomics experiment.
                        **Regular with majors** is identical to regular, with the addition, that the quantity of a few peptides
                        is greatly enhanced, i.e. in plasma samples.'''

ionCurrent_discript = 'Ion current for acquiring spectra.'

acquisition_discript = 'Ways to acquire ions.'
acquisition_rb_discript = '''In **usual MS1** all ions are accumulated simultaneously.
                          In **BoxCar** acquisition ions are selected in narrow “packages” according to their _m/z_.'''


 #Block2 MS1 parameters
resolution_descript = 'Allows changing the mass spectral resolution of MS1 spectra.'
AGC_discript = '''**Automatic Gate Control (AGC)** target is the maximum total number of 
               ions that can be collected prior to ion detection in the Orbitrap.'''
IT_descript = '''The **maximum injection time** is the time allowed to be spent to accumulate ions and potentially 
              reach the AGC target value.'''
IT_button_description = 'Set maximum injection time for MS1 spectra.'

#Block3 MS2 parameters
resolutionMS2_descript = 'Allows changing the mass spectral resolution of MS2 spectra.'
IT_MS2_descript = '''The **maximum injection time** is the time allowed to be spent to 
                  accumulate ions and potentially reach the AGC target value for MS2 spectra.'''
IT_MS2_button_descript = 'Set maximum injection time for MS2 spectra.'               
topN_discript = 'Controls the number of fragmentation events (MS2 scans) that is scheduled for each parent ion scan.'
IT_button_descript = 'Set maximum injection time for MS2 spectra.'
parallel_descript = '''Acquire MS1 and MS2 in parallel.'''