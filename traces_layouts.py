#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 13:06:14 2020

@author: julia
"""
import numpy as np
import plotly.graph_objs as go
import math

def get_main_layout(x_range, y_range):
    '''
    Main plot layout

    Parameters
    ----------
    x_range : tuple
        range of x-coordinates - [min, max].
    y_range : tuple
        range of y-coordinates.
    '''
    return go.Layout(showlegend=True,
                     margin={'t':30},
                     xaxis={'title': 'm/z', 'range':x_range},
                     yaxis={'title': 'Intensity', 'range': y_range})

def get_dr_tr(row):
    '''
    Single trace in dynamic range plot
    '''
    return go.Scatter(x=row['x'],
                      y=row['y'],
                      line={'width': 7, 'color':row['color']},
                      mode='lines+text',
                      text=['{:.2f}'.format(np.log10(row['x'][0] / row['x'][1])), row['text']],        
                      textposition=['middle right', 'middle left']) 

def get_dr_layout(dr_df):
    '''
    Dynamic range plot layout

    Parameters
    ----------
    dr_df : pandas.DataFrame
        all information for Dynamic Range plot.
    '''
    return go.Layout(margin={'t': 0,
                             'l': 10,
                             'b': 40},
                     xaxis={'title': 'Abundance',
                            'type': 'log',
                            'exponentformat': 'power',
                            'range': [np.log10(dr_df.at['Peptide','x'][1]) - 2, 
                                      np.log10(dr_df.at['Peptide','x'][0]) + 1] },
                     yaxis={'visible': False,
                            'range': [-1, len(dr_df)]},
                     showlegend=False,
                     width= 400,
                     height=180,
                     hovermode=False)

def get_obsPep_tr(observed_peptides, observed_color, missing_color):
    '''
    Traces for observed peptides plot

    Parameters
    ----------
    observed_peptides : float
        percentage of peptides observed, should be between 0 and 100.
    observed_color : string
        color code for observed peptides bar.
    missing_color : string
        color code for missing peptides bar.
    '''
    return [go.Bar(x=[0],
                   y=[observed_peptides],
                   width=1,
                   orientation='v',
                   name='% observed peptides',
                   text=str(observed_peptides),
                   textposition='inside',
                   marker_color=observed_color),
            
            go.Bar(x=[0],
                   y=[100 - observed_peptides],
                   width=1,
                   orientation='v',
                   name='% missing peptides',
                   marker_color=missing_color)
            ]

def get_obsPep_layout():
    '''
    Observed peptides plot layout
    '''
    return go.Layout(margin={'t':10,
                             'l':40,
                             'r':10,
                             'b': 10},
                     xaxis={'visible': False},
                     yaxis={'title': '% visible peptides',
                            'range': [0, 100]},
                     showlegend=False,
                     barmode='stack',
                     width=100,
                     height=180,
                     hovermode=False)

def get_theta_range(theta_range, step=1):
    '''
    Generate array with points in a range [start, stop] with step distance
    between points. Both ends are included.

    Parameters
    ----------
    theta_range : tuple
        range defined as [start, stop].
    step : float
        distance between points

    Returns
    -------
    np.ndarray
        range points.

    '''
    step = -1 * step if theta_range[1] < theta_range[0] else step
    theta = np.arange(theta_range[0], theta_range[1], step)
    return np.append(theta, theta_range[1])
    
def get_cycle_grid_tr():
    '''
    Circular grid of cycle time plot
    '''
    return go.Scatterpolar(r=[0.5] * 120 + [0.7] * 120 + [0.9] * 120,
                           theta=np.concatenate([np.linspace(90, 450, 120)]*3),
                           mode='lines',
                           line={'width': 1, 'color':'#cccccc'},
                           showlegend=False,
                           hoverinfo='skip')

def get_cycle_annotations_tr(cycletime):
    '''
    Annotations of circular grid
    
    Parameters
    ----------
    cycletime : float
        length of the duty cycle.
    '''
    return go.Scatterpolar(r=[0.0, 0.57, 0.77, 0.97],
                           theta=[90] * 4,
                           mode='text',
                           text=['{:.3f} sec'.format(cycletime/1000),
                                 'Ion Trap', 'Orbitrap', 'Ion Accumulation'],
                           showlegend=False,
                           textfont={'size': [18] + [11] * 3},
                           textposition='middle center',
                           hoverinfo='skip')

def get_cycle_text_tr(ms1_scan_text, ms2_scan_text):
    '''
    Number of MS1 and MS2 scans, located at the cycle plot

    Parameters
    ----------
    ms1_scan_text : string
        Text with MS1 scan information.
    ms2_scan_text : string
        Text with MS2 scan information.
    '''
    r1 = 1.1
    theta1 = -115
    theta2 = -125
    r2 = r1 * math.cos(abs(theta1) * np.pi / 180) / math.cos(abs(theta2) * np.pi / 180)
    text_trace = go.Scatterpolar(r=[r1, r2],
                                 theta=[theta1, theta2],
                                 mode='text',
                                 text= [ms1_scan_text, ms2_scan_text],
                                 showlegend=False,
                                 textfont={'size': [14] * 2},
                                 textposition='bottom right',
                                 hoverinfo='skip')
    return text_trace
    
def get_cycle_tr(row):
    '''
    Single trace in a circular plot
    '''
    return go.Scatterpolar(r=[row['r']] * len(row['theta']),
                           theta=row['theta'] ,
                           mode=row['mode'],
                           text=row['text'],
                           name=row['name'],
                           showlegend=row['showlegend'],
                           line={'width': row['line_width'], 'color':row['line_color']},
                           textposition='bottom right',
                           hoverinfo='text')



def get_cycle_layout():
    '''
    Cycle Time plot layout
    '''
    return go.Layout(polar={'radialaxis': {'visible':False},
                            'angularaxis': {'visible':False} },
                     showlegend=True,
                     legend={'x': 0.95,
                             'y': 0.85},
                     margin={'l':0,
                             'r':250,
                             'b':0,
                             't':0,
                             'pad':4})
