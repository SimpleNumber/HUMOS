#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 13:06:14 2020

@author: julia
"""
import numpy as np
import plotly.graph_objs as go
import math

#layout for main plot
def get_main_layout(x_range, y_range):
    return go.Layout(showlegend=True,
                     margin={'t':30},
                     xaxis={'title': 'm/z', 'range':x_range},
                     yaxis={'title': 'Intensity', 'range': y_range})

#dynamic range trace
def get_dr_tr(row):
    trace = go.Scatter(x=row['x'],
                       y=row['y'],
                       line={'width': 7, 'color':row['color']},
                       mode='lines+text',
                       text=['{:.2f}'.format(np.log10(row['x'][0] / row['x'][1])), row['text']],        
                       textposition=['middle right', 'middle left']) 
    return trace

#layout for dynamic range plot    
def get_dr_layout(dr_df):
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

#observed peptides trace
def get_obsPep_tr(observed_peptides, observed_color, missing_color):
    traces = [go.Bar(x=[0],
                     y=[observed_peptides],
                     width=1,
                     orientation='v',
                     name='% observed peptides',
                     text=str(observed_peptides),
                     textposition='inside',
                     marker_color=observed_color
                    ),
              go.Bar(x=[0],
                     y=[100 - observed_peptides],
                     width=1,
                     orientation='v',
                     name='% missing peptides',
                     marker_color=missing_color
                    )
              ]   
    return traces

#observed peptides layout
def get_obsPep_layout():
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

def get_theta_ranges(theta_range):
    theta = np.arange(theta_range[0], theta_range[1], -0.2)
    return np.append(theta, theta_range[1])

    
def get_cycle_grid_tr():
    grid_trace = go.Scatterpolar(r=[0.5] * 50 + [0.7] * 50 + [0.9] * 50,
                    theta=list(np.linspace(90, 450, 50)) * 3,
                    mode="lines",
                    line={'width': 1, 'color':'#cccccc'},
                    showlegend=False,
                    hoverinfo='skip',)
    return grid_trace

def get_cycle_annotations_tr(cycletime):
    annotation_trace = go.Scatterpolar(r=[0.0, 0.57, 0.77, 0.97],
                                       theta=[90] * 4,
                                       mode="text",
                                       text=['{:.3f} sec'.format(cycletime/1000),
                                             'Ion Trap', 'Orbitrap', 'Ion Accumulation'],
                                       showlegend=False,
                                       textfont={"size": [18] + [11] * 3},
                                       textposition='middle center',
                                       hoverinfo='text')
    return annotation_trace

def get_cycle_text_tr(ms1_scan_text, ms2_scan_text):
    # r1 = 1.1
    # theta1 = -15
    # theta2 = -25
    # r2 = r1 * math.cos(abs(theta1) * np.pi / 180) / math.cos(abs(theta2) * np.pi / 180)
    text_trace = go.Scatterpolar(r=[1.05, 1.2],
                                 theta=[-90, -90],
                                 mode='text',
                                 text= [ms1_scan_text, ms2_scan_text],
                                 showlegend=False,
                                 textfont={"size": [15] * 2},
                                 textposition='bottom center',
                                 hoverinfo='none')
    return text_trace


    
def get_cycle_tr(row):
    trace = go.Scatterpolar(r=[row['r']] * len(row['theta']),
                            theta=row['theta'] ,
                            mode=row['mode'],
                            text=row['text'],
                            name=row['name'],
                            showlegend=row['showlegend'],
                            line={'width': row['line_width'], 'color':row['line_color']},
                            textposition='bottom right',
                            hoverinfo='text')
    return trace



def get_cycle_layout():
    return go.Layout(polar = {'radialaxis': {'visible':False},
                              'angularaxis': {'visible':False}
                             },
                     showlegend=True,
                     legend={'x': 0.95,
                             'y': 0.85},
                     margin={'l':0,
                             'r':300,
                             'b':0,
                             't':0,
                             'pad':4}
                     ) 