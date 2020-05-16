#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 16 13:06:14 2020

@author: julia
"""
import numpy as np
import pandas as pd
import mechanics, params, toolTips 
import plotly.graph_objs as go

def get_theta_ranges(theta_range):
    theta = np.arange(theta_range[0], theta_range[1], -0.2)
    return np.append(theta, theta_range[1])

def lightening_color(rgb_color, coef=0.5):
    r, g, b = [int(i) for i in rgb_color[4:-1].split(',')]
    return 'rgb({},{},{})'.format(min(int(r + r * coef), 255),
                                  min(int(g + g * coef), 255),
                                  min(int(b + b * coef), 255))
    

def get_cycle_layout():
    rs = [0.5, 0.7, 0.9]
    return go.Layout(polar = {'radialaxis': {'angle':90, 
                                                     'range':[0, 1.2], 
                                                     'showticklabels':False, # 'tickvals':rs, 'textvals':rs,    
                                                     'tickvals':rs,#'tick':['IT', 'OT', 'Ion Accumulation'],
                                                     'visible':False, 
                                                     'color':'#cccccc',
                                                     'showline':True},
                                      'angularaxis': {'showticklabels': False, 
                                                    'ticks':'',
                                                    'visible':False }},
                            showlegend=True,
                            legend= {'x': 0.95,
                                   'y': 0.85},
                            margin={'l':0,
                                  'r':300,
                                  'b':0,
                                  't':0,
                                  'pad':4}
                            ) 

def get_cycle_annotations_tr(cycletime):
    annotation_trace = go.Scatterpolar(r=[0.0,  0.55, 0.75, 0.95],
                    theta=[90] + [90]*5,#'Acquisition in Ion Trap',
                    mode="text",
                    text=['{} sec'.format(cycletime/1000),
                          'Ion Trap','Orbitrap', 'Ion Accumulation'],
                    showlegend=False,
                    textfont={"size": [18]*1+[11]*3},
                    textposition='middle center',
                    hoverinfo='text'
                    )
    return annotation_trace

def get_cycle_text_tr(ms1_scan_n, ms2_scan_n):
    text_trace = go.Scatterpolar(r=[1.1, 1.2, 0.74, 0.94],
                           theta=[-15, -30] ,
                           mode='text',
                           text=['MS1 Scans in {} minutes: {}'.format(params.LC_time, ms1_scan_n),
                                 'MS2 Scans in {} minutes: {}'.format(params.LC_time, ms2_scan_n),],
                           showlegend=False,
                           textfont={"size": [15]*2},
                           textposition='bottom right',
                           hoverinfo='text'
                           )
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
                                        hoverinfo='text'
                                        )
    return trace

def get_cycle_grid_tr():
    grid_trace = go.Scatterpolar(r=[0.5] * 50 + [0.7] * 50 + [0.9] * 50,
                    theta=list(np.linspace(90, 450, 50)) * 3,
                    mode="lines",
                    line={'width': 1, 'color':'#cccccc'},
                    showlegend=False,
                    hoverinfo='skip',)
    return grid_trace