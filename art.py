#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file contains supportive/backend functions related to visual
representation
"""
import math
import params
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from colorsys import rgb_to_hsv, hsv_to_rgb
from plotly.colors import convert_colors_to_same_type
from plotly.express.colors import qualitative

def get_zoom(relayout_data, min_x, max_x, min_y, max_y):
    '''
    Reads and preserves zoom information from a plotly graph
    '''
    x_range = []
    y_range = []
    
    if 'xaxis.range[0]' in relayout_data.keys():
        x_range = [ relayout_data['xaxis.range[0]'],
                    relayout_data['xaxis.range[1]'] ]
    else:
        x_range = [min_x, max_x]
        
    if 'yaxis.range[0]' in (relayout_data.keys()):
            y_range= [ relayout_data['yaxis.range[0]'],
                       relayout_data['yaxis.range[1]'] ]
    else:
        y_range = [min_y, max_y]
        
    return x_range, y_range

def make_table(real_ats, real_agcs, labels, resolution):
    '''
    Create a table with acquisition parameters 

    Parameters
    ----------
    real_ats : list
        ion accumulation times per scan.
    real_agcs : list
        number of collected ions per scan.
    labels : list
        labels for scans.
    resolution : int
        used resolution.

    Returns
    -------
    df : pandas.DataFrame
        acquisition parameters table.

    '''
    real_sts = [max(acc_time, params.transients[resolution]) for acc_time in real_ats]
    df = pd.DataFrame([real_ats, real_agcs, real_sts], index = ["AT", "AGC", "ST"])
    df.loc['ST', :] = df.loc['ST', :].map('{:.2f}'.format)
    df.loc['AT', :] = df.loc['AT', :].map('{:.2f}'.format)
    df.loc['AGC', :] = df.loc['AGC', :].map('{:.1e}'.format)
    df.columns = labels
    df.insert(0, ' ', ['Ion accumulation time, ms', 'Accumulated charges', 'Scan time, ms'])
    return df

def tabletodf(data):
    '''
    Parse the table from HTML components format to pandas.DataFrame

    Parameters
    ----------
    data : dict
        table structure as returned by Dash, has to be `Table` type.

    Raises
    ------
    Exception
        if the type of element is not Table.

    Returns
    -------
    pandas.DataFrame
        representation of Dash table.

    '''
    #helper functions to parse Dash Table element
    def getContent(row):
        content = []
        for child in row['props']['children']:
            if child['type'] == 'Td' or child['type'] == 'Th':
                content.append(child['props']['children'])
        
        return content
    
    def getRows(data):
        rows = []
        for child in data['props']['children']:
            if child['type'] == 'Tr':
                rows.append(getContent(child))
        
        return rows
    
    #function body
    if data['type'] == 'Table':
        for child in data['props']['children']:
            if child['type'] == 'Thead':
                headers = getContent(child['props']['children'][0])
            elif child['type'] == 'Tbody':
                data = getRows(child)
                
        return pd.DataFrame(data, columns=headers)
    else:
        raise Exception("Not a Table")
    
def lightening_color(rgb_color):
    '''
    Lighten the color tone

    Parameters
    ----------
    rgb_color : str
        string representation of color in the following format
        'rgb(r, g, b)', where r, g, b are integers from 0 to 255.

    Returns
    -------
    str
        string representation of lightened color in the same format
        'rgb(r, g, b)', where r, g, b are integers from 0 to 255.

    '''
    r, g, b = [int(i) / 255 for i in rgb_color[4:-1].split(',')]
    hsv_color = list(rgb_to_hsv(r,g,b))
    hsv_color[1] *= 0.5 #desaturate 50%
    
    if hsv_color[1] == 0: # gray tones (magic stuff)
        hsv_color[2] = min(1.0, hsv_color[2] * 1.7) 
    else:
        hsv_color[2] = min(1.0, hsv_color[2] * 1.2)
        
    r, g, b = [int(i * 255) for i in hsv_to_rgb(*hsv_color)]
    
    return 'rgb({}, {}, {})'.format(r, g, b)

def get_colors(n_scans):
    '''
    Create color palette used in the tool

    Parameters
    ----------
    n_scans : int
        number of BoxCar Scans.

    Returns
    -------
    colors : list
        color codes in tuple type.

    '''
    colors = ['rgb(171, 226, 251)', qualitative.Dark2[-1], qualitative.D3[0]]
    
    #colors forBoxCar plots  
    c = qualitative.D3[1:3] + qualitative.Antique
    additional =  c * (n_scans // len(c)) + c [:n_scans % len(c)] #cycling palette
    colors += additional

    return convert_colors_to_same_type(colors)[0]

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
                     margin={'t': 30,
                             'l': 50},
                     xaxis={'title': 'm/z', 'range':x_range},
                     yaxis={'title': 'Abundance', 'range': y_range},
                     legend={'yanchor': 'top', 'y': 0.99,
                             'xanchor': 'right', 'x':0.99})

def get_dynrange_trace(row):
    '''
    Single trace in dynamic range plot
    '''
    return go.Scatter(x=row['x'],
                      y=row['y'],
                      line={'width': 7, 'color':row['color']},
                      mode='lines+text',
                      text=['{:.2f}'.format(np.log10(row['x'][0] / row['x'][1])), row['text']],        
                      textposition=['middle right', 'middle left']) 

def get_dynrange_layout(dr_df):
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
                            'range': [-1, len(dr_df)] },
                     showlegend=False,
                     width= 400,
                     height=180,
                     hovermode=False)

def get_obsPep_trace(observed_peptides, observed_color, missing_color):
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
                   marker_color=missing_color) ]

def get_obsPep_layout():
    '''
    Observed peptides plot layout
    '''
    return go.Layout(margin={'t': 10,
                             'l': 40,
                             'r': 10,
                             'b': 10},
                     xaxis={'visible': False},
                     yaxis={'title': '% detected peptides',
                            'range': [0, 100]},
                     showlegend=False,
                     barmode='stack',
                     width=100,
                     height=180,
                     hovermode=False)

def get_range(theta_range, step=1):
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
    
def get_cycle_grid():
    '''
    Circular grid of cycle time plot with annotations
    '''
    return [ go.Scatterpolar(r=[0.5] * 120 + [0.7] * 120 + [0.9] * 120,
                             theta=np.concatenate([np.linspace(0, 360, 120)] * 3),
                             mode='lines',
                             line={'width': 1,
                                   'color': '#cccccc'},
                             showlegend=False,
                             hoverinfo='skip'),
        
             go.Scatterpolar(r=[0.57, 0.77, 0.97],
                             theta=[0] * 3,
                             mode='text',
                             text=['Ion Trap', 'Orbitrap', 'Ion Accumulation'],
                             showlegend=False,
                             textfont={'size': [11] * 3},
                             textposition='middle center',
                             hoverinfo='skip') ]

def get_cycle_texts(cycletime, topN, ms1_scan_text, ms2_scan_text):
    '''
    Text objects located at the cycle plot

    Parameters
    ----------
    cycletime : float
        length of duty cycle
    topN : int
        number of MS2 spectra per cycle
    ms1_scan_text : string
        Text with MS1 scan information.
    ms2_scan_text : string
        Text with MS2 scan information.
    '''
    
    #calculte the location of text elements
    theta1 = 130
    theta2 = 136
    r1 = 1.05
    r2 = r1 * math.sin(theta1 * math.pi / 180) / math.sin(theta2 * math.pi / 180)
    
    text_trace = go.Scatterpolar(r=[0.07, 0.07, r1, r2],
                                 theta=[0, 180, theta1, theta2],
                                 mode='text',
                                 text= ['{:.3f} sec'.format(cycletime/1000),
                                        '#MS2: {}'.format(topN),
                                        ms1_scan_text, ms2_scan_text],
                                 showlegend=False,
                                 textfont={'size': [18] * 2 + [14] * 2},
                                 textposition=['middle center'] * 2 + ['bottom right'] * 2,
                                 hoverinfo='skip')
    return text_trace
    
def get_cycle_trace(row):
    '''
    Single trace in a circular plot
    '''
    return go.Scatterpolar(r=[row['r']] * len(row['theta']),
                           theta=row['theta'] ,
                           mode=row['mode'],
                           text=row['text'],
                           name=row['name'],
                           showlegend=row['showlegend'],
                           line={'width': row['line_width'],
                                 'color': row['line_color']},
                           textposition='bottom right',
                           hoverinfo='text')

def get_cycle_layout():
    '''
    Cycle Time plot layout
    '''
    return go.Layout(polar={'radialaxis': {'visible': False},
                            'angularaxis': {'rotation': 90, #start at the top
                                            'direction': 'clockwise', #increase angle clockwise
                                            'visible': False} },
                     showlegend=True,
                     legend={'x': 0.95,
                             'y': 0.85},
                     margin={'l': 0,
                             'r': 200,
                             'b': 0,
                             't': 0,
                             'pad': 0})

def get_ppp_trace(tX, tY, tC, sX, sY, sC, peptide):
    '''
    Data traces for points-per-peak plot
    
    Parameters
    ----------
    tX : `np.array`
        theoretical trace, x-array
    tY : `np.array`
        theoretical trace, y-array
    tC : string
        theoretical trace, color
    sX : `np.array`
        observed trace, x-array
    sY : `np.array`
        observed trace, y-array
    sC : string
        observed trace, color
    peptide : `pd.Series`
        `ion_data` element for the most abundant peak in all mass traces
    '''
    texts = [peptide['sequence'],
             'm/z {:.2f} {}+'.format(peptide['mz'], peptide['z'])]
    
    return [go.Scatter(x=tX, #theoretical points
                       y=tY,
                       mode='lines',
                       name='Theoretical',
                       fill='tozeroy',
                       line={'color': tC}),
            
            go.Scatter(x=sX, #sampling points
                       y=sY,
                       mode='lines+markers',
                       name='Observed',
                       fill='tozeroy',
                       line={'color': sC}),
            
            go.Scatter(x=[19.5, 19.5], #peptide info
                       y=[0.64, 0.58],
                       mode='text',
                       showlegend=False,
                       text=texts,
                       textposition='bottom left',
                       textfont={'size': 14}),
                       
            go.Scatter(x=[19.5], #points per peak info
                       y=[0.75],
                       mode='text',
                       showlegend=False,
                       text="Points/peak: {}".format(sX.shape[0] - 2), # explude border points
                       textposition='top left',
                       textfont={'size': 14}) ]

def get_ppp_layout():
    '''
    Points-per-peak plot layout
    '''
    return go.Layout(showlegend=True,
                     legend={'xanchor': 'right', 'x':0.99,
                             'yanchor': 'top', 'y': 0.99,},
                     margin={'l': 0,
                             'r': 0,
                             'b': 40,
                             't': 20},
                     xaxis={'range': [0, 20],
                            'title': 'RT (s)',
                            'showticklabels': False,
                            'zeroline': False},
                     yaxis={'range': [0, 1.01],
                            'ticks': ''},
                     hovermode=False)
