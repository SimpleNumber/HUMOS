#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is web-app frontend
'''

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import mechanics, params
from dash.dependencies import Input, Output, State

#fixed data for resolution graph and space-charge effect graph
tmt_spectrum =  np.array([[127.12476, 1],[127.13108, 2]])
agc_spectrum = np.array([[1277.13108, 1]])

#generate DataFrame with theoretical ion currents for all models
ion_data = mechanics.get_ion_data(params.peptide_collection_size)
mechanics.normalize_ion_currents(ion_data, params.TIC, params.low_mass, params.high_mass)
boxes = mechanics.get_boxes(params.low_mass, params.high_mass, params.nBoxes, params.nScans, params.box_overlap)
mechanics.add_boxes(ion_data, boxes)

#interface
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
                meta_tags=[{'name': 'robots',
                            'content': 'noindex, nofollow'}])

app.title = 'HUMOS'
app.layout = html.Div([
    #header part
    html.Div([
        html.H1('HUMOS: How to Understand My Orbitrap Spectrum?', style={'flex-grow': '1'}),
        html.Img(src=mechanics.load_image('./helps/humos_logo.png'),
                 style={'height': '80px', 'padding-left': '2rem', 'padding-right': '2rem',
                        'transform': 'rotate(-10deg) skewY(4deg)'}),
             ], style={'display': 'flex'}),
    
    #AGC info table
    html.Div(dash_table.DataTable(
            id='table',
            style_cell_conditional=[{
                                     'if': {'column_id': ''},
                                     'fontWeight': 'bold'
                                    }],
            style_header={
                    'backgroundColor': 'white',
                    'fontWeight': 'bold',
                    'borderBottom': '1px solid black'
                         },
                                                        
            style_cell={
                    'boxShadow': '0 0',
                    'minWidth': '180px',
                       },
            ),
            style={ 'display':'flex', 'font-size':'16px', 'font':'CorierNew'}),
    
    #simulated mass spectrum        
    dcc.Graph(id='main-graph'),   
    
    #model parameters switches
    html.Div([
            html.Div([
                    html.H5('Peptide distribution'),
                    html.Div(dcc.RadioItems(
                            id='distribution',
                            options=[
                                    {'label': 'Equimolar', 'value': 'equal'},
                                    {'label': 'Regular', 'value': 'lognormal'},
                                    {'label': 'Regular with majors', 'value': 'lognormal-major'}
                                    ],
                            value='lognormal', 
                            ),style={'width': '80%','padding-left':'3%', 'padding-right':'10%',
                                    'font':'CorierNew'}),
                    ], style={'width':'400px'}),
                    
            html.Div([
                    html.H5('Resolution'),
                    html.Div([dcc.Slider(    
                            id='resolution-slider',
                            min=0,
                            max=len(params.resolutions_list)-1,
                            value=2,
                            marks={i: str(resolution) for i,resolution in enumerate(params.resolutions_list)},
                            step=1,
                            ),], style={'width': '80%','padding-left':'3%', 'padding-right':'10%', 'padding-bottom':'1%'}),
                    html.Br(),
                    html.H5('AGC Target'),
                    html.Div([dcc.Slider(
                                id='AGC-slider',
                                min=0,
                                max=len(params.agc_list)-1,
                                value=2,
                                marks={i: '{:.0e}'.format(agc) for i, agc in enumerate(params.agc_list)},
                                )], style={'width': '80%','padding-left':'3%', 'padding-right':'10%', 'padding-bottom':'1%'}),
                    html.Br(),
                    
                            html.H5('Max Injection Time (ms)'),
                            dcc.Input(id='mit-box', type='number',size='25', value=100),
                            html.Button('set', id='it-button'),

                    ],style={'width':'400px'}),
            
             html.Div([
                    html.H5('Acquiring Method'),
                    html.Div(dcc.RadioItems(id='method-choice',
                            options=[
                                    {'label': 'BoxCar', 'value': 'bc'},
                                    {'label': 'Usual MS1', 'value': 'ms1'},
                                    ],
                            value='ms1', 
                            ), ),
                    html.H5('TopN'),
                    html.Div([dcc.Slider(
                                id='topN-slider',
                                min=1,
                                max=40,
                                value=15,
                                marks={5*i: '{}'.format(5*i) for i in range(1,9)},
                                tooltip={i: "top" for i in range(1, 41)},
                                )], style={'width': '80%','padding-left':'3%', 'padding-right':'10%'}),
                    html.Br(),
                    html.Div([
                                html.P(id='cycletime'),
                                html.P(id='ms1-scan-n'),
                                html.P(id='ms2-scan-n')
                            
                            ], style={'width': '80%','padding-left':'3%', 'padding-right':'1%'})
                    ], style={'width':'400px'}),       
                             
            ], style={ 'display':'flex', 'flex-wrap': 'wrap', 'padding-bottom': '4rem'}),
    
    #smaller graphs
    html.Div([
            html.Div([
                    html.Center([
                            html.H5('Mass Spectral Resolution'),
                            html.P('The graph shows two adjacent TMT 10-plex reporter ions',
                                   style={'font-style': 'italic'}),
                            dcc.Graph(id='resolution-graph')
                            ]),
                    ],
                    style={'width':'600px', 'height':'525px'}),
            html.Div([
                    html.Center([
                            html.H5('AGC influence on mass accuracy'),
                            html.P('No calibration for space-charge effect applied',
                                   style={'font-style': 'italic'}),
                            dcc.Graph(id='accuracy-graph')
                            ])
                    ], 
                    style={'width':'600px', 'height':'525px'}),
            
            
        ],style={ 'display':'flex', 'flex-wrap': 'wrap'}),

    #footer part            
    html.Div([
        html.Img(src=mechanics.load_image('./helps/sdu_logo.png'),
                    style={'height': '30px', 'padding-top': '4rem'}),
            ], style={'textAlign': 'center'}),
    html.Div([
        html.P('Department of Biochemistry and Molecular Biology, University of Southern Denmark'),
        html.P(['Do you have any questions and suggestions about HUMOS? Contact us via ',
			    html.A('Github', href='https://github.com/SimpleNumber/HUMOS'), 
			   u' write to vgor (\u0430t) bmb.sdu.dk or juliabubis (\u0430t) gmail.com'
               ])
            ], style={'textAlign': 'center'}),
        ], style={'marginLeft':25, 'marginRight':25, 'marginBottom': 25, 'marginTop': 25}
    )

def get_zoom(relayout_data, min_x, max_x, min_y, max_y):
    #reads and preserves zoom information from a plotly graph
    x_range = []
    y_range = []
    if 'xaxis.range[0]' in relayout_data.keys():
        x_range =[relayout_data['xaxis.range[0]'],
                    relayout_data['xaxis.range[1]']
                ]
    else:
        x_range = [min_x, max_x]
    if 'yaxis.range[0]' in (relayout_data.keys()):
            y_range= [
                    relayout_data['yaxis.range[0]'],
                    relayout_data['yaxis.range[1]']
                ]
    else:
        y_range = [min_y, max_y]
    return x_range, y_range

def update_figure(selected_resolution, selected_agc, distribution, mit_clicked, 
                  method, relayout_data, max_it, topN):
    #respond to changes in modeling parameters
    boxCar = (method == 'bc')
    resolution = params.resolutions_list[selected_resolution]
    agc = params.agc_list[selected_agc]
    centroid_spectrum, real_st, real_agc = mechanics.get_full_spectrum(ion_data, distribution, agc, max_it)
    real_agcs = [real_agc]
    real_sts = [real_st]
    main_spectrum = mechanics.get_profile_spectrum(centroid_spectrum, resolution)
    
    if relayout_data == None:
        x_range = [min(centroid_spectrum[:,0]), max(centroid_spectrum[:,0])]
        y_range = [0, max(main_spectrum[1])]
    else:
        x_range, y_range = get_zoom(relayout_data,
                                    min(centroid_spectrum[:,0]),
                                    max(centroid_spectrum[:,0]),
                                    0,
                                    max(main_spectrum[1])
                                    )
        
    main_traces = [go.Scatter(x=main_spectrum[0], y=main_spectrum[1],
                              name='MS1 spectrum'
                              ),]
    labels_bc = []
    if boxCar:
        bc_spectra = mechanics.get_boxcar_spectra(ion_data, distribution,
                                                       agc, max_it, params.nBoxes, params.nScans)
        labels_bc = ['BoxCar scan {}'.format(i) for i in range(1, len(bc_spectra) +1 )]
        for bc_index, bc_label in zip(range(len(bc_spectra)), labels_bc):
            bc_spectrum = mechanics.get_profile_spectrum(bc_spectra[bc_index][0], resolution)
            main_traces.append(go.Scatter(x=bc_spectrum[0],y=bc_spectrum[1], 
                                          name=bc_label))
            real_agcs.append(bc_spectra[bc_index][2])
            real_sts.append(bc_spectra[bc_index][1])
        
    table = mechanics.make_table(real_sts, real_agcs, ['MS1'] + labels_bc, resolution)
    resolution_spectrum = mechanics.get_profile_spectrum(tmt_spectrum, resolution, points=51)
    resolution_traces = [go.Scatter(x=resolution_spectrum[0],
                                    y=resolution_spectrum[1],
                                    name=' '.join(['R =', str(resolution)])),
                         go.Scatter(x=[tmt_spectrum[0,0], tmt_spectrum[0,0]],
                                    y=[0, tmt_spectrum[0,1]],
                                    text='TMT 127N', 
                                    name='',
                                    mode='lines'
                                    ),
                         go.Scatter(x=[tmt_spectrum[1,0], tmt_spectrum[1,0]],
                                    y=[0, tmt_spectrum[1,1]],
                                    text='TMT 127C', 
                                    name='',
                                    mode='lines')]
    agc_spectrum_theoretical = mechanics.get_profile_spectrum(agc_spectrum, resolution, points=51)
    agc_mass_experimental = list(map(lambda x: [mechanics.charge_space_effect(x[0], agc), x[1]],
                                         agc_spectrum))
    agc_spectrum_experimental = mechanics.get_profile_spectrum(agc_mass_experimental, resolution, points=51)
    agc_traces = [go.Scatter(x=agc_spectrum_theoretical[0], 
                             y=agc_spectrum_theoretical[1],
                             text ='theoretical spectrum', 
                             name=''),
                  go.Scatter(x=[agc_spectrum[0,0], agc_spectrum[0,0]],
                             y=[0, agc_spectrum[0,1]],
                             text='{}'.format(agc_spectrum[0,0]), 
                             mode='lines',
                             name='',
                             line= {'color':'#1f77b4'},
                             ),

                  go.Scatter(x=agc_spectrum_experimental[0], 
                             y=agc_spectrum_experimental[1],
                             text ='experimental spectrum', 
                             line= {'color':'#ff7f0e'},
                             name=''),]
    
    return [
    [{"name": i, "id": i } for i in table.columns],
          table.to_dict('records'),

            {
        'data': main_traces,
        'layout': go.Layout(
                showlegend=True,
                margin={'t':30},
            xaxis={'title': 'm/z', 'range':x_range},
            yaxis={'title': 'Intensity', 'range': y_range}
        )
    
    },
    {
        'data': resolution_traces,
        'layout': go.Layout( 
                            margin={'t':10},
                            showlegend=False,
                            xaxis={'title': 'm/z'},
                            yaxis={'title': 'Intensity'},
        )
    
    },
     {
        'data': agc_traces,
        'layout': go.Layout(
                        margin={'t':10},
                            showlegend=False,
                            xaxis={ 'title': 'm/z',
                                   'range':[agc_spectrum[0,0]-0.15,agc_spectrum[0,0]+0.25,]},
                                   
                            yaxis={'title': 'Intensity'},
        )
    
    }
    ]

def update_ms_counts(topN, method, data, selected_resolution ):
    #update only counts of MS spectra, i.e. no changes to main spectrum applied
    boxCar = (method == 'bc')
    resolution = params.resolutions_list[selected_resolution]
    if data == None:
       return 'Select topN', '', ''
    else:
        data = pd.DataFrame(data)
        data = data.iloc[:, 1:].apply(pd.to_numeric)
        if boxCar:
            cycletime, ms1_scan_n, ms2_scan_n = mechanics.get_MS_counts('boxcar', data.iloc[0,:], 
                                                         topN, params.ms2length, params.LC_time, resolution)
        else:
            cycletime, ms1_scan_n, ms2_scan_n = mechanics.get_MS_counts('full', data.iloc[0,0], topN, 
                                                             params.ms2length, params.LC_time, resolution)
            
    return  'MS Cycle length: {:.3f} sec'.format(cycletime * 1e-3),\
            'MS1 Scans in {} minutes: {}'.format(params.LC_time, ms1_scan_n),\
            'MS2 Scans in {} minutes: {}'.format(params.LC_time, ms2_scan_n)

app.callback(
    [Output('table', 'columns'),
     Output('table', 'data'),
     Output('main-graph', 'figure'), 
     Output('resolution-graph', 'figure'), 
     Output('accuracy-graph', 'figure'),],  
    [Input('resolution-slider', 'value'), 
     Input('AGC-slider', 'value'),
     Input('distribution', 'value'),
     Input('it-button', 'n_clicks'),
     Input('method-choice', 'value')],    
     [State('main-graph', 'relayoutData'),
      State('mit-box', 'value'),
      State('topN-slider', 'value')])(update_figure)

app.callback(
    [Output('cycletime', 'children'),
     Output('ms1-scan-n', 'children'),
     Output('ms2-scan-n', 'children')],    
    [Input('topN-slider', 'value'), 
     Input('method-choice', 'value'),
     Input('table','data'),
     Input('resolution-slider', 'value')])(update_ms_counts)

server = app.server

if __name__ == '__main__':
    app.run_server()
