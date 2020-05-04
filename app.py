#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is web-app frontend
'''

import dash
import dash_table
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import mechanics, params, toolTips
from dash.dependencies import Input, Output, State


#fixed data for resolution graph and space-charge effect graph
tmt_spectrum =  np.array([[127.12476, 1],[127.13108, 2]])
agc_spectrum = np.array([[1277.13108, 1]])

#generate DataFrame with theoretical ion currents for all models
ion_data = mechanics.get_ion_data(params.peptide_collection_size)
mechanics.normalize_ion_currents(ion_data, params.TIC[2], params.low_mass, params.high_mass)
boxes = mechanics.get_boxes(params.low_mass, params.high_mass, params.nBoxes, params.nScans, params.box_overlap)
mechanics.add_boxes(ion_data, boxes)

#default TIC
TIC = params.TIC[2]

#interface
#TODO license informrmation for CSS scheme?
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets,
                meta_tags=[{'name': 'robots',
                            'content': 'noindex, nofollow'}])

app.title = 'HUMOS'

block_style = {'width':'400px'}
small_panel_style = {'width': '80%','padding-left':'1%', 'padding-right':'10%'}
big_panel_style = { 'display':'flex', 'flex-wrap': 'wrap', 'padding-bottom': '4rem', 'justify-content': 'space-around'}
figure_style = {'width':'500px', 'height':'425px', 'padding-bottom': '4rem'}

def table_dynRange_html():
 #AGC info table and Dynamic range graph
     return html.Div([
            html.Div([
                    html.H5('Information table', id='table-header'),
                    #html.Div(id='table'),
                    dash_table.DataTable(
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
                            'minWidth': '120px',
                                },
                    )
                    ], style={'height': '250px'}),
            html.Div([
                    html.H5('Dynamic range', id='dynamic-range-header'),
                    html.Div([
                            dcc.Graph(id='dynamic-range-bar'),
                            dcc.Graph(id='observed-peptides'),
                            dbc.Tooltip( "Peptides observed in MS1 spectra",
                                        target="observed-peptides")
                            ],
                            style={'display':'flex', 'flex-wrap': 'wrap'}
                        ),

                    ],
                    style={'height': '250px'}),
            toolTips.text_tooltip(toolTips.table_descript, 'table-header'),
            toolTips.text_tooltip(toolTips.dynRange_descript, 'dynamic-range-header')],
            style={'display':'flex',
                   'font':'CorierNew',
                   'flex-wrap': 'wrap',
                   'padding-bottom': '0rem',
                   'padding-left':'1%',
                   'padding-right':'10%',
                   'justify-content': 'space-between'})

def block1_html():
    #Block1 distribution and Acquisition method
    return html.Div([
                    html.H5('Peptide distribution', id='peptide-distr-header'),
                    html.Div(dcc.RadioItems(
                            id='distribution',
                            options=[
                                    {'label': 'Equimolar', 'value': 'equal', },
                                    {'label': 'Regular', 'value': 'lognormal'},
                                    {'label': 'Regular with majors', 'value': 'lognormal-major'}
                                    ],
                            value='lognormal'),
                            style=small_panel_style
                               ),
                    html.H5('Ion Current (charge/sec)', id='ion-current-header'),
                    html.Div(dcc.Slider(
                                id='ionFlux',
                                min=0,
                                max=len(params.TIC) - 1,
                                value=2,
                                marks={i: '{:1.0e}'.format(v) for i, v in enumerate(params.TIC)},
                                step=1),
                            style=small_panel_style
                            ),
                    html.H5('Acquisition Method', id='acquisition-meth-header'),
                    html.Div(dcc.RadioItems(
                            id='method-choice',
                            options=[
                                    {'label': 'BoxCar', 'value': 'bc'},
                                    {'label': 'Usual MS1', 'value': 'ms1'},
                                    ],
                            value='ms1'),
                            style=small_panel_style
                            ),
                    toolTips.text_tooltip(toolTips.pep_distr_descript, 'peptide-distr-header'),

                    toolTips.text_tooltip(toolTips.ionCurrent_discript, 'ion-current-header'),

                    toolTips.text_tooltip(toolTips.acquisition_discript, 'acquisition-meth-header'),

                    ],
                    style=block_style)

def block2_html():
    #Block2 MS1 parameters
    return html.Div([
                    html.H5('MS1 Resolution', id='MS1-resolution-header'),
                    html.Div([dcc.Slider(
                            id='resolution-slider',
                            min=1,
                            max=len(params.resolutions_list) - 1,
                            value=4,
                            marks={i: str(resolution) for i, resolution in enumerate(params.resolutions_list)},
                            step=1)
                            ],
                            style=small_panel_style),

                    html.H5('MS1 AGC Target', id='MS1-AGC-header'),
                    html.Div([dcc.Slider(
                                id='AGC-slider',
                                min=0,
                                max=len(params.agc_list)-1,
                                value=2,
                                marks={i: '{:.0e}'.format(agc) for i, agc in enumerate(params.agc_list)},
                                )
                             ],
                             style=small_panel_style),

                    html.H5('MS1 Max Injection Time (ms)', id='MS1-IT-header'),
                    html.Div([
                        dcc.Input(id='mit-box', type='number',size='20', value=100),
                        html.Button('set', id='it-button'),
                        ],
                        style=small_panel_style
                        ),
                    toolTips.text_tooltip(toolTips.resolution_descript, 'MS1-resolution-header'),
                    toolTips.text_tooltip(toolTips.AGC_discript, 'MS1-AGC-header'),
                    toolTips.text_tooltip(toolTips.IT_descript,'MS1-IT-header'),
                    ],
                style=block_style)

def block3_html():
    #Block3 MS2 parameters
    return html.Div([
                    html.H5('MS2 Resolution', id='MS2-resolution-header'),
                    html.Div([dcc.Slider(
                            id='resolution-ms2-slider',
                            min=0,
                            max=len(params.resolutions_list) - 1,
                            value=2,
                            marks={i: str(resolution) for i,resolution in enumerate(params.resolutions_list)},
                            step=1,
                            ),], style=small_panel_style),

                    html.H5('MS2 Max Injection Time (ms)', id='IT-MS2-header'),
                    html.Div([
                            dcc.Input(id='mit-ms2-box', type='number',size='20', value=30),
                            html.Button('set', id='it-ms2-button')
                            ],
                            style=small_panel_style
                        ),

                    html.H5('TopN', id='topN-header'),
                    html.Div([dcc.Slider(
                                id='topN-slider',
                                min=1,
                                max=40,
                                value=15,
                                marks={5*i: '{}'.format(5*i) for i in range(1,9)},
                                tooltip={'placement': 'top'},
                                )],
                        style=small_panel_style),
                     html.Div(dcc.Checklist( id='paral-checklist',
                            options=[{'label': 'Parallelization', 'value': 'on'},],
                                    value=['on'],
                            style={'padding-bottom': '1rem'}), id='paral-checklist-target'),
                    html.Div([
                                html.P(id='cycletime', ),
                                html.P(id='ms1-scan-n'),
                                html.P(id='ms2-scan-n')

                            ], style=small_panel_style),
                    toolTips.text_tooltip(toolTips.resolutionMS2_descript,'MS2-resolution-header'),
                    toolTips.text_tooltip(toolTips.IT_descript,'IT-MS2-header'),
                    toolTips.text_tooltip(toolTips.topN_discript,'topN-header'),
                    toolTips.text_tooltip(toolTips.parallel_descript,'paral-checklist-target')
                    ], style=block_style)

def res_fig_html():
     #smaller graphs
    return html.Div([
                    html.Center([
                            html.H5('Mass Spectral Resolution'),
                            html.P('The graph shows two adjacent TMT 10-plex reporter ions',
                                   style={'font-style': 'italic'}),
                            dcc.Graph(id='resolution-graph')
                            ]),
                    ],
                    style=figure_style)

def AGC_fig_html():
    return html.Div([
                    html.Center([
                            html.H5('AGC influence on mass accuracy'),
                            html.P('No calibration for space-charge effect applied',
                                   style={'font-style': 'italic'}),
                            dcc.Graph(id='accuracy-graph')
                            ])
                    ],
                    style=figure_style)

app.layout = html.Div([
    #header part
    html.Div([
        html.H1('HUMOS: How to Understand My Orbitrap Spectrum?', style={'flex-grow': '1'}),
        html.Img(id='logo', src='/assets/humos_logo.png',
                 style={'height': '80px',
                        'padding-left': '2rem',
                        'padding-right': '2rem',
                        'transform': 'rotate(-10deg) skewY(4deg)'}),
            toolTips.logo_tooltip()
             ], style={'display': 'flex'}),
    table_dynRange_html(),
    #simulated mass spectrum
    dcc.Graph(id='main-graph'),

    #model parameters switches
    html.Div([
            #Block1 distribution and Acquisition
            block1_html(),
            #Block2 MS1 parameters
            block2_html(),
             #Block3 MS2 parameters
            block3_html(),
            ], style=big_panel_style),

    #smaller figures
    html.Div([
            res_fig_html(),
            AGC_fig_html()

        ],style=big_panel_style),

    #footer part
    html.Div([
        html.Img(src='/assets/sdu_logo.png',
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
                  method, ionFlux, relayout_data, max_it):
    #respond to changes in modeling parameters
    boxCar = (method == 'bc')
    resolution = params.resolutions_list[selected_resolution]
    agc = params.agc_list[selected_agc]
    ionFlux = params.TIC[ionFlux]

    global TIC

    if TIC != ionFlux:
        mechanics.scale_ion_currents(ion_data, ionFlux)
        TIC = ionFlux
    
    dyn_range = {}
    dyn_range['Peptide'] = [ion_data['ic_' + distribution].max(),
                            ion_data['ic_' + distribution].min()]

    centroid_spectrum, real_st, real_agc, peptides, max_int, min_int = \
        mechanics.get_full_spectrum(ion_data, distribution, agc, max_it)
    
    if max_int > 0 and min_int > 0:
        dyn_range['MS1'] = [max_int, min_int]

    real_agcs = [real_agc]
    real_sts = [real_st]
    main_spectrum = mechanics.get_profile_spectrum(centroid_spectrum, resolution)

    if relayout_data == None:
        x_range = [min(main_spectrum[0]), max(main_spectrum[0])]
        y_range = [0, max(main_spectrum[1])]
    else:
        x_range, y_range = get_zoom(relayout_data,
                                    min(main_spectrum[0]),
                                    max(main_spectrum[0]),
                                    0,
                                    max(main_spectrum[1]))

    main_traces = [go.Scatter(x=main_spectrum[0], y=main_spectrum[1],
                                  name='MS1 spectrum')]
    
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
            peptides.update(bc_spectra[bc_index][3])

            dyn_range[bc_label] = [bc_spectra[bc_index][4],
                                   bc_spectra[bc_index][5]]
            
            if bc_spectra[bc_index][4] > 0:
                max_int = max(max_int, bc_spectra[bc_index][4])
            
            if bc_spectra[bc_index][5] > 0:
                min_int = min(min_int, bc_spectra[bc_index][5])
    
    dyn_range['Spectrum'] = [max_int, min_int]

    observed_peptides = np.round(100 * len(peptides) / len(ion_data["sequence"].unique()),1)

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
    agc_mass_experimental = np.array([[mechanics.charge_space_effect(x[0], agc), x[1]]
                                         for x in agc_spectrum])
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
    
    dynRange_traces = [go.Scatter(x=drange,
                              y=[i, i],
                              line={'width': 7},
                              mode='lines+text',
                              text=['{:.2f}'.format(np.log10(drange[0] / drange[1])), label],
                              textposition=['middle right', 'middle left']
                             ) for i, (label, drange) in enumerate(dyn_range.items())]

    obsPeptides_traces = [go.Bar(x=[0],
                              y=[observed_peptides],
                              width=1,
                              orientation='v',
                              name='% observed peptides',
                              text=str(observed_peptides),
                              textposition='inside',
                              marker_color=['#a7e2f9']
                             ),
                          go.Bar(x=[0],
                              y=[100 - observed_peptides],
                              width=1,
                              orientation='v',
                              name='% missing peptides',
                              marker_color=['#0576b0']
                             )
                       ]

    return [
    #dbc.Table.from_dataframe(table, bordered=True, hover=True), 
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

    },
    {
        'data': dynRange_traces,
        'layout': go.Layout(
                        margin={'t':0,
                                'l':10,
                                'b': 40},
                        xaxis={'title': 'Abundance',
                               'type': 'log',
                               'exponentformat': 'power',
                               'range': [np.log10(dyn_range['Peptide'][1]) - 2,
                                         np.log10(dyn_range['Peptide'][0]) + 1]},
                        yaxis={'visible': False,
                               'range': [-1, len(dyn_range)]},
                        showlegend=False,
                        width= 400,
                        height=140,
                        hovermode=False
        )

    },
    {
        'data': obsPeptides_traces,
        'layout': go.Layout(
                        margin={'t':10,
                                'l':40,
                                'r':10,
                                'b': 10},
                        xaxis={'visible': False},
                        yaxis={'title': '% visible peptides',
                               'range': [0, 100]},
                        showlegend=False,
                        barmode='stack',
                        width=100,
                        height=140,
                        hovermode=False
        )

    }
    ]

def update_ms_counts(topN, method, data, selected_resolution, ms2_resolution, parallel, mit_clicked,  mit_ms2 ):
    #update only counts of MS spectra, i.e. no changes to main spectrum applied
    boxCar = (method == 'bc')
    parallel = True if len(parallel) > 0 else False
    ms2_resolution = params.resolutions_list[ms2_resolution]
    resolution = params.resolutions_list[selected_resolution]
    if data == None:
       return 'Select topN', '', ''
    else:
        # print(data.children)
        data = pd.DataFrame(data)
        data = data.iloc[:, 1:].apply(pd.to_numeric)
        if boxCar:
            cycletime, ms1_scan_n, ms2_scan_n = mechanics.get_MS_counts('boxcar', data.iloc[0,:],
                                                         topN, (mit_ms2, ms2_resolution), params.LC_time,
                                                         resolution, parallel=parallel)
        else:
            cycletime, ms1_scan_n, ms2_scan_n = mechanics.get_MS_counts('full', data.iloc[0,0], topN,
                                                             (mit_ms2, ms2_resolution), params.LC_time,
                                                             resolution, parallel=parallel)

    return  'MS Cycle length: {:.3f} sec'.format(cycletime * 1e-3),\
            'MS1 Scans in {} minutes: {}'.format(params.LC_time, ms1_scan_n),\
            'MS2 Scans in {} minutes: {}'.format(params.LC_time, ms2_scan_n)

app.callback(
    [Output('table', 'columns'),
     Output('table', 'data'),
     #Output('table', 'children'),
     Output('main-graph', 'figure'),
     Output('resolution-graph', 'figure'),
     Output('accuracy-graph', 'figure'),
     Output('dynamic-range-bar','figure'),
     Output('observed-peptides', 'figure')],
    [Input('resolution-slider', 'value'),
     Input('AGC-slider', 'value'),
     Input('distribution', 'value'),
     Input('it-button', 'n_clicks'),
     Input('method-choice', 'value'),
     Input('ionFlux', 'value')],
     [State('main-graph', 'relayoutData'),
      State('mit-box', 'value'),
      ])(update_figure)

app.callback(
    [Output('cycletime', 'children'),
     Output('ms1-scan-n', 'children'),
     Output('ms2-scan-n', 'children')],
    [Input('topN-slider', 'value'),
     Input('method-choice', 'value'),
     Input('table','data'),
     Input('resolution-slider', 'value'),
     Input('resolution-ms2-slider', 'value'),
     Input('paral-checklist', 'value'),
     Input('it-ms2-button','n_clicks')],
    [State('mit-ms2-box', 'value')])(update_ms_counts)

server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)
