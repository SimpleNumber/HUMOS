#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This is web-app frontend
'''

import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import numpy as np
import pandas as pd
import mechanics, params, toolTips
import traces_layouts as tl
from plotly.express.colors import qualitative
from dash.dependencies import Input, Output, State

#color scheme
colors = mechanics.get_colors(params.nScans)

#fixed data for resolution graph and space-charge effect graph
tmt_spectrum =  np.array([[127.12476, 1],[127.13108, 2]])
agc_spectrum = np.array([[1277.13108, 1]])

#generate DataFrame with theoretical ion currents for all models
ion_data = mechanics.get_ion_data(params.peptide_collection_size)
mechanics.normalize_ion_currents(ion_data, params.low_mass, params.high_mass)
boxes = mechanics.get_boxes(params.low_mass, params.high_mass, params.nBoxes, params.nScans, params.box_overlap)
mechanics.add_boxes(ion_data, boxes)

#interface
block_style = {'width':'400px'}
small_panel_style = {'width': '80%','padding-left':'1%', 'padding-right':'10%', 'margin-top':'1rem', 'margin-bottom':'1rem'}
big_panel_style = { 'display':'flex', 'flex-wrap': 'wrap', 'padding-bottom': '4rem', 'justify-content': 'space-around'}
res_figure_style = {'width':'500px', 'height':'425px', 'padding-bottom': '4rem'}
cycle_figure_style = {'width':'650px', 'height':'425px', 'padding-bottom': '4rem'}
info_style = {'height': '15px', 'padding-bottom':'5px', 'padding-left':'3px', 'display':'inline'}
i_src = '/assets/info.png'

def table_dynRange_html():
 #AGC info table and Dynamic range graph
     return html.Div([
                html.Div([
                        html.H6('Information table', id='table-header'),
                        html.Img(id='i-table', src=i_src, style=info_style),
                        html.Div(id='table')
                        ], style={'height': '240px'}),
                html.Div([
                        html.H6('Dynamic range', id='dynamic-range-header'),
                        html.Img(id='i-dynamic-range', src=i_src, style=info_style),
                        html.Div([
                                dcc.Graph(id='dynamic-range-bar', config={'displayModeBar': False}),
                                dcc.Graph(id='observed-peptides', config={'displayModeBar': False})
                                ],
                                style={'display':'flex', 'flex-wrap': 'wrap'}
                                ),
    
                        ],
                    style={'height': '240px'}),
                toolTips.text_tooltip(toolTips.table_descript, 'table-header'),
                toolTips.text_tooltip(toolTips.table_descript, 'i-table'),
                toolTips.text_tooltip(toolTips.dynRange_descript, 'dynamic-range-header'),
                toolTips.text_tooltip(toolTips.dynRange_descript, 'i-dynamic-range'),
                toolTips.text_tooltip(toolTips.obsPeptides_descript, 'observed-peptides')],
                style={'display':'flex',
                       'flex-wrap': 'wrap',
                       'padding-bottom': '0rem',
                       'padding-left':'1%',
                       'padding-right':'10%',
                       'justify-content': 'space-between'})

def block1_html():
    #Block1 distribution and Acquisition method
    return html.Div([
                    html.H6('Peptide distribution', id='peptide-distr-header'),
                    html.Img(id='i-peptide-distr', src=i_src, style=info_style),
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
                    html.H6('Ion Current (charge/sec)', id='ion-current-header'),
                    html.Img(id='i-ion-current', src=i_src, style=info_style),
                    html.Div(dcc.Slider(
                                id='ionFlux',
                                min=0,
                                max=len(params.TIC) - 1,
                                value=2,
                                marks={i: '{:1.0e}'.format(v) for i, v in enumerate(params.TIC)},
                                step=1),
                            style=small_panel_style
                            ),
                    html.H6('Acquisition Method', id='acquisition-meth-header'),
                    html.Img(id='i-acquisition-meth', src=i_src, style=info_style),
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
                    toolTips.text_tooltip(toolTips.pep_distr_descript, 'i-peptide-distr'),
                    toolTips.text_tooltip(toolTips.ionCurrent_discript, 'ion-current-header'),
                    toolTips.text_tooltip(toolTips.ionCurrent_discript, 'i-ion-current'),
                    toolTips.text_tooltip(toolTips.acquisition_discript, 'acquisition-meth-header'),
                    toolTips.text_tooltip(toolTips.acquisition_discript, 'i-acquisition-meth'),
                    ],
                    style=block_style)

def block2_html():
    #Block2 MS1 parameters
    return html.Div([
                    html.H6('MS1 Resolution', id='MS1-resolution-header'),
                    html.Img(id='i-ms1-res', src=i_src, style=info_style),
                    html.Div([dcc.Slider(
                            id='resolution-slider',
                            min=1,
                            max=len(params.resolutions_list) - 1,
                            value=4,
                            marks={i: str(resolution) for i, resolution in enumerate(params.resolutions_list)},
                            step=1)
                            ],
                            style=small_panel_style),

                    html.H6('MS1 AGC Target', id='MS1-AGC-header'),
                    html.Img(id='i-ms1-agc', src=i_src, style=info_style),
                    html.Div([dcc.Slider(
                                id='AGC-slider',
                                min=0,
                                max=len(params.agc_list)-1,
                                value=2,
                                marks={i: '{:.0e}'.format(agc) for i, agc in enumerate(params.agc_list)},
                                )
                             ],
                             style=small_panel_style),

                    html.H6('MS1 Max Injection Time (ms)', id='MS1-IT-header'),
                    html.Img(id='i-ms1-mit', src=i_src, style=info_style),
                    html.Div([
                        dcc.Input(id='mit-box', type='number',size='20', value=100),
                        html.Button('set', id='it-button'),
                        ],
                        style=small_panel_style
                        ),
                    toolTips.text_tooltip(toolTips.resolution_descript, 'MS1-resolution-header'),
                    toolTips.text_tooltip(toolTips.resolution_descript, 'i-ms1-res'),
                    toolTips.text_tooltip(toolTips.AGC_discript, 'MS1-AGC-header'),
                    toolTips.text_tooltip(toolTips.AGC_discript, 'i-ms1-agc'),
                    toolTips.text_tooltip(toolTips.IT_descript,'MS1-IT-header'),
                    toolTips.text_tooltip(toolTips.IT_descript,'i-ms1-mit'),
                    ],
                style=block_style)

def block3_html():
    #Block3 MS2 parameters
    return html.Div([
                    html.H6('MS2 Resolution', id='MS2-resolution-header'),
                    html.Img(id='i-ms2-resolution', src=i_src, style=info_style),
                    html.Div([dcc.Slider(
                            id='resolution-ms2-slider',
                            min=0,
                            max=len(params.resolutions_list) - 1,
                            value=2,
                            marks={i: str(resolution) for i,resolution in enumerate(params.resolutions_list)},
                            step=1,
                            ),], style=small_panel_style),

                    html.H6('MS2 Max Injection Time (ms)', id='IT-MS2-header'),
                    html.Img(id='i-ms2-mit', src=i_src, style=info_style),
                    html.Div([
                            dcc.Input(id='mit-ms2-box', type='number',size='20', value=30),
                            html.Button('set', id='it-ms2-button')
                            ],
                            style=small_panel_style
                        ),

                    html.H6('TopN', id='topN-header'),
                    html.Img(id='i-topN', src=i_src, style=info_style),
                    html.Div([dcc.Slider(
                                id='topN-slider',
                                min=1,
                                max=40,
                                value=15,
                                marks={5*i: '{}'.format(5*i) for i in range(1,9)},
                                )],
                        style=small_panel_style),
                    dcc.Checklist(id='paral-checklist',
                                             options=[{'label': 'Parallelization', 'value': 'on'},],
                                             value=['on'],
                                             labelStyle={'display': 'inline-block'},
                                             style={'padding-bottom': '1rem', 'display':'inline' }),
                    html.Img(id='i-paral', src=i_src, style=info_style,),

                    toolTips.text_tooltip(toolTips.resolutionMS2_descript,'MS2-resolution-header'),
                    toolTips.text_tooltip(toolTips.resolutionMS2_descript,'i-ms2-resolution'),
                    toolTips.text_tooltip(toolTips.IT_descript,'IT-MS2-header'),
                    toolTips.text_tooltip(toolTips.IT_descript,'i-ms2-mit'),
                    toolTips.text_tooltip(toolTips.topN_discript,'topN-header'),
                    toolTips.text_tooltip(toolTips.topN_discript,'i-topN'),
                    toolTips.text_tooltip(toolTips.parallel_descript,'i-paral'),
                    ], style=block_style)

def res_fig_html():
     #smaller graphs
    return html.Div([
                    html.Center([
                            html.H6('Mass Spectral Resolution'),
                            html.P('The graph shows two adjacent TMT 10-plex reporter ions',
                                   style={'font-style': 'italic'}),
                            dcc.Graph(id='resolution-graph')
                            ]),
                    ],
                    style=res_figure_style)

def cycle_time_html():
    #cycle time plot etc
    return html.Div([
                    html.Center([
                            html.H6('Cycle Time'),
                            dcc.Graph(id='cycle-time-graph')
                            ]),
                    ],
                    style=cycle_figure_style)

#Main window
app = dash.Dash(__name__,
                meta_tags=[{'name': 'robots',
                            'content': 'noindex, nofollow'}])

app.title = 'HUMOS'

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
    
    #upper part - info table, dynamic range plot, observed peptides
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

    #lower part - resolution plot, cycle time plot
    html.Div([
            res_fig_html(),
            cycle_time_html()
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

    mechanics.scale_ion_currents(ion_data, ionFlux)
    
    centroid_spectrum, real_st, real_agc, peptides, max_int, min_int = \
               mechanics.get_full_spectrum(ion_data, distribution, agc, max_it)
    
    dr_df = pd.DataFrame({'text':['Peptide'], 
                          'x': [[ion_data['ic_' + distribution].max(),
                                 ion_data['ic_' + distribution].min()]],
                          'color': colors[1]})
    
    #Check if MS1 spectrum was empty
    if max_int > 0 and min_int > 0:
        dr_df = dr_df.append(pd.DataFrame({'text': 'MS1',
                                           'x': [[max_int, min_int]],
                                           'color': colors[2]}),
                               ignore_index=True)
    
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
        labels_bc = ['BoxCar scan {}'.format(i) for i in range(1, params.nScans + 1)]
        dr_df = dr_df.append(pd.DataFrame({'text': labels_bc,
                                           'x':list(bc_spectra[:, 4:6]),
                                           'color':colors[3:]}),
                                         ignore_index=True)

        for bc_index, bc_label in enumerate(labels_bc):
            bc_spectrum = mechanics.get_profile_spectrum(bc_spectra[bc_index][0], resolution)
            main_traces.append(go.Scatter(x=bc_spectrum[0],
                                          y=bc_spectrum[1],
                                          name=bc_label))
            
            real_agcs.append(bc_spectra[bc_index][2])
            real_sts.append(bc_spectra[bc_index][1])
            peptides.update(bc_spectra[bc_index][3])

            #check if BoxCar spectrum was empty
            if bc_spectra[bc_index][4] > 0:
                max_int = max(max_int, bc_spectra[bc_index][4])
            
            if bc_spectra[bc_index][5] > 0:
                min_int = min(min_int, bc_spectra[bc_index][5])

    dr_df = dr_df.append(pd.DataFrame({'text': 'Spectrum',
                                       'x': [[max_int, min_int]],
                                       'color': colors[0]}),
                           ignore_index=True)

    observed_peptides = np.round(100 * len(peptides) / len(ion_data["sequence"].unique()), 1)

    table = mechanics.make_table(real_sts, real_agcs, ['MS1'] + labels_bc, resolution)
    
    dr_df['y'] = [[i,i] for i in range(0, len(dr_df))]
    dr_df['trace'] = dr_df.apply(tl.get_dr_tr, axis=1)
    dr_df.index = dr_df['text']

    obsPeptides_traces = tl.get_obsPep_tr(observed_peptides, colors[0], colors[1])    

    return [dbc.Table.from_dataframe(table),
            {'data': main_traces, 'layout': tl.get_main_layout(x_range, y_range)},
            {'data': list(dr_df['trace']), 'layout': tl.get_dr_layout(dr_df)},        
            {'data': obsPeptides_traces, 'layout': tl.get_obsPep_layout()}]
            
def update_ms_counts(topN, method, data, selected_resolution, ms2_resolution, 
                     parallel, mit_clicked,  mit_ms2 ):
    #update counts of MS spectra and cycle time graph, no changes to main spectrum applied
    boxCar = (method == 'bc')
    parallel = True if len(parallel) > 0 else False #value is a list of selected options
    ms2_resolution = params.resolutions_list[ms2_resolution]
    resolution = params.resolutions_list[selected_resolution]
    
    if data == None:
       return 'Select topN', '', ''
    else:
        data = mechanics.tabletodf(data)
        data = data.iloc[:, 1:].apply(pd.to_numeric)
        if boxCar:
            cycletime, ms1_scan_n, ms2_scan_n, queues = mechanics.get_MS_counts('boxcar', data.iloc[0,:],
                                                         resolution, topN, ms2_resolution, mit_ms2,
                                                         params.LC_time, parallel=parallel)

        else:
            cycletime, ms1_scan_n, ms2_scan_n, queues = mechanics.get_MS_counts('full', data.iloc[0,0], 
                                                         resolution, topN, ms2_resolution, mit_ms2,
                                                         params.LC_time, parallel=parallel)
    
    ms1_scan_text = 'MS1 Scans in {} minutes: {}'.format(params.LC_time, ms1_scan_n)
    ms2_scan_text = 'MS2 Scans in {} minutes: {}'.format(params.LC_time, ms2_scan_n)
    
    iat =['Accumulation ' + i for i in list(data.columns) + ['MS2 ' + str(j) for j in range(1, topN + 1)]]
    ot = ['Acquisition '+ i for i in list(data.columns)]
    ot_names = ['Acquisition ' + i for i in list(data.columns)]
    
    main_colors = colors[2: 2 + len(data.columns)]
    theta_start = 90 - pd.Series(queues['IS'][::2]) / cycletime * 360
    theta_end = 90 - pd.Series(queues['IS'][1::2]) / cycletime * 360
    theta_start = theta_start.append(90 - pd.Series(queues['OT'][::2]) / cycletime * 360, ignore_index=True)
    theta_end = theta_end.append(90 - pd.Series(queues['OT'][1::2]) / cycletime * 360, ignore_index=True)
    it = []
    it_names= []
    if len(queues['IT']) > 1:
        it = ['Acquisition MS2 ' + str(i) for i in range(1, topN + 1)]
        it_names = ['Acquisition MS2 '] * topN
        main_colors += [qualitative.Dark2[5]] * topN
        theta_start = theta_start.append(90 - pd.Series(queues['IT'][::2]) / cycletime * 360, ignore_index=True)
        theta_end = theta_end.append(90 - pd.Series(queues['IT'][1::2]) / cycletime * 360, ignore_index=True)
    else:
        ot += ['Acquisition MS2 ' + str(i) for i in range(1, topN + 1)]
        ot_names += ['Acquisition MS2 '] * topN
        main_colors += [qualitative.Dark2[-1]] * topN
#    print( main_colors)
    aqc_names = [ 'Accumulation ' + i for i in list(data.columns) + ['MS2 '] * topN]
    aqc_show_legend = [True] * len(data.columns) + [True] + [False] * (topN - 1)

#    print(len(iat + ot + it))
#    print(len(aqc_names + ot_names + it_names))
#    print(len(aqc_show_legend))
#    print(len([mechanics.lightening_color(i) for i in main_colors] + main_colors))
#    print(len(aqc_names + ot_names + it_names))
#    print(len(aqc_names + ot_names + it_names))
    
    cycle_df = pd.DataFrame({'text': iat + ot + it, 
                           'mode': 'lines',
                           'name': aqc_names + ot_names + it_names,
                           'hoverinfo': 'text',
                           'hoveron': 'fills',
                           'showlegend': aqc_show_legend * 2,
                           'r': [0.9] * len(iat) + [0.7] * len(ot) + [0.5] * len(it),
                           'line_color': [mechanics.lightening_color(i) for i in main_colors] + main_colors,
                           'line_width': 13,
                           'start': theta_start ,
                           'end': theta_end ,
                           })

    cycle_df['theta'] = cycle_df.loc[:,['start', 'end']].apply(tl.get_theta_ranges, axis=1)
    cycle_df['traces'] = cycle_df.apply(tl.get_cycle_tr, axis=1)
    cycle_traces = [tl.get_cycle_grid_tr()]
    cycle_traces += list(cycle_df['traces'])
    cycle_traces.append(tl.get_cycle_annotations_tr(cycletime))
    cycle_traces.append(tl.get_cycle_text_tr(ms1_scan_text, ms2_scan_text))    
    return  [{'data': cycle_traces,'layout': tl.get_cycle_layout()}]
                      

def update_resolution_graph(selected_resolution):
    #change only the resolution graph
    resolution = 1000 if selected_resolution == 0 else params.resolutions_list[selected_resolution]
    #ion trap have resolution ~1000 (0.1 Th) at 100
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

    return [ {'data': resolution_traces,
              'layout': go.Layout(
                            margin={'t':10},
                            showlegend=False,
                            xaxis={'title': 'm/z'},
                            yaxis={'title': 'Intensity'})}]

app.callback(
    [Output('table', 'children'),
     Output('main-graph', 'figure'),
     Output('dynamic-range-bar','figure'),
     Output('observed-peptides', 'figure')],
    [Input('resolution-slider', 'value'),
     Input('AGC-slider', 'value'),
     Input('distribution', 'value'),
     Input('it-button', 'n_clicks'),
     Input('method-choice', 'value'),
     Input('ionFlux', 'value')],
     [State('main-graph', 'relayoutData'),
      State('mit-box', 'value')])(update_figure)

app.callback(
    [
     Output('cycle-time-graph', 'figure')],
    [Input('topN-slider', 'value'),
     Input('method-choice', 'value'),
     Input('table','children'),
     Input('resolution-slider', 'value'),
     Input('resolution-ms2-slider', 'value'),
     Input('paral-checklist', 'value'),
     Input('it-ms2-button','n_clicks')],
    [State('mit-ms2-box', 'value')])(update_ms_counts)

app.callback(
    [Output('resolution-graph', 'figure')],
    [Input('resolution-ms2-slider', 'value')])(update_resolution_graph)

server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)
