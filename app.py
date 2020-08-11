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
import art, mechanics, params, tooltips
from plotly.express.colors import qualitative
from dash.dependencies import Input, Output, State

#color scheme
colors = art.get_colors(params.nScans)

#fixed data for resolution graph
tmt_spectrum =  np.array([[127.12476, 1],[127.13108, 2]])

#generate DataFrame with theoretical ion currents for all models
ion_data = mechanics.get_ion_data(params.peptide_collection_size)
mechanics.normalize_ion_currents(ion_data, params.low_mass, params.high_mass)
boxes = mechanics.get_boxes(params.low_mass, params.high_mass, params.nBoxes, params.nScans, params.box_overlap)
mechanics.add_boxes(ion_data, boxes)

### Interface building blocks ###
block_style = {'width':'400px'}
small_panel_style = {'padding-left': '2%', 'padding-right': '2%', 'margin-top': '1rem', 'margin-bottom': '1rem'}
big_panel_style = {'display': 'flex', 'flex-wrap': 'wrap', 'padding': '0 1% 2rem 1%', 'justify-content': 'space-between'}
main_graph_style ={'flex': '1 1 800px', 'min-width': '400px'}
res_figure_style = {'width': '600px', 'height': '450px', 'padding-bottom': '4rem'}
cycle_figure_style = {'width': '600px', 'height': '450px', 'padding-bottom': '4rem'}
info_style = {'height': '15px', 'padding-bottom': '5px', 'padding-left': '3px', 'display': 'inline'}
header_style = {'display': 'inline', 'font-size': '2rem', 'margin-bottom': '1rem', 'margin-top': '1rem'}
ppp_figure_style = {'width':'300px', 'padding-bottom': '4rem', 'padding-right': '1rem'}
i_src = '/assets/info.png'

def table_dynRange_html():
    '''
    AGC info table and Dynamic range graph
    '''
    return html.Div([
                html.Div([
                        html.H6('Information Table', id='table-header'),
                        html.Img(id='i-table', src=i_src, style=info_style),
                        html.Div(id='table')
                        ], style={'flex-grow': '1'}),
                html.Div([
                        html.H6('Dynamic Range', id='dynamic-range-header'),
                        html.Img(id='i-dynamic-range', src=i_src, style=info_style),
                        dcc.Graph(id='dynamic-range-bar', config={'displayModeBar': False})
                        ], style={'height': '240px'}),
                html.Div([
                        dcc.Graph(id='observed-peptides', config={'displayModeBar': False})
                        ],
                        id='i-observed-peptides',
                        style={'height': '210px', 'padding-top': '30px'}),
                
                tooltips.text_tooltip(tooltips.info_table, 'table-header'),
                tooltips.text_tooltip(tooltips.info_table, 'i-table'),
                tooltips.text_tooltip(tooltips.dynamic_range, 'dynamic-range-header'),
                tooltips.text_tooltip(tooltips.dynamic_range, 'i-dynamic-range'),
                tooltips.text_tooltip(tooltips.observed_peptides, 'i-observed-peptides')],
                style={'display':'flex',
                       'flex-wrap': 'wrap',
                       'padding-bottom': '2rem',
                       'padding-left':'1%',
                       'padding-right':'1%',
                       'justify-content': 'space-between'})

def block_global_html():
    '''
    Global parameters: Distribution, TIC, and Acquisition method
    '''
    return html.Div([
                    html.H6('Peptide Distribution', id='peptide-distr-header'),
                    html.Img(id='i-peptide-distr', src=i_src, style=info_style),
                    html.Div(dcc.RadioItems(
                            id='distribution',
                            options=[
                                    {'label': 'Equimolar', 'value': 'equal', },
                                    {'label': 'Regular', 'value': 'lognormal'},
                                    {'label': 'Regular with majors', 'value': 'lognormal-major'}
                                    ],
                            value='lognormal'),
                            style=small_panel_style),
    
                    html.H6('Total Ion Current (ion/sec)', id='ion-current-header'),
                    html.Img(id='i-ion-current', src=i_src, style=info_style),
                    html.Div(dcc.Slider(
                                id='ionFlux',
                                min=0,
                                max=len(params.TIC) - 1,
                                value=4,
                                marks={i: '{:1.0e}'.format(v) for i, v in enumerate(params.TIC)},
                                step=1),
                            style=small_panel_style),
    
                    html.H6('Acquisition Method', id='acquisition-meth-header'),
                    html.Img(id='i-acquisition-meth', src=i_src, style=info_style),
                    html.Div(dcc.RadioItems(
                            id='method-choice',
                            options=[
                                    {'label': 'BoxCar', 'value': 'bc'},
                                    {'label': 'Usual MS1', 'value': 'ms1'},
                                    ],
                            value='ms1'),
                            style=small_panel_style),
                    
                    tooltips.text_tooltip(tooltips.peptide_distribution, 'peptide-distr-header'),
                    tooltips.text_tooltip(tooltips.peptide_distribution, 'i-peptide-distr'),
                    tooltips.text_tooltip(tooltips.ion_current, 'ion-current-header'),
                    tooltips.text_tooltip(tooltips.ion_current, 'i-ion-current'),
                    tooltips.text_tooltip(tooltips.acquisition, 'acquisition-meth-header'),
                    tooltips.text_tooltip(tooltips.acquisition, 'i-acquisition-meth'),
                    ],
                    style=block_style)

def block_MS1_html():
    '''
    MS1-related parameters
    '''
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
                    
                    tooltips.text_tooltip(tooltips.resolution, 'MS1-resolution-header'),
                    tooltips.text_tooltip(tooltips.resolution, 'i-ms1-res'),
                    tooltips.text_tooltip(tooltips.AGC, 'MS1-AGC-header'),
                    tooltips.text_tooltip(tooltips.AGC, 'i-ms1-agc'),
                    tooltips.text_tooltip(tooltips.MaxIT, 'MS1-IT-header'),
                    tooltips.text_tooltip(tooltips.MaxIT, 'i-ms1-mit'),
                    ],
                    style=block_style)

def block_MS2_html():
    '''
    MS2-related parameters
    '''
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
                            ),
                            ], style=small_panel_style),

                    html.H6('MS2 Max Injection Time (ms)', id='IT-MS2-header'),
                    html.Img(id='i-ms2-mit', src=i_src, style=info_style),
                    html.Div([
                            dcc.Input(id='mit-ms2-box', type='number',size='20', value=30),
                            html.Button('set', id='it-ms2-button')
                            ],
                            style=small_panel_style),

                    dcc.RadioItems(id='topN-topSpeed-choice', 
                                   options=[
                                           {'label': 'TopN', 'value': 'topN'},
                                           {'label': 'TopSpeed', 'value': 'topSpeed'}
                                           ],
                                   value='topN',
                                   labelStyle=header_style,
                                   style={'display': 'inline'}),
                    html.Img(id='i-topN', src=i_src, style=info_style),
                    html.Div([dcc.Slider(id='topN-slider',
                                         min=0,
                                         max=40,
                                         value=15,
                                         marks={i: '{}'.format(i) for i in range(0, 41, 5)},
                                         tooltip={'placement': 'bottom'})],
                             style=small_panel_style),
                    
                    dcc.Checklist(id='paral-checklist',
                                  options=[{'label': 'Parallelization', 'value': 'on'},],
                                  value=['on'],
                                  labelStyle={'display': 'inline-block'},
                                  style={'padding-bottom': '1rem', 'display':'inline'}),
                    html.Img(id='i-paral', src=i_src, style=info_style),

                    tooltips.text_tooltip(tooltips.resolutionMS2, 'MS2-resolution-header'),
                    tooltips.text_tooltip(tooltips.resolutionMS2, 'i-ms2-resolution'),
                    tooltips.text_tooltip(tooltips.MaxIT, 'IT-MS2-header'),
                    tooltips.text_tooltip(tooltips.MaxIT, 'i-ms2-mit'),
                    tooltips.text_tooltip(tooltips.topN, 'topN-topSpeed-choice'),
                    tooltips.text_tooltip(tooltips.topN, 'i-topN'),
                    tooltips.text_tooltip(tooltips.parallel, 'i-paral'),
                    ],
                    style=block_style)

def res_plot_html():
    '''
    Mass Spectral Resolution plot
    '''
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
    '''
    Cycle time plot
    '''
    return html.Div([
                    html.Center([
                            html.H6('Cycle Time', id='cycle-time-header'),
                            html.Img(id='i-cycle-time', src=i_src, style=info_style),
                            dcc.Graph(id='cycle-time-graph')
                            ]),
                        
                            tooltips.text_tooltip(tooltips.cycle_time, 'cycle-time-header'),
                            tooltips.text_tooltip(tooltips.cycle_time, 'i-cycle-time')
                    ],
                    style=cycle_figure_style)

### End interface building blocks ###

### Main window ###
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
        tooltips.logo_tooltip()
             ], style={'display': 'flex', 'padding-bottom': '1rem'}),
    
    #upper part - info table, dynamic range plot, observed peptides
    table_dynRange_html(),
    
    #simulated mass spectrum and points-pep-peak graph
    html.Div([
            dcc.Graph(id='main-graph', style=main_graph_style),
            html.Div([
                    html.H6('Peptide Elution Profile', id='ppp-header'),
                    html.Img(id='i-ppp-graph', src=i_src, style=info_style),
                    dcc.Graph(id='ppp-graph', config={'displayModeBar': False}),
                    
                    tooltips.text_tooltip(tooltips.ppp, 'i-ppp-graph'),
                    tooltips.text_tooltip(tooltips.ppp, 'ppp-header')
                    ], style=ppp_figure_style)
           
            ], style=big_panel_style),

    #model parameters switches
    html.Div([
            #Block distribution, TIC and Acquisition
            block_global_html(),
            #Block MS1 parameters
            block_MS1_html(),
            #Block MS2 parameters
            block_MS2_html(),
            ], style=big_panel_style),

    #lower part - resolution plot, cycle time plot
    html.Div([
            res_plot_html(),
            cycle_time_html()
            ], style=big_panel_style),

    #footer part
    html.Div([
        html.Img(src='/assets/sdu_logo.png',
                 style={'height': '30px'}),
             ], style={'textAlign': 'center'}),
                 
    html.Div([
        html.P('Department of Biochemistry and Molecular Biology, University of Southern Denmark'),
        html.P(['Do you have any questions and suggestions about HUMOS? Contact us via ',
			    html.A('Github', href='https://github.com/SimpleNumber/HUMOS'),
			   u' write to vgor (\u0430t) bmb.sdu.dk or juliabubis (\u0430t) gmail.com'
               ])
            ], style={'textAlign': 'center'}),
    
        ], style={'margin': '25px'}
        
    ) #end layout
### End main window ###

### Callback functions ###
def update_figure(selected_resolution, selected_agc, distribution, mit_clicked,
                  method, ionFlux, relayout_data, max_it):
    '''
    Update of the main graph, dynamic range graph and information table
    '''
    
    boxCar = (method == 'bc')
    resolution = params.resolutions_list[selected_resolution]
    agc = params.agc_list[selected_agc]
    ionFlux = params.TIC[ionFlux]

    mechanics.scale_ion_currents(ion_data, ionFlux) #apply TIC value
    
    centroid_spectrum, real_st, real_agc, peptides, max_int, min_int = \
               mechanics.get_full_spectrum(ion_data, distribution, agc, max_it)
    
    #DataFrame with dynamic range information
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

    #save zoom region or use default
    if relayout_data == None:
        x_range = [min(main_spectrum[0]) - 0.5, max(main_spectrum[0]) + 0.5]
        y_range = [0, max(main_spectrum[1]) * 1.01]
    else:
        x_range, y_range = art.get_zoom(relayout_data,
                                    min(main_spectrum[0]) - 0.5,
                                    max(main_spectrum[0]) + 0.5,
                                    0,
                                    max(main_spectrum[1]) * 1.01)
    
    main_traces = [go.Scatter(x=main_spectrum[0],
                              y=main_spectrum[1],
                              name='MS1 spectrum')]
    
    #process BoxCar spectra
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
    
    #information table
    table = art.make_table(real_sts, real_agcs, ['MS1'] + labels_bc, resolution)
    
    #finalize dynamic range DataFrame
    dr_df['y'] = [[i, i] for i in range(0, len(dr_df))]
    dr_df.index = dr_df['text']
    
    return [dbc.Table.from_dataframe(table),
            {'data': main_traces, 'layout': art.get_main_layout(x_range, y_range)},
            {'data': dr_df.apply(art.get_dynrange_trace, axis=1).tolist(),
             'layout': art.get_dynrange_layout(dr_df)},        
            {'data': art.get_obsPep_trace(observed_peptides, colors[0], colors[1]),
             'layout': art.get_obsPep_layout()}]
            
def update_ms_counts(topN, method, data, selected_resolution, ms2_resolution, 
                     parallel, mit_clicked, main_graph, mit_ms2, top_mode):
    '''
    Update counts of MS spectra, cycle time graph and ponts-per-peak plot
    '''
    
    boxCar = (method == 'bc')
    parallel = True if len(parallel) > 0 else False #value is a list of selected options
    ms2_resolution = params.resolutions_list[ms2_resolution]
    resolution = params.resolutions_list[selected_resolution]
    
    if data == None:
       return None #void return, before table data is ready

    #parse infromation table
    data = art.tabletodf(data) 
    data = data.iloc[:, 1:].apply(pd.to_numeric)
    
    #cycletime calculation paramters
    ccParam = {'resolution' : resolution,
               'ms2resolution': ms2_resolution,
               'ms2IT': mit_ms2,
               'LC_time': params.LC_time, 
               'parallel': parallel}
    
    #translate topN and topSpeed
    if top_mode == 'topN':
        ccParam['topN'] = topN
    elif top_mode == 'topSpeed':
        ccParam['topSpeed']= topN * 1000

    #translate scan method
    if boxCar:
        ccParam['scan_method'] = 'boxcar'
        ccParam['acc_time'] = data.iloc[0,:]
    else:
        ccParam['scan_method'] = 'full'
        ccParam['acc_time'] = data.iloc[0,0]

    #perform calculation
    cycletime, topN, ms1_scan_n, ms2_scan_n, queues = mechanics.get_MS_counts(**ccParam)
    
    ms1_scan_text = 'MS1 Scans in {} minutes: {}'.format(params.LC_time, ms1_scan_n)
    ms2_scan_text = 'MS2 Scans in {} minutes: {}'.format(params.LC_time, ms2_scan_n)
    
    #preparing data for cycle plot
    #Add all data for accumulation traces;
    #labels are hovering above blocks, names are used in the legend
    ia_labels = ['Accumulation ' + i for i in list(data.columns) + ['MS2 ' + str(j) for j in range(1, topN + 1)]]
    ia_names = ['Accumulation ' + i for i in list(data.columns) + ['MS2'] * topN]
    
    #Add MS1 and BoxCar scans to OT traces
    ot_labels = ['Acquisition ' + i for i in list(data.columns)] 
    ot_names = ot_labels[:]
    
    #Converting linear coordinates to circular
    theta_start = queues['IS'][:, 0] / cycletime * 360
    theta_end = queues['IS'][:, 1] / cycletime * 360
    theta_start = np.concatenate((theta_start, queues['OT'][:, 0] / cycletime * 360))
    theta_end = np.concatenate((theta_end, queues['OT'][:, 1] / cycletime * 360))
    
    if queues['IT'].shape[0] > 0: #IT was used
        it_labels = ['Acquisition MS2 ' + str(i) for i in range(1, topN + 1)]
        it_names = ['Acquisition MS2'] * topN
        theta_start = np.concatenate((theta_start, queues['IT'][:, 0] / cycletime * 360))
        theta_end = np.concatenate((theta_end, queues['IT'][:, 1] / cycletime * 360))
    else:
        it_labels = []
        it_names = []
        ot_labels += ['Acquisition MS2 ' + str(i) for i in range(1, topN + 1)]
        ot_names += ['Acquisition MS2'] * topN
    
    #colors of traces
    main_colors = colors[2: 2 + len(data.columns)] + [qualitative.Dark2[-1]] * topN
    
    #select information to be shown in the legend
    #show names for MS1 and BoxCar (data.columns) and one label for MS2 (if there any MS2)
    #repeat twice, first for accumulation traces, second for acquisition traces
    if topN > 0:
        show_legend = ([True] * (len(data.columns) + 1) + [False] * (topN - 1)) * 2
    else:
        show_legend = ([True] * (len(data.columns))) * 2

    #create DataFrame with all information
    cycle_df = pd.DataFrame({'text': ia_labels + ot_labels + it_labels,
                             'mode': 'lines',
                             'name': ia_names + ot_names + it_names,
                             'hoverinfo': 'text',
                             'hoveron': 'fills',
                             'showlegend': show_legend,
                             'r': [0.9] * len(ia_labels) + [0.7] * len(ot_labels) + [0.5] * len(it_labels),
                             'line_color': [art.lightening_color(i) for i in main_colors] + main_colors,
                             'line_width': 13,
                             'start': theta_start,
                             'end': theta_end,
                            })

    cycle_df['theta'] = cycle_df.loc[:, ['start', 'end']].apply(art.get_range, axis=1)
    
    #collecting traces
    cycle_traces = art.get_cycle_grid()
    cycle_traces += cycle_df.apply(art.get_cycle_trace, axis=1).tolist()
    cycle_traces.append(art.get_cycle_texts(cycletime, topN, ms1_scan_text, ms2_scan_text))
    
    ##Points per peak plot
    ppp_data = [] #default data is empty
    
    #Detect highest peak in all spectra
    maxMz = 0
    maxInt = 0
    for trace in main_graph['data']:
        maxI = np.argmax(trace['y'])
        if trace['y'][maxI] > maxInt:
            maxInt = trace['y'][maxI]
            maxMz = trace['x'][maxI]
    
    if maxInt > 0: #non-empty spectrum
        topPeptide = ion_data.loc[(ion_data['mz'] - maxMz).abs().idxmin(), :]
        
        #parameters of LC peak
        center = 3
        width = 2
        top = 1
        
        #theoretical LC peak
        tRT = np.linspace(0, 10, 100)
        tProfile = mechanics.get_LC_profile(center, top, width, tRT)
        
        #sampling LC peak
        sRT = np.arange(-cycletime/2000, 10, cycletime/1000)
        sRT = np.append(sRT, sRT[-1] + cycletime/1000) #last point
        sProfile = mechanics.get_LC_profile(center, top, width, sRT)
    
        ppp_data = art.get_ppp_trace(tRT, tProfile, colors[1],#theoretical
                                     sRT, sProfile, colors[0],#sampling
                                     topPeptide)
    
    return  [{'data': cycle_traces, 'layout': art.get_cycle_layout()},
             {'data': ppp_data, 'layout': art.get_ppp_layout()}]
                      
def update_resolution_graph(selected_resolution):
    '''
    Update resolution graph
    '''
    #ion trap has resolution ~1000 (0.1 Th) at 100
    resolution = 1000 if selected_resolution == 0 else params.resolutions_list[selected_resolution]
    
    resolution_spectrum = mechanics.get_profile_spectrum(tmt_spectrum, resolution, points=51)
    resolution_traces = [go.Scatter(x=resolution_spectrum[0],
                                    y=resolution_spectrum[1],
                                    name=' '.join(['R = ', str(resolution)])),
        
                         go.Scatter(x=[tmt_spectrum[0,0], tmt_spectrum[0,0]],
                                    y=[0, tmt_spectrum[0,1]],
                                    text='TMT 127N',
                                    name='',
                                    mode='lines',
                                    line={'color': qualitative.Dark2[5]}),
                                    
                         go.Scatter(x=[tmt_spectrum[1,0], tmt_spectrum[1,0]],
                                    y=[0, tmt_spectrum[1,1]],
                                    text='TMT 127C',
                                    name='',
                                    mode='lines',
                                    line={'color': qualitative.Dark2[6]}
                                    ) ]

    return [ {'data': resolution_traces,
              'layout': go.Layout(margin={'t': 10,
                                          'l': 50},
                                  showlegend=False,
                                  xaxis={'title': 'm/z'},
                                  yaxis={'title': 'Abundance'})
              } ]

def update_top_slider(choice_value):
    '''
    Switch between TopN and TopSpeed sliders
    '''
    if choice_value == 'topN':
        return (0, #min
                40, #max
                1, #step
                15, #value
                {i: '{}'.format(i) for i in range(0, 41, 5)}) #marks
    elif choice_value == 'topSpeed':
        return (0.0, #min
                5.0, #max
                0.05, #step
                2.0, #value
                {i: '{}s'.format(i) for i in range(0, 6)}) #marks
    else:
        raise ValueError('Unknown value ({}) in TopN-TopSpeed choice'.format(choice_value))
### End callback functions ###

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
    [Output('cycle-time-graph', 'figure'),
     Output('ppp-graph', 'figure')],
    [Input('topN-slider', 'value'),
     Input('method-choice', 'value'),
     Input('table','children'),
     Input('resolution-slider', 'value'),
     Input('resolution-ms2-slider', 'value'),
     Input('paral-checklist', 'value'),
     Input('it-ms2-button','n_clicks')],
    [State('main-graph', 'figure'),
     State('mit-ms2-box', 'value'),
     State('topN-topSpeed-choice', 'value')])(update_ms_counts)

app.callback(
    [Output('resolution-graph', 'figure')],
    [Input('resolution-ms2-slider', 'value')])(update_resolution_graph)

app.callback(
        [Output('topN-slider', 'min'),
         Output('topN-slider', 'max'),
         Output('topN-slider', 'step'),
         Output('topN-slider', 'value'),
         Output('topN-slider', 'marks')],
        [Input('topN-topSpeed-choice', 'value')])(update_top_slider)

server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)
