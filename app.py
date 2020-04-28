
import urllib.parse
import base64
import io

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

import pandas as pd

from preprocess import preprocess, calculate_pairwise_similarities
from constants import functional_groups, reaction_classes

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.DARKLY],
    meta_tags=[
        {
            'name': 'viewport',
            'content': 'width=device-width, initial-scale=1, maximum-scale=1.0, user-scalable=no',
        }
    ],
)

server = app.server

func_group_options = [{'label': group, 'value': value} for group, value in functional_groups.items()]
reaction_class_options = [{'label': group, 'value': value} for group, value in reaction_classes.items()]

header = html.Div(
    [
        html.H2('molecular DB builder.'),
        html.Div(id="molecule-count"),
    ],
    className='header'
)

reaction_class_selection = html.Div(
    className='group-selector',
    children=[
        html.Label('Select Reaction Classes', id='rxn-selector-descr'),
        html.Div(
            children=dcc.Dropdown(
                className='dropdown',
                id='rxn-class-select',
                options=reaction_class_options,
                value=[],
                multi=True,
                searchable=True,
                placeholder='Select classes...'
            ),
        ),
        dbc.Tooltip(
            'Molecules compatible with these reaction classes will be added.',
            target='rxn-selector-descr',
        ),
    ],
)

func_group_selection = html.Div(
    className='group-selector',
    children=[
        html.Label('Select Functional Groups', id='fgroup-selector-descr'),
        html.Div(
            children=dcc.Dropdown(
                className='dropdown',
                id='fgroup-class-select',
                options=func_group_options,
                value=[],
                multi=True,
                searchable=True,
                placeholder='Select groups...'
            ),
        ),
        dbc.Tooltip(
            'Molecules containing these functional groups will be added.',
            target='fgroup-selector-descr',
        ),
    ],
)

func_group_exclusion = html.Div(
    className='group-selector',
    children=[
        html.Label('Exclude Functional Groups', id='fgroup-exclude-descr'),
        html.Div(
            children=dcc.Dropdown(
                className='dropdown',
                id='fgroup-class-exclude',
                options=func_group_options,
                value=[],
                multi=True,
                searchable=True,
                placeholder='Select groups...'
            ),
        ),
        dbc.Tooltip(
            'Molecules containing these functional groups will be excluded.',
            target='fgroup-exclude-descr',
        ),
    ],
)

slider_molwt = html.Div(
    [
        html.Div(id='slider-molwt-output', className='slider-annotation'),
        dcc.Slider(
            id='slider-molwt',
            min=0,
            max=1000,
            step=0.1,
            value=1000,
        ),
    ]
)

slider_logp = html.Div(
    [
        html.Div(id='slider-logp-output', className='slider-annotation'),
        dcc.Slider(
            id='slider-logp',
            min=0,
            max=20,
            step=0.1,
            value=20,
        ),
    ]
)

upload_button = dcc.Upload(
    dbc.Button(
        [
            html.Img(src='./assets/upload.svg', className='button-img'),
            'Import'
        ],
        size='lg', className='button-import'),
    id='upload-data',
    multiple=True,
)


export_button = dbc.Button(
    [
        html.Img(src='./assets/download.svg', className='button-img'),
        html.A('Export', id='download-link', download='rawdata.csv', href='', target='_blank'),
    ],
    size='lg', outline=True, className='button-export'
)

explainer = html.Div(
    [
        html.Hr(),
        html.Div(
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H4('1. Import'),
                            html.Img(src='./assets/upload-blk1.svg', className='explainer-img'),
                            html.P('a .csv of SMILES strings to be filtered'),
                        ],
                        lg=4,

                    ),
                    dbc.Col(
                        [
                            html.H4('2. Filter'),
                            html.Img(src='./assets/filter.svg', className='explainer-img'),
                            html.P('to your specifications using the selection tools'),
                        ],
                        lg=4,

                    ),
                    dbc.Col(
                        [
                            html.H4('3. Export'),
                            html.Img(src='./assets/download.svg', className='explainer-img'),
                            html.P('the filtered database'),
                        ],
                        lg=4,

                    ),
                ]
            ),
            className='explainer-section',
        )
    ]
)

##########
# LAYOUT #
##########

app.layout = html.Div(
    [
        header,
        dbc.Row(
            [
                dbc.Col(
                    [
                        reaction_class_selection,
                        func_group_selection,
                        func_group_exclusion,

                        html.Div(
                            [
                                slider_molwt,
                                slider_logp
                            ],
                            className='slider-container'
                        ),
                        html.Div(
                            [
                                upload_button,
                                export_button
                            ],
                            className='button-container'
                        )
                    ],
                    lg=4,
                    className='column-left'
                ),
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.Div([html.H6('LogP')], className='plot-header', id='plot-header-1'),
                                dbc.Tooltip(
                                    'Calculated LogP distribution of DB building blocks.',
                                    target='plot-header-1'),
                                dcc.Loading(
                                    type='default',
                                    children=[dcc.Graph(id='updating-graph1')],
                                    className='spinner'
                                ),
                            ],
                            className='plot-container'),
                        html.Div(
                            [
                                html.Div([html.H6('Molecular Weight')], className='plot-header', id='plot-header-2'),
                                dbc.Tooltip(
                                    'Calculated molwt distribution of DB building blocks.',
                                    target='plot-header-2'
                                ),
                                dcc.Loading(
                                    type='default',
                                    children=[dcc.Graph(id='updating-graph2')],
                                ),
                            ],
                            className='plot-container'),
                    ],
                    lg=4,
                    className='column'
                ),
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.Div([html.H6('Pairwise Similarity', className='plot-header-text')],
                                         className='plot-header', id='plot-header-3'),
                                dbc.Tooltip(
                                    '''Calculated pairwise Tanimoto similarity distribution of DB building blocks
                                    (Computed only for filtered dataset).''',
                                    target='plot-header-3'),
                                dcc.Loading(
                                    type='default',
                                    children=[dcc.Graph(id='updating-graph3')],
                                    className='spinner'
                                ),
                            ],
                            className='plot-container'),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H6('Functional Groups', className='plot-header-text')
                                    ],
                                    className='plot-header', id='plot-header-4'),
                                dbc.Tooltip('Count of functional groups in DB building blocks.',
                                            target='plot-header-4'),
                                dcc.Loading(
                                    type='default',
                                    children=[dcc.Graph(id='updating-graph4')],
                                    className='spinner'
                                ),
                            ],
                            className='plot-container'),
                    ],
                    lg=4,
                    className='column'
                )
            ],
            className='column-container'
        ),
        html.Div(id='export-data', style={'display': 'none'}),
        explainer
    ],
    className='master-div'
)


###########
# FIGURES #
###########

def generate_histgram_content(x, x_filtered, title):

    data_array = [
        {
            'x': x_filtered,
            'type': 'histogram',
            'marker': {'color': '#9656a1'},
            'name': 'filtered'
        },
    ]
    if x is not None:
        data_array = [
            {
                'x': x,
                'type': 'histogram',
                'opacity': 1.0,
                'marker': {'color': '#9e9e9e'},
                'name': 'original'
            },
        ] + data_array

    return {
        'data': data_array,
        'layout': {
            'autosize': True,
            'automargin': True,
            'showlegend': True,
            'paper_bgcolor': 'rgb(0, 0, 0, 0)',
            'plot_bgcolor': '#f7f7f7',
            'barmode': 'overlay',
            'bargap': '0.01',
            'font': {'color': 'grey'},
            'height': '230',
            'margin': dict(l=50, r=50, b=50, t=50, pad=0),
            'yaxis': {'gridcolor': '#dedede'},
            'legend': {
                'x': 0.8,
                'y': 1.0,
                'bgcolor': '#f7f7f7',
                'borderwidth': 0
            }
        }
    }


def generate_bargraph_content(x, x_filtered, title):

    y_data_filtered = [x_filtered[fg].sum() for fg in list(functional_groups.keys())]
    y_data = [x[fg].sum() for fg in list(functional_groups.keys())]
    fgs = list(functional_groups.keys())

    data_array = [
        {
            'x': fgs,
            'y': y_data_filtered,
            'type': 'bar',
            'marker': {'color': '#9656a1'},
            'name': 'filtered'
        },
    ]
    if x is not None:
        data_array = [
            {
                'x': fgs,
                'y': y_data,
                'type': 'bar',
                'opacity': 1.0,
                'marker': {'color': '#9e9e9e'},
                'name': 'original'
            },
        ] + data_array

    return {
        'data': data_array,
        'layout': {
            'autosize': True,
            'automargin': True,
            'showlegend': False,
            'paper_bgcolor': 'rgb(0, 0, 0, 0)',
            'plot_bgcolor': '#f7f7f7',
            'barmode': 'overlay',
            # 'bargap': '0.03',
            'font': {'color': 'grey'},
            'height': '230',
            'margin': dict(l=50, r=50, b=50, t=50, pad=0),
            'yaxis': {'gridcolor': '#dedede'},
            'legend': {
                'x': 0.8,
                'y': 1.0,
                'bgcolor': '#f7f7f7',
                'borderwidth': 0
            }
        }
    }


###########
# SLIDERS #
###########

'''MolWt Slider'''
@app.callback(
    dash.dependencies.Output('slider-molwt-output', 'children'),
    [dash.dependencies.Input('slider-molwt', 'value')])
def update_molwt_output(value):
    return 'MolWt Cutoff: {:.2f}'.format(value)


'''LogP Slider'''
@app.callback(
    dash.dependencies.Output('slider-logp-output', 'children'),
    [dash.dependencies.Input('slider-logp', 'value')])
def update_logp_output(value):
    return 'LogP Cutoff: {:.2f}'.format(value)


#########################
# UPLOADING & FILTERING #
#########################

@app.callback(
    [
        Output('updating-graph1', 'figure'),
        Output('updating-graph2', 'figure'),
        Output('updating-graph3', 'figure'),
        Output('updating-graph4', 'figure'),
        Output('export-data', 'children'),
        Output('download-link', 'href'),
        Output('molecule-count', 'children')
    ],
    [
        Input('upload-data', 'contents'),
        Input('rxn-class-select', 'value'),
        Input('fgroup-class-select', 'value'),
        Input('fgroup-class-exclude', 'value'),
        Input('slider-molwt', 'value'),
        Input('slider-logp', 'value'),
    ],
    [
        State('upload-data', 'filename')
    ])
def update_output(file_contents, active_rxns, active_fgroups, inactive_fgroups, molwt_cutoff, logp_cutoff, file_name):

    if file_contents is not None:

        file_name = file_name[0]
        file_contents = file_contents[0]
        _, content_string = file_contents.split(',')
        decoded = base64.b64decode(content_string)

        try:
            if 'csv' in file_name:
                # Assume that the user uploaded a CSV file
                df = pd.read_csv(
                    io.StringIO(decoded.decode('utf-8')))
            elif 'xls' in file_name:
                # Assume that the user uploaded an excel file
                df = pd.read_excel(io.BytesIO(decoded))
        except Exception as e:
            print(e)
            return html.Div([
                'There was an error processing this file - did you upload a .csv or excel file?'
            ])

        # preprocess uploaded data
        data = preprocess(df)

    else:
        data = preprocess(pd.DataFrame({'smiles': []}))

    # include reaction classes and functional groups
    filtered_data = pd.DataFrame(columns=data.columns)
    for group in (active_fgroups + active_rxns):
        filtered_data = pd.concat([filtered_data, data[data[group] > 0]])

    # exclude functional groups
    filtered_data = filtered_data.drop_duplicates()
    for group in inactive_fgroups:
        filtered_data = filtered_data[filtered_data[group] == 0]

    # filter by molwt
    filtered_data = filtered_data[filtered_data.molwt < molwt_cutoff]

    # filter by logp
    filtered_data = filtered_data[filtered_data.logp < logp_cutoff]

    # update export link
    export_data = filtered_data['smiles']
    csv_string = export_data.to_csv(index=False, encoding='utf-8')
    csv_string = 'data:text/csv;charset=utf-8,' + urllib.parse.quote(csv_string)

    # calculate pairwise similarities of filtered data only
    similarities = calculate_pairwise_similarities(filtered_data)

    # retrieve figure contents
    logp_figure = generate_histgram_content(data['logp'], filtered_data['logp'], 'LogP')
    molwt_figure = generate_histgram_content(data['molwt'], filtered_data['molwt'], 'MolWt')
    similarity_figure = generate_histgram_content(None, similarities, 'Pairwise Similarity')
    fg_figure = generate_bargraph_content(data, filtered_data, 'Functional Groups')

    return [
        logp_figure,
        molwt_figure,
        similarity_figure,
        fg_figure,
        filtered_data.to_json(date_format='iso', orient='split'),
        csv_string,
        f'molecule count: {len(filtered_data.smiles)}',
    ]


if __name__ == '__main__':
    app.run_server(debug=True)
