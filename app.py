import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["NUMEXPR_MAX_THREADS"] = "8"

import dash
from dash import html, dcc, Input, Output, State, callback
import plotly.graph_objects as go
import networkx as nx
from rcbgk import *
import numpy as np
import json
from graphviz import Digraph
import base64
from io import BytesIO
import random

def create_image_figure(image_data, title=""):
    """Create a figure with proper aspect ratio"""
    return {
        'data': [],
        'layout': {
            'images': [{
                'source': image_data,
                'xref': "paper",
                'yref': "paper",
                'x': 0,
                'y': 1,
                'sizex': 1,
                'sizey': 1,
                'sizing': "contain",  # Changed from stretch to contain
                'layer': "below"
            }],
            'xaxis': {'visible': False, 'range': [0, 1]},
            'yaxis': {'visible': False, 'range': [0, 1], 'scaleanchor': 'x'},
            'width': 800,
            'height': 600,
            'margin': {'l': 0, 'r': 0, 't': 30, 'b': 0},
            'title': title
        }
    }

# Initialize the Dash app
app = dash.Dash(__name__, suppress_callback_exceptions=True)

# Define styles
CONTENT_STYLE = {
    "margin-left": "2rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}

INPUT_STYLE = {
    "width": "200px",
    "margin": "10px",
    "display": "inline-block"
}

BUTTON_STYLE = {
    "margin": "10px",
    "backgroundColor": "#119DFF",
    "color": "white"
}

# App layout
app.layout = html.Div([
    # Header
    html.H1("Causal Discovery and Analysis Tool", 
            style={"textAlign": "center", "color": "#2c3e50", "marginBottom": "30px"}),
    
    # Main content
    html.Div([
        dcc.Tabs([
            # Tab 1: Graph Generation & Conversion
            dcc.Tab(label='Graph Generation & Conversion', children=[
                html.Div([
                    # Random Chordal Graph Generation
                    html.Div([
                        html.H3("Generate Random Chordal Graph"),
                        html.Div([
                            html.Label("Number of Vertices:"),
                            dcc.Input(
                                id='chordal-vertices',
                                type='number',
                                value=10,
                                min=3,
                                max=50,
                                style=INPUT_STYLE
                            ),
                            html.Label("Number of Edges:"),
                            dcc.Input(
                                id='chordal-edges',
                                type='number',
                                value=15,
                                min=3,
                                style=INPUT_STYLE
                            ),
                            html.Button(
                                'Generate Chordal Graph',
                                id='generate-chordal-btn',
                                n_clicks=0,
                                style=BUTTON_STYLE
                            ),
                        ]),
                        dcc.Loading(
                            id="loading-chordal",
                            children=[
                                html.Div(id='chordal-results'),
                                dcc.Graph(id='chordal-graph')
                            ]
                        )
                    ], style={'marginBottom': '30px'}),

                    # Random DAG Generation
                    html.Div([
                        html.H3("Generate Random DAG"),
                        html.Div([
                            html.Label("Number of Vertices:"),
                            dcc.Input(
                                id='dag-vertices',
                                type='number',
                                value=10,
                                min=3,
                                max=50,
                                style=INPUT_STYLE
                            ),
                            html.Label("Edge Probability:"),
                            dcc.Input(
                                id='edge-prob',
                                type='number',
                                value=0.3,
                                min=0,
                                max=1,
                                step=0.1,
                                style=INPUT_STYLE
                            ),
                            html.Button(
                                'Generate DAG',
                                id='generate-dag-btn',
                                n_clicks=0,
                                style=BUTTON_STYLE
                            ),
                        ]),
                        dcc.Loading(
                            id="loading-dag",
                            children=[
                                html.Div(id='dag-results'),
                                dcc.Graph(id='dag-graph')
                            ]
                        )
                    ], style={'marginBottom': '30px'}),

                    # Graph Conversion
                    html.Div([
                        html.H3("Graph Conversion"),
                        html.Div([
                            html.Label("Input Graph:"),
                            dcc.Textarea(
                                id='convert-input',
                                placeholder="Enter graph as dictionary",
                                value="{0:{1,2}, 1:{0,2}, 2:{0,1}}",
                                style={'width': '100%', 'height': 100}
                            ),
                            html.Label("Conversion Type:"),
                            dcc.Dropdown(
                                id='conversion-type',
                                options=[
                                    {'label': 'DAG to CPDAG', 'value': 'dag2cpdag'},
                                ],
                                value='dag2cpdag'
                            ),
                            html.Button(
                                'Convert/Analyze',
                                id='convert-analyze-btn',
                                n_clicks=0,
                                style=BUTTON_STYLE
                            ),
                        ]),
                        # Additional input for chordal comparison
                        html.Div(
                            id='chordal-compare-container',
                            children=[
                                html.Label("Second Graph:"),
                                dcc.Textarea(
                                    id='second-graph-input',
                                    placeholder="Enter second graph for comparison",
                                    value="{0:{1,2}, 1:{0,2}, 2:{0,1}}",
                                    style={'width': '100%', 'height': 100}
                                ),
                            ],
                            style={'display': 'none'}
                        ),
                        dcc.Loading(
                            id="loading-convert",
                            children=[
                                html.Div(id='conversion-results'),
                                dcc.Graph(id='conversion-graph-input'),
                                dcc.Graph(id='conversion-graph-output')
                            ]
                        )
                    ])
                ], style=CONTENT_STYLE)
            ]),

            # Tab 2: Background Knowledge Analysis
            dcc.Tab(label='Background Knowledge Analysis', children=[
                html.Div([
                    # Background Knowledge Generation
                    html.Div([
                        html.H3("Generate Background Knowledge"),
                        html.Div([
                            html.Label("Input CPDAG or MPDAG:"),
                            dcc.Textarea(
                                id='bgk-cpdag-input',
                                placeholder="Enter CPDAG or MPDAG as dictionary (e.g., {0:{1,2}, 1:{0,2}, 2:{0,1}})",
                                value="{0:{1,2}, 1:{0,2}, 2:{0,1}}",
                                style={'width': '100%', 'height': 100}
                            ),
                            html.Label("Number of Constraints:"),
                            dcc.Input(
                                id='num-constraints',
                                type='number',
                                value=3,
                                min=1,
                                style=INPUT_STYLE
                            ),
                            html.Label("Background Knowledge Type:"),
                            dcc.Dropdown(
                                id='bgk-type',
                                options=[
                                    {'label': 'Direct Causal Relations', 'value': 'direct'},
                                    {'label': 'Ancestral Relations (Including Indirect)', 'value': 'ancestral'},
                                    {'label': 'Non-causal Relations', 'value': 'non-ancestral'}
                                ],
                                value='direct'
                            ),
                            html.Button(
                                'Generate BGK',
                                id='generate-bgk-btn',
                                n_clicks=0,
                                style=BUTTON_STYLE
                            ),
                        ]),
                        dcc.Loading(
                            id="loading-bgk",
                            children=[
                                html.Div(id='bgk-results'),
                                dcc.Graph(id='bgk-graph')
                            ]
                        )
                    ], style={'marginBottom': '30px'}),

                    # CPDAG to MPDAG Conversion
                    html.Div([
                        html.H3("CPDAG or MPDAG plus Knowledge to MPDAG Conversion"),
                        html.Div([
                            html.Label("Input CPDAG or MPDAG:"),
                            dcc.Textarea(
                                id='mpdag-input',
                                placeholder="Enter CPDAG or MPDAG as dictionary",
                                value='{"0":{"1","2"}, "1":{"0","2"}, "2":{"0","1"}}',
                                style={'width': '100%', 'height': 100}
                            ),
                            html.H4("Background Knowledge"),
                            html.P([
                                "Background knowledge format: [(x,y,t)], where:",
                                html.Br(),
                                "• t=1: x is a cause of y (x → y)",
                                html.Br(),
                                "• t=0: x is NOT a cause of y (x ↛ y)"
                            ]),
                            dcc.Textarea(
                                id='mpdag-bgk-input',
                                placeholder="Enter background knowledge",
                                value="""[
    ('0', '1', 1),  # 0 is a cause of 1
    ('1', '2', 1),  # 1 is a cause of 2
    ('1', '0', 0),  # 1 is NOT a cause of 0
    ('2', '0', 0)   # 2 is NOT a cause of 0
]""",
                                style={'width': '100%', 'height': 150}
                            ),
                            html.Button(
                                'Convert to MPDAG',
                                id='mpdag-convert-btn',
                                n_clicks=0,
                                style=BUTTON_STYLE
                            ),
                        ]),
                        dcc.Loading(
                            id="loading-mpdag",
                            children=[
                                html.Div(id='mpdag-results'),
                                dcc.Graph(id='mpdag-graph')
                            ]
                        )
                    ])
                ], style=CONTENT_STYLE)
            ]),

            # Tab 3: Causal Inference
            dcc.Tab(label='Causal Inference', children=[
                html.Div([
                    # Step 1: Input/Generate DAG with weights
                    html.Div([
                        html.H3("Step 1: Input DAG with Weights"),
                        html.P([
                            "First, input a DAG with edge weights learned from observational data or generate a random one for testing. ",
                            "This DAG represents a possible causal structure in the Markov equivalence class of the true graph. ",
                            "The weights could be learned from observational data."
                        ], style={"marginBottom": "15px"}),
                        
                        # DAG input
                        html.Label("Input DAG:"),
                        dcc.Textarea(
                            id='causal-dag-input',
                            placeholder="Enter DAG as dictionary (e.g., {'0':{'1','2'}, '1':{'2'}, '2':{}})",
                            value='{"0": ["1", "2", "4"], "1": ["2", "3"], "2": ["3"], "4": [], "3": []}',
                            style={'width': '100%', 'height': 100}
                        ),
                        
                        # Weights input
                        html.Label("Edge Weights:"),
                        dcc.Textarea(
                            id='causal-weights-input',
                            placeholder='Enter weights as dictionary (e.g., {"0,1": 1.5, "1,2": 0.8})',
                            value='{"0,1": 1.54, "0,2": 1.6, "0,4": 1.52, "1,2": 1.77, "1,3": 1.75, "2,3": 1.49}',
                            style={'width': '100%', 'height': 100}
                        ),
                        
                        # 添加更新按钮
                        html.Button(
                            'Update DAG and CPDAG',
                            id='update-dag-btn',
                            n_clicks=0,
                            style=BUTTON_STYLE
                        ),
                        
                        # Random generation parameters
                        html.Div([
                            html.Label("Number of Nodes:"),
                            dcc.Input(
                                id='random-dag-nodes',
                                type='number',
                                value=5,
                                min=2,
                                style=INPUT_STYLE
                            ),
                            html.Label("Edge Probability:"),
                            dcc.Input(
                                id='random-dag-edge-prob',
                                type='number',
                                value=0.3,
                                min=0.1,
                                max=0.9,
                                step=0.1,
                                style=INPUT_STYLE
                            ),
                        ]),
                        
                        # Random generation button
                        html.Button(
                            'Generate Random DAG with Weights',
                            id='generate-causal-dag-btn',
                            n_clicks=0,
                            style=BUTTON_STYLE
                        ),
                        
                        # Visualization
                        dcc.Loading(
                            id="loading-causal-dag",
                            children=[
                                html.Div(id='causal-dag-results'),
                                dcc.Graph(id='causal-dag-graph'),
                                dcc.Store(id='mcov-store'),
                                dcc.Store(id='dag-store'),
                                dcc.Store(id='cpdag-store')
                            ]
                        )
                    ], style={'marginBottom': '30px'}),
                    
                    # Step 2: Background Knowledge and Causal Effect Inference
                    html.Div([
                        html.H3("Step 2: Background Knowledge and Causal Effect Inference"),
                        html.P([
                            "Input background knowledge to obtain MPDAG, then select source and target nodes for causal effect estimation."
                        ], style={"marginBottom": "15px"}),
                        
                        # Background Knowledge Input
                        html.Label("Background Knowledge:"),
                        html.P([
                            "Format: [(x,y,t)], where t=1: x is a cause of y (x → y), t=0: x is NOT a cause of y (x ↛ y)"
                        ], style={"color": "#666666", "fontSize": "0.9em"}),
                        dcc.Textarea(
                            id='causal-bgk-input',
                            placeholder="Enter background knowledge (e.g., [('0','1',1), ('1','2',0)])",
                            value="[('2','3',1)]",
                            style={'width': '100%', 'height': 100}
                        ),
                        
                        # Source and Target Selection
                        html.Div([
                            html.Label("Source Node (x):"),
                            dcc.Input(
                                id='causal-source',
                                type='text',
                                value='1',
                                style=INPUT_STYLE
                            ),
                            html.Label("Target Node (y):"),
                            dcc.Input(
                                id='causal-target',
                                type='text',
                                value='3',
                                style=INPUT_STYLE
                            ),
                        ], style={'marginTop': '15px'}),
                        
                        # Analysis Button
                        html.Button(
                            'Analyze Causal Effects',
                            id='analyze-causal-btn',
                            n_clicks=0,
                            style=BUTTON_STYLE
                        ),
                        
                        # Results Display
                        dcc.Loading(
                            id="loading-causal-analysis",
                            children=[
                                html.Div(id='causal-analysis-results'),
                                dcc.Graph(id='causal-mpdag-graph')
                            ]
                        )
                    ], style={'marginBottom': '30px'})
                ])
            ])
        ])
    ])
])



# Callback for graph size evaluation
@callback(
    [Output('evaluation-results', 'children'),
     Output('evaluation-graph', 'figure')],
    [Input('evaluate-btn', 'n_clicks')],
    [State('graph-input', 'value'),
     State('graph-type-input', 'value')]
)
def evaluate_graph_size(n_clicks, graph_input, graph_type):
    if n_clicks == 0:
        return "", {}
    
    try:
        # Parse input graph
        G = eval(graph_input)
        
        if graph_type == 'dag':
            # Convert DAG to CPDAG first
            G_nx = nx.DiGraph(G)
            cpdag = dag_to_cpdag(G_nx)
            G = {n: set(cpdag.neighbors(n)) for n in cpdag.nodes()}
        
        # Evaluate size
        size = evalsize(G)
        
        # Create visualization using plot_graph from rcbgk
        image_data = plot_graph(nx.Graph(G), 
                              title="Input Graph")
        
        # Create figure with the image
        fig = create_image_figure(image_data, "Input Graph")
        
        # Create results text
        results = html.Div([
            html.H4("Graph Analysis Results:", style={"marginBottom": "10px"}),
            html.P(f"Graph size: {size}"),
            html.P(f"Number of vertices: {len(G)}"),
            html.P(f"Number of edges: {sum(len(v) for v in G.values())//2}"),
            html.P(f"Is chordal: {nx.is_chordal(nx.Graph(G))}")
        ])
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

# Callback for IDA analysis
@callback(
    [Output('ida-results', 'children'),
     Output('ida-graph', 'figure')],
    [Input('ida-btn', 'n_clicks')],
    [State('ida-graph-input', 'value'),
     State('ida-weights-input', 'value'),
     State('ida-source', 'value'),
     State('ida-target', 'value')],
    prevent_initial_call=True
)
def perform_ida_analysis(n_clicks, graph_input, weights_input, source_node, target_node):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # 解析输入
        input_dict = eval(graph_input)
        source_node = int(source_node)
        target_node = int(target_node)
        
        # 转换为NetworkX图
        mpdag = nx.DiGraph()
        mpdag.add_nodes_from(input_dict.keys())
        for node, neighbors in input_dict.items():
            for neighbor in neighbors:
                mpdag.add_edge(node, neighbor)
        
        # 添加权重
        if weights_input:
            weights = eval(weights_input)
            for edge_str, w in weights.items():
                # 从字符串解析出节点对
                u, v = edge_str.split(',')
                if mpdag.has_edge(u, v):
                    mpdag[u][v]['weight'] = w
        else:
            # 使用默认权重
            for u, v in mpdag.edges():
                mpdag[u][v]['weight'] = 1.0
        
        # 计算协方差矩阵
        mcov = true_cov(mpdag)
        
        # 执行IDA分析
        effects = bgk_ida(source_node, target_node, mcov, mpdag)
        
        # 创建可视化
        edge_colors = {}
        for u, v in mpdag.edges():
            if u == source_node and v == target_node:
                edge_colors[(u, v)] = '#2ecc71'
        
        image_data = plot_weighted_graph(mpdag, 
                                       title="MPDAG with IDA Analysis",
                                       edge_colors=edge_colors)
        
        fig = create_image_figure(image_data, "MPDAG with IDA Analysis")
        
        # 创建结果文本
        results = html.Div([
            html.H4("IDA Analysis Results:", style={"marginBottom": "10px"}),
            html.P(f"Source node: {source_node}"),
            html.P(f"Target node: {target_node}"),
            html.P("Edge weights:"),
            html.Pre(json.dumps({f"({u},{v})": mpdag[u][v]['weight'] 
                               for u, v in mpdag.edges()}, indent=2)),
            html.P("Possible causal effects:"),
            html.Ul([html.Li(f"{effect:.4f}") for effect in effects])
        ])
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

# Callback to show/hide BGK input based on conversion type
@callback(
    Output('bgk-input-container', 'style'),
    [Input('graph-type', 'value')]
)
def toggle_bgk_input(graph_type):
    if graph_type == 'cpdag2mpdag':
        return {'display': 'block'}
    return {'display': 'none'}

# Update the DAG generation callback
@callback(
    [Output('dag-results', 'children'),
     Output('dag-graph', 'figure')],
    [Input('generate-dag-btn', 'n_clicks')],
    [State('dag-vertices', 'value'),
     State('edge-prob', 'value')],
    prevent_initial_call=True
)
def generate_dag(n_clicks, vertices, edge_prob):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # Generate random DAG
        chordal_graph = generate_undirected_chordal(vertices, edge_prob)
        dag = nx.DiGraph()
        dag.add_nodes_from(chordal_graph.nodes())
        for u, v in chordal_graph.edges():
            if u < v:
                dag.add_edge(u, v)
            else:
                dag.add_edge(v, u)
        
        # Convert to dictionary format
        dag_dict = {n: set(dag.successors(n)) for n in dag.nodes()}
        
        # Create visualization using plot_graph from rcbgk
        image_data = plot_graph(dag, title="Generated DAG")
        
        # Create figure with the image
        fig = create_image_figure(image_data, "Generated DAG")
        
        # Create results text
        results = html.Div([
            html.H4("DAG Properties:", style={"marginBottom": "10px"}),
            html.P(f"Number of vertices: {len(dag.nodes())}"),
            html.P(f"Number of edges: {len(dag.edges())}"),
            html.P(f"Is acyclic: {nx.is_directed_acyclic_graph(dag)}"),
            html.P("DAG as dictionary:"),
            html.Pre(json.dumps(dag_dict, default=list))
        ])
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

# Callback for graph conversions
@callback(
    [Output('conversion-results', 'children'),
     Output('conversion-graph-input', 'figure'),
     Output('conversion-graph-output', 'figure')],
    [Input('convert-analyze-btn', 'n_clicks')],
    [State('conversion-type', 'value'),
     State('convert-input', 'value')],
    prevent_initial_call=True
)
def convert_and_analyze_graph(n_clicks, conv_type, graph_input):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # 解析输入图
        input_dict = eval(graph_input)
        input_graph = nx.DiGraph()
        input_graph.add_nodes_from(input_dict.keys())
        for node, neighbors in input_dict.items():
            for neighbor in neighbors:
                input_graph.add_edge(node, neighbor)
        
        # 转换为CPDAG
        output_graph = dag_to_cpdag(input_graph)
        
        # 绘制DAG和CPDAG
        dag_image = plot_graph(input_graph, title="Input DAG")
        cpdag_image = plot_graph(output_graph, title="Output CPDAG")
        
        # 创建figures
        fig_in = create_image_figure(dag_image, "Input DAG")
        fig_out = create_image_figure(cpdag_image, "Output CPDAG")
        
        # 创建结果文本
        results = html.Div([
            html.H4("Conversion Results:", style={"marginBottom": "10px"}),
            html.P(f"Input graph: {len(input_graph.nodes())} vertices, {len(input_graph.edges())} edges"),
            html.P(f"Output graph: {len(output_graph.nodes())} vertices, {len(output_graph.edges())} edges"),
            html.P("Output graph as dictionary:"),
            html.Pre(json.dumps({str(n): list(output_graph.neighbors(n)) for n in output_graph.nodes()}, default=list))
        ])
        
        return results, fig_in, fig_out
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}, {}

# Callback for CPDAG to MPDAG conversion
@callback(
    [Output('mpdag-results', 'children'),
     Output('mpdag-graph', 'figure')],
    [Input('mpdag-convert-btn', 'n_clicks')],
    [State('mpdag-input', 'value'),
     State('mpdag-bgk-input', 'value')],
    prevent_initial_call=True
)
def convert_to_mpdag(n_clicks, cpdag_input, bgk_input):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # Parse input CPDAG and background knowledge
        input_dict = eval(cpdag_input)
        bgk = eval(bgk_input) if bgk_input else []
        
        # Create a mapping between original labels and numeric indices
        nodes = sorted(input_dict.keys())
        label_to_index = {node: i for i, node in enumerate(nodes)}
        index_to_label = {i: node for node, i in label_to_index.items()}
        
        # Convert dictionary to NetworkX DiGraph
        input_graph = nx.DiGraph()
        input_graph.add_nodes_from(range(len(nodes)))
        
        # Add edges maintaining CPDAG structure
        for node, neighbors in input_dict.items():
            node_idx = label_to_index[node]
            for neighbor in neighbors:
                neighbor_idx = label_to_index[neighbor]
                # Only add edge in one direction for directed edges
                if node in input_dict.get(neighbor, []):
                    # Undirected edge - add both directions
                    input_graph.add_edge(node_idx, neighbor_idx)
                    input_graph.add_edge(neighbor_idx, node_idx)
                else:
                    # Directed edge - add only one direction
                    input_graph.add_edge(node_idx, neighbor_idx)
        
        # Convert background knowledge to numeric indices
        numeric_bgk = []
        for x, y, t in bgk:
            x_idx = label_to_index[x] if isinstance(x, str) else x
            y_idx = label_to_index[y] if isinstance(y, str) else y
            numeric_bgk.append((x_idx, y_idx, t))
        
        # Convert CPDAG and background knowledge to MPDAG
        dcc = transCbgk(input_graph, numeric_bgk)
        mpdag = construct_mpdag(input_graph, dcc)
        
        # Convert back to original labels
        output_graph = nx.DiGraph()
        output_graph.add_nodes_from(nodes)
        
        # Add edges while preserving direction
        for u, v in mpdag.edges():
            u_label = index_to_label[u]
            v_label = index_to_label[v]
            output_graph.add_edge(u_label, v_label)
        
        # Create visualization
        image_data = plot_graph(output_graph, title="MPDAG with Background Knowledge")
        
        # Create figure
        fig = create_image_figure(image_data, "MPDAG with Background Knowledge")
        
        # Create results text
        results = html.Div([
            html.H4("MPDAG Conversion Results:", style={"marginBottom": "10px"}),
            html.P(f"Number of vertices: {len(output_graph.nodes())}"),
            html.P(f"Number of edges: {len(output_graph.edges())}"),
            html.P("Applied background knowledge:"),
            html.Ul([
                html.Li([
                    f"({x}, {y}): ",
                    html.Span(
                        {
                            1: "Cause",
                            0: "Non-cause"
                        }[t],
                        style={"color": {
                            1: "#2ecc71",  # 绿色表示因果关系
                            0: "#e74c3c"   # 红色表示非因果关系
                        }[t]}
                    )
                ]) for x, y, t in bgk
            ]),
            html.P("MPDAG as dictionary:"),
            html.Pre(json.dumps({str(n): list(output_graph.neighbors(n)) 
                               for n in output_graph.nodes()}, default=list))
        ])
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

# Callback for background knowledge generation
@callback(
    [Output('bgk-results', 'children'),
     Output('bgk-graph', 'figure')],
    [Input('generate-bgk-btn', 'n_clicks')],
    [State('bgk-cpdag-input', 'value'),
     State('num-constraints', 'value'),
     State('bgk-type', 'value')],
    prevent_initial_call=True
)
def generate_background_knowledge(n_clicks, cpdag_input, num_constraints, bgk_type):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # Parse input CPDAG
        input_dict = eval(cpdag_input)
        
        # Convert dictionary to NetworkX graph
        cpdag = nx.Graph(input_dict)  # 使用无向图因为CPDAG是无向的
        
        # Map bgk_type to numeric value
        type_map = {
            'direct': 1,      # 直接因果关系
            'ancestral': 3,   # 祖先关系（包括间接因果关系）
            'non-ancestral': 2  # 非因果关系
        }
        
        # Generate background knowledge
        bgk = genConsistentCbgk(cpdag, num_constraints, type_map[bgk_type])
        
        # Create visualization
        image_data = plot_graph_with_bgk(cpdag, bgk, title="CPDAG with Background Knowledge")
        
        # Create figure
        fig = create_image_figure(image_data, "CPDAG with Background Knowledge")
        
        # Create results text
        results = html.Div([
            html.H4("Generated Background Knowledge:", style={"marginBottom": "10px"}),
            html.P(f"Generated {len(bgk)} background knowledge constraints:"),
            html.Ul([
                html.Li([
                    f"({x}, {y}): ",
                    html.Span(
                        {
                            1: "Cause",
                            0: "Non-cause"
                        }[t],
                        style={"color": {
                            1: "#2ecc71",  # 绿色表示因果关系
                            0: "#e74c3c"   # 红色表示非因果关系
                        }[t]}
                    )
                ]) for x, y, t in bgk
            ]),
            # 添加背景知识的代码格式显示
            html.P("Background knowledge in code format:"),
            html.Pre(
                "[\n" + 
                "\n".join([f"    ('{x}', '{y}', {t}),  # {x} {'is' if t==1 else 'is NOT'} a cause of {y}" 
                          for x, y, t in bgk]) + 
                "\n]"
            ),
            html.P("Background knowledge can be used in CPDAG to MPDAG conversion.")
        ])
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

# Add background knowledge analysis section
html.Div([
    html.H3("Background Knowledge Analysis"),
    html.Div([
        html.Label("Input DAG:"),
        dcc.Textarea(
            id='bgk-dag-input',
            placeholder="Enter DAG as dictionary",
            value="{0:{1,2}, 1:{2,3}, 2:{3}, 3:{}}",
            style={'width': '100%', 'height': 100}
        ),
        html.Label("Background Knowledge:"),
        dcc.Textarea(
            id='bgk-analysis-input',
            placeholder="Enter background knowledge as list of tuples",
            value="[(0,1,1), (1,2,2), (2,3,3)]",
            style={'width': '100%', 'height': 100}
        ),
        html.Label("Analysis Type:"),
        dcc.Dropdown(
            id='bgk-analysis-type',
            options=[
                {'label': 'Check Compatibility', 'value': 'compatibility'},
                {'label': 'Check Equivalence', 'value': 'equivalence'},
                {'label': 'Decompose BGK', 'value': 'decompose'}
            ],
            value='compatibility'
        ),
        html.Button(
            'Analyze BGK',
            id='analyze-bgk-btn',
            n_clicks=0,
            style=BUTTON_STYLE
        ),
    ]),
    dcc.Loading(
        id="loading-bgk-analysis",
        children=[
            html.Div(id='bgk-analysis-results'),
            dcc.Graph(id='bgk-analysis-graph')
        ]
    )
], style={'marginBottom': '30px'})

# Add callback for background knowledge analysis
@callback(
    [Output('bgk-analysis-results', 'children'),
     Output('bgk-analysis-graph', 'figure')],
    [Input('analyze-bgk-btn', 'n_clicks')],
    [State('bgk-dag-input', 'value'),
     State('bgk-analysis-input', 'value'),
     State('bgk-analysis-type', 'value')],
    prevent_initial_call=True
)
def analyze_background_knowledge(n_clicks, dag_input, bgk_input, analysis_type):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # Parse inputs
        dag_dict = eval(dag_input)
        bgk = eval(bgk_input)
        
        # Convert to NetworkX graph
        dag = nx.DiGraph(dag_dict)
        
        if analysis_type == 'compatibility':
            # Check if background knowledge is compatible with DAG
            is_compatible = check_bgk_compatibility(dag, bgk)
            
            # Create visualization
            edge_colors = {}
            for x, y, t in bgk:
                if (x, y) in dag.edges():
                    edge_colors[(x, y)] = '#2ecc71' if is_compatible else '#e74c3c'
            
            image_data = plot_graph(dag, 
                                  title="DAG with Background Knowledge",
                                  edge_colors=edge_colors)
            
            results = html.Div([
                html.H4("Compatibility Analysis Results:", style={"marginBottom": "10px"}),
                html.P("Background knowledge is " + 
                      ("compatible" if is_compatible else "not compatible") +
                      " with the DAG."),
                html.P("Incompatible constraints:") if not is_compatible else None,
                html.Ul([
                    html.Li(f"({x}, {y}, {t})") 
                    for x, y, t in bgk 
                    if not check_single_constraint(dag, x, y, t)
                ]) if not is_compatible else None
            ])
            
        elif analysis_type == 'equivalence':
            # Check if two sets of background knowledge are equivalent
            bgk2 = generate_equivalent_bgk(bgk)  # Generate an equivalent set
            is_equivalent = check_bgk_equivalence(bgk, bgk2)
            
            results = html.Div([
                html.H4("Equivalence Analysis Results:", style={"marginBottom": "10px"}),
                html.P("Generated equivalent background knowledge:"),
                html.Pre(json.dumps(bgk2, indent=2)),
                html.P("The two sets are " + 
                      ("equivalent" if is_equivalent else "not equivalent"))
            ])
            
            image_data = plot_graph(dag, title="DAG")
            
        else:  # decompose
            # Decompose background knowledge into direct causes and subclauses
            direct_causes, subclauses = decompose_bgk(bgk)
            
            results = html.Div([
                html.H4("Decomposition Results:", style={"marginBottom": "10px"}),
                html.P("Direct causal relationships:"),
                html.Ul([html.Li(f"({x}, {y})") for x, y in direct_causes]),
                html.P("Subclauses:"),
                html.Ul([html.Li(str(clause)) for clause in subclauses])
            ])
            
            image_data = plot_graph(dag, title="DAG")
        
        # Create figure
        fig = create_image_figure(image_data, "Background Knowledge Analysis")
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

# Add helper functions for background knowledge analysis

def check_bgk_compatibility(dag, bgk):
    """Check if background knowledge is compatible with DAG."""
    return all(check_single_constraint(dag, x, y, t) for x, y, t in bgk)

def check_single_constraint(dag, x, y, t):
    """Check if a single background knowledge constraint is compatible with DAG."""
    if t == 1:  # 因果关系
        return nx.has_path(dag, x, y)
    elif t == 0:  # 非因果关系
        return not nx.has_path(dag, x, y)
    return False

def generate_equivalent_bgk(bgk):
    """Generate an equivalent set of background knowledge."""
    # This is a simplified version that just adds transitive relations
    result = bgk.copy()
    direct_causes = [(x, y) for x, y, t in bgk if t == 1]
    
    # Create a graph from direct causes
    G = nx.DiGraph()
    G.add_edges_from(direct_causes)
    
    # Add transitive relations
    for x in G.nodes():
        for y in G.nodes():
            if x != y and nx.has_path(G, x, y):
                result.append((x, y, 2))  # Add ancestor relation
    
    return result

def check_bgk_equivalence(bgk1, bgk2):
    """Check if two sets of background knowledge are equivalent."""
    # Convert both sets to graphs
    G1 = nx.DiGraph()
    G2 = nx.DiGraph()
    
    # Add edges from direct causes
    for x, y, t in bgk1:
        if t == 1:
            G1.add_edge(x, y)
    for x, y, t in bgk2:
        if t == 1:
            G2.add_edge(x, y)
    
    # Check if the transitive closures are the same
    tc1 = nx.transitive_closure(G1)
    tc2 = nx.transitive_closure(G2)
    
    return nx.is_isomorphic(tc1, tc2)

def decompose_bgk(bgk):
    """Decompose background knowledge into direct causes and subclauses."""
    direct_causes = []
    subclauses = []
    
    for x, y, t in bgk:
        if t == 1:
            direct_causes.append((x, y))
        else:
            subclauses.append((x, y, t))
    
    return direct_causes, subclauses

# Add graph type dropdown to size evaluation section
html.Div([
    html.Label("Graph Type:"),
    dcc.Dropdown(
        id='graph-type-input',
        options=[
            {'label': 'CPDAG', 'value': 'cpdag'},
            {'label': 'DAG', 'value': 'dag'}
        ],
        value='cpdag',
        style={"width": "200px", "marginBottom": "10px"}
    )
])

# 添加随机权重生成的回调函数
@callback(
    Output('ida-weights-input', 'value'),
    [Input('generate-random-weights-btn', 'n_clicks'),
     State('ida-graph-input', 'value')],
    prevent_initial_call=True
)
def generate_random_weights(n_clicks, graph_input):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # 解析图结构
        graph_dict = eval(graph_input)
        
        # 为每条边生成随机权重，使用字符串作为键
        weights = {}
        for node, neighbors in graph_dict.items():
            for neighbor in neighbors:
                # 使用字符串格式 "node,neighbor" 作为键
                key = f"{node},{neighbor}"
                weights[key] = round(random.uniform(0.5, 2.0), 2)
        
        return json.dumps(weights, indent=2)
    except Exception as e:
        return str(e)

def generate_random_dag(n_nodes, edge_prob):
    """Generate a random DAG with string node labels.
    
    Args:
        n_nodes: Number of nodes
        edge_prob: Probability of edge creation
        
    Returns:
        NetworkX DiGraph with string node labels
    """
    # 生成随机有向图，确保无环
    G = nx.gnp_random_graph(n_nodes, edge_prob, directed=True)
    # 将节点按顺序编号，确保只保留从小到大的边以避免环
    G = nx.DiGraph([(str(u), str(v)) for (u,v) in G.edges() if u < v])
    # 添加孤立节点（如果有的话）
    G.add_nodes_from(str(i) for i in range(n_nodes))
    return G

# 修改原有的回调函数
@callback(
    [Output('causal-dag-results', 'children'),
     Output('causal-dag-graph', 'figure'),
     Output('causal-dag-input', 'value'),
     Output('causal-weights-input', 'value'),
     Output('mcov-store', 'data'),
     Output('dag-store', 'data'),
     Output('cpdag-store', 'data')],
    [Input('generate-causal-dag-btn', 'n_clicks'),
     Input('update-dag-btn', 'n_clicks')],
    [State('causal-dag-input', 'value'),
     State('causal-weights-input', 'value'),
     State('random-dag-nodes', 'value'),
     State('random-dag-edge-prob', 'value')],
    prevent_initial_call=False  # 允许初始加载时触发
)
def update_causal_dag(gen_clicks, update_clicks, dag_input, weights_input, n_nodes, edge_prob):
    ctx = dash.callback_context
    
    try:
        if not ctx.triggered:
            # 初始加载，使用默认值
            button_id = 'initial'
        else:
            button_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        if button_id == 'generate-causal-dag-btn':
            # 生成新的随机DAG
            G = generate_random_dag(n_nodes, edge_prob)
            
            # 添加随机权重
            for u, v in G.edges():
                G[u][v]['weight'] = round(random.uniform(0.5, 2.0), 2)
        else:
            # 使用输入的DAG和权重（包括初始加载和更新按钮）
            dag_dict = json.loads(dag_input)
            weights_dict = json.loads(weights_input)
            
            G = nx.DiGraph()
            G.add_nodes_from(dag_dict.keys())
            for node, neighbors in dag_dict.items():
                for neighbor in neighbors:
                    G.add_edge(node, neighbor)
                    
            for edge_str, weight in weights_dict.items():
                u, v = edge_str.split(',')
                G[u][v]['weight'] = weight
        
        # 计算协方差矩阵（新增）
        mcov = true_cov(G)
        
        # 创建CPDAG（新增）
        cpdag = dag_to_cpdag(G)
        
        # 使用plot_weighted_graph画DAG（显示权重），使用plot_graph画CPDAG
        dag_image = plot_weighted_graph(G, title="")
        cpdag_image = plot_graph(cpdag, title="")
        
        # 创建两个独立的图形显示
        fig = {
            'data': [],
            'layout': {
                'grid': {'rows': 1, 'columns': 2, 'pattern': 'independent'},
                'images': [
                    {
                        'source': dag_image,
                        'xref': "x",
                        'yref': "y",
                        'x': 0,
                        'y': 1,
                        'sizex': 1,
                        'sizey': 1,
                        'sizing': "contain",
                        'layer': "below"
                    },
                    {
                        'source': cpdag_image,
                        'xref': "x2",
                        'yref': "y2",
                        'x': 0,
                        'y': 1,
                        'sizex': 1,
                        'sizey': 1,
                        'sizing': "contain",
                        'layer': "below"
                    }
                ],
                'xaxis': {'visible': False, 'range': [0, 1], 'domain': [0, 0.48]},
                'yaxis': {'visible': False, 'range': [0, 1]},
                'xaxis2': {'visible': False, 'range': [0, 1], 'domain': [0.52, 1]},
                'yaxis2': {'visible': False, 'range': [0, 1]},
                'width': 1200,
                'height': 500,
                'margin': {'l': 20, 'r': 20, 't': 80, 'b': 20},  # 增加顶部边距
                'showlegend': False,
                'annotations': [
                    {
                        'text': "DAG",
                        'xref': "paper",
                        'yref': "paper",
                        'x': 0.24,
                        'y': 1.1,     # 将标题位置向上移动
                        'showarrow': False,
                        'font': {'size': 20}
                    },
                    {
                        'text': "CPDAG",
                        'xref': "paper",
                        'yref': "paper",
                        'x': 0.76,
                        'y': 1.1,     # 将标题位置向上移动
                        'showarrow': False,
                        'font': {'size': 20}
                    }
                ]
            }
        }
        
        # 准备存储数据
        mcov_list = mcov.tolist()
        dag_dict = {str(n): list(G.successors(n)) for n in G.nodes()}
        cpdag_dict = {str(n): list(cpdag.successors(n)) for n in cpdag.nodes()}
        weights_dict = {f"{u},{v}": G[u][v]['weight'] for u, v in G.edges()}
        
        # 将DAG和权重信息合并到一个字典中
        dag_store = {
            'structure': dag_dict,
            'weights': weights_dict
        }
        
        # 创建结果文本
        results = html.Div([
            html.H4("DAG Properties:", style={"marginBottom": "10px"}),
            html.P(f"Number of nodes: {len(G.nodes())}, Number of edges: {len(G.edges())}")
        ])
        
        # 返回所有输出
        return (
            results, 
            fig,
            json.dumps(dag_dict),           # 原有的DAG结构输入框更新
            json.dumps(weights_dict),        # 原有的权重输入框更新
            mcov_list,                      # 新增：存储协方差矩阵
            dag_store,                      # 新增：存储DAG信息
            cpdag_dict                      # 新增：存储CPDAG信息
        )
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}, "{}", "{}", [], {}, {}  # 出错时返回空值

# Step 2: Background Knowledge and Causal Effect Inference
html.Div([
    html.H3("Step 2: Background Knowledge and Causal Effect Inference"),
    html.P([
        "Input background knowledge to obtain MPDAG, then select source and target nodes for causal effect estimation."
    ], style={"marginBottom": "15px"}),
    
    # Background Knowledge Input
    html.Label("Background Knowledge:"),
    html.P([
        "Format: [(x,y,t)], where t=1: x is a cause of y (x → y), t=0: x is NOT a cause of y (x ↛ y)"
    ], style={"color": "#666666", "fontSize": "0.9em"}),
    dcc.Textarea(
        id='causal-bgk-input',
        placeholder="Enter background knowledge (e.g., [('0','1',1), ('1','2',0)])",
        value="[('2','3',1)]",
        style={'width': '100%', 'height': 100}
    ),
    
    # Source and Target Selection
    html.Div([
        html.Label("Source Node (x):"),
        dcc.Input(
            id='causal-source',
            type='text',
            value='1',
            style=INPUT_STYLE
        ),
        html.Label("Target Node (y):"),
        dcc.Input(
            id='causal-target',
            type='text',
            value='3',
            style=INPUT_STYLE
        ),
    ], style={'marginTop': '15px'}),
    
    # Analysis Button
    html.Button(
        'Analyze Causal Effects',
        id='analyze-causal-btn',
        n_clicks=0,
        style=BUTTON_STYLE
    ),
    
    # Results Display
    dcc.Loading(
        id="loading-causal-analysis",
        children=[
            html.Div(id='causal-analysis-results'),
            dcc.Graph(id='causal-mpdag-graph')
        ]
    )
], style={'marginBottom': '30px'})

# Add callback for background knowledge analysis
@callback(
    [Output('causal-analysis-results', 'children'),
     Output('causal-mpdag-graph', 'figure')],
    [Input('analyze-causal-btn', 'n_clicks')],
    [State('causal-bgk-input', 'value'),
     State('causal-source', 'value'),
     State('causal-target', 'value'),
     State('cpdag-store', 'data'),
     State('mcov-store', 'data'),
     State('dag-store', 'data')],  # 添加 DAG 信息
    prevent_initial_call=True
)
def analyze_causal_effects(n_clicks, bgk_input, source, target, cpdag_dict, mcov_list, dag_store):
    if not n_clicks:
        raise dash.exceptions.PreventUpdate
    
    try:
        # 解析输入
        bgk = eval(bgk_input)
        mcov = np.array(mcov_list)
        
        # 创建CPDAG
        cpdag = nx.DiGraph()
        cpdag.add_nodes_from(cpdag_dict.keys())
        for node, neighbors in cpdag_dict.items():
            for neighbor in neighbors:
                cpdag.add_edge(node, neighbor)
        
        # 转换背景知识为DCC
        dcc = transCbgk(cpdag, bgk)
        
        # 构造MPDAG
        mpdag = construct_mpdag(cpdag, dcc)
        
        # 从原始DAG计算真实因果效应
        true_dag = nx.DiGraph()
        true_dag.add_nodes_from(dag_store['structure'].keys())
        for node, neighbors in dag_store['structure'].items():
            for neighbor in neighbors:
                true_dag.add_edge(node, neighbor)
        
        # 创建权重矩阵
        n = len(true_dag.nodes())
        weight_matrix = np.zeros((n, n))
        for edge_str, weight in dag_store['weights'].items():
            u, v = edge_str.split(',')
            weight_matrix[int(u), int(v)] = weight
        
        # 计算协方差矩阵
        mcov = true_cov(true_dag, weight_matrix)
        
        # 计算真实因果效应
        true_effect = bgk_ida(source, target, mcov, true_dag)[0]
        
        # 创建可视化
        # 设置边的颜色（背景知识）
        edge_colors = {}
        for x, y, t in bgk:
            edge_colors[(x, y)] = '#2ecc71' if t == 1 else '#e74c3c'
        
        # 设置节点的颜色（source和target）
        node_colors = {}
        node_colors[source] = '#3498db'  # source节点蓝色
        node_colors[target] = '#e67e22'  # target节点橙色
        
        image_data = plot_graph(mpdag, 
                               title="MPDAG with Background Knowledge",
                               edge_colors=edge_colors,
                               node_colors=node_colors)
        
        fig = create_image_figure(image_data)
        
        # 创建结果文本
        results = html.Div([
            html.H4("Causal Analysis Results:", style={"marginBottom": "10px"}),
            html.P([
                "Source node: ",
                html.Span(source, style={"color": "#3498db", "fontWeight": "bold"})
            ]),
            html.P([
                "Target node: ",
                html.Span(target, style={"color": "#e67e22", "fontWeight": "bold"})
            ]),
            html.P("Possible causal effects based on MPDAG:"),
            html.Ul([html.Li(f"{effect:.4f}") for effect in bgk_ida(source, target, mcov, cpdag, mpdag=mpdag, bgk=bgk)]),
            html.Div([
                html.P([
                    "True causal effect (if the original DAG represents the true causal structure): ",
                    html.Span(f"{true_effect:.4f}", style={"fontWeight": "bold"})
                ])
            ], style={"marginTop": "15px", "padding": "10px", "backgroundColor": "#f8f9fa"})
        ])
        
        return results, fig
        
    except Exception as e:
        return html.Div([
            html.P("Error:", style={"color": "red"}),
            html.P(str(e))
        ]), {}

#if __name__ == '__main__':   
    #app.run_server(debug=True) 
if __name__ == '__main__':
    port = int(os.environ.get("PORT", 8050))
    app.run_server(debug=False, host='0.0.0.0', port=port)
