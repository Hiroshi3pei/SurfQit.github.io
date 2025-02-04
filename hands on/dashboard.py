import dash
from dash import html, dcc
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import os
import glob
import re
import itertools
from typing import Dict, List, Any, Optional
import math
from ase.io import read,write
from ase import Atom
import nglview
import ipywidgets as widgets
from IPython.display import display
import asyncio
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.colors
import visualize

class DashApp:
    """Class of Dash Management"""
    def __init__(self, elements, total_atoms, visualization_manager: visualize.VisualizationManager, initial_plots: Optional[List[str]] = None):
        self.elements = elements
        self.total_atoms = total_atoms
        self.viz_manager = visualization_manager
        self.app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
        self.plot_configs = self.viz_manager.get_available_plots()
        self.initial_plots = self._validate_plot_selection(initial_plots or [])
        self.setup_layout()
        self.setup_callbacks()
    
    def _validate_plot_selection(self, plots: List[str]) -> List[str]:
        available_plots = set(self.viz_manager.get_available_plots().keys())
        return [p for p in plots if p in available_plots]
    
    def setup_layout(self):
        INITIAL_RATIOS = [self.total_atoms, 0, 0]
        # Sidebar component with Plot Selector, Sliders and Dropdowns
        sidebar = html.Div([dbc.Row([html.H4('Settings', style={"font-size": "16px"})], style={"height": "5vh"}),
                            dbc.Row([html.Div([html.H5('Element Selection', style={"text-align": "center", "font-size": "14px", "margin-bottom": "15px"}),
                                               *[html.Div([dbc.Row([dbc.Col([html.Label(element, style={"font-size": "12px"})],width=2, style={"text-align": "right", "padding-right": "0px"}),
                                                                    dbc.Col([html.Div(dcc.Slider(id=f'slider-{element.lower()}',
                                                                                                 min=0,
                                                                                                 max=100,#self.total_atoms,
                                                                                                 step=1,
                                                                                                 value=(ratio/self.total_atoms)*100,#ratio,
                                                                                                 marks={i: str(i) for i in range(0, 101, 20)},
                                                                                                 included=True,
                                                                                                 tooltip={"placement": "bottom", "always_visible": False}),
                                                                                      style={'width': '100%', 'margin': '0 5px'})], width=8),
                                                                    dbc.Col([html.Div(id=f'value-{element.lower()}', style={'text-align': 'right', "font-size": "12px"})], width=2, style={"padding-left": "0px"}),], className='mb-1')])
                                                for element, ratio in zip(self.elements, INITIAL_RATIOS)],
                                               html.Div(id='composition-display', style={'margin-top': '20px', "margin-bottom": "30px", "font-size": "14px"}),
                                               html.H5('Plot Selection', style={"text-align": "center", "font-size": "14px", "margin-bottom": "15px"}),
                                               dcc.Dropdown(id='plot-selector',
                                                            options=[{'label': k, 'value': k} for k in self.plot_configs.keys()],
                                                            value=self.initial_plots,
                                                            multi=True,
                                                            style={"font-size": "12px", "margin-bottom": "30px"}),
                                               html.H5('RDF Analysis', style={"text-align": "center", "font-size": "14px", "margin-bottom": "15px"}),
                                               dbc.Row([dbc.Col([html.Label("Center Element", style={"font-size": "12px"}),
                                                                 dcc.Dropdown(id='center-element-dropdown',
                                                                              options=[{'label': elem, 'value': elem} for elem in self.elements],
                                                                              value=self.elements[0],
                                                                              style={"font-size": "12px"})], width=12, className='mb-3'),
                                                        dbc.Col([html.Label("Target Elements", style={"font-size": "12px"}),
                                                                 dcc.Dropdown(id='target-elements-dropdown',
                                                                              options=[{'label': elem, 'value': elem} for elem in self.elements],
                                                                              value=[self.elements[0]],
                                                                              multi=True,
                                                                              style={"font-size": "12px"})], width=12)])])], style={"height": "95vh", "overflow-y": "auto"})], className='bg-light')
        # Content area for graphs
        content = html.Div(id='plot-container', style={"height": "100vh"})
        # Main layout
        self.app.layout = html.Div([dbc.Container([dbc.Row([dbc.Col([html.H1(className="mb-4")])]),
                                    dbc.Row([dbc.Col(sidebar, width=3, className='bg-light'),
                                             dbc.Col(content, width=9)])], fluid=True, className="vh-100 p-0")])
    
    def _get_plot_dimensions(self, n_plots: int) -> tuple:
        """get the number of figure for each line"""
        if n_plots <= 3:
            return 1, n_plots
        else:
            rows = (n_plots + 2) // 3
            if n_plots > 3:
                # 下の段は最大2つ
                lower_row_plots = min(2, n_plots - 3)
                return rows + (1 if lower_row_plots > 0 else 0), min(3, n_plots - lower_row_plots)
            return rows, min(3, n_plots)

    def _get_column_width(self, n_plots: int, plot_index: int) -> int: 
        """get the width of column""" 
        if n_plots <= 3: 
            return math.floor(12 / n_plots) 
        else: 
            if plot_index < 3: 
                return math.floor(12 / 3) 
            else: 
                remaining_plots = n_plots - 3 
                plots_per_row = min(2, remaining_plots)
                return math.floor(12 / plots_per_row)

    def _get_plot_size(self, plot_name: str, n_rows: int, is_single_plot: bool, row_index: int) -> tuple:
        base_height = 1400 // n_rows
        
        if row_index == 0:
            height_factor = 1  
        else: 
            height_factor = 1.2 

        size_configs = {
            "model": {
                "single": (800, int(500 * height_factor)),
                "multi": (int(800 * height_factor), int(600 * height_factor))},
            "composition": {
                "single": (800, int(500 * height_factor)),
                "multi": (int(800 * height_factor), int(600 * height_factor))},
            "rdf": {
                "single": (800, int(500 * height_factor)),
                "multi": (int(800 * height_factor), int(600 * height_factor))},
            "relativeposition": {
                "single": (800, int(500 * height_factor)),
                "multi": (int(800 * height_factor), int(600 * height_factor))},
            "energy": {
                "single": (800, int(800 * height_factor)),
                "multi": (int(800 * height_factor), int(600 * height_factor))}}
        
        default_size = (800, int(base_height * height_factor))
        plot_config = size_configs.get(plot_name, {"single": default_size, "multi": default_size})
        size = plot_config["single"] if is_single_plot else plot_config["multi"]
        if size[1] > base_height * height_factor and not is_single_plot:
            ratio = size[0] / size[1]
            size = (int(base_height * height_factor * ratio), int(base_height * height_factor))
        return size

    def _organize_plots(self, selected_plots: List[str], ratios: List[int], center_element: str = None, target_elements: List[str] = None) -> List[List[dict]]:
        if not selected_plots:
            return []

        n_plots = len(selected_plots)
        n_rows, plots_per_row = self._get_plot_dimensions(n_plots)
        is_single_plot = n_plots == 1

        rows = []
        current_row = []
        plot_index = 0

        for plot_name in selected_plots:
            figure = self._generate_plot(plot_name, 
                                         ratios, 
                                         center_element=center_element, 
                                         target_elements=target_elements,
                                         size=self._get_plot_size(plot_name, n_rows, is_single_plot, len(rows)))

            col_width = self._get_column_width(n_plots, plot_index)
            plot_info = {'figure': figure,
                         'width': col_width,
                         'id': f'plot-{len(rows)}-{len(current_row)}'}
            
            current_row.append(plot_info)
            plot_index += 1
            if len(current_row) == (3 if plot_index <= 3 else 2) or plot_name == selected_plots[-1]:
                if len(current_row) < (3 if plot_index <= 3 else 2):
                    remaining_width = 12 - sum(plot['width'] for plot in current_row)
                    additional_width = remaining_width // len(current_row)
                    for plot in current_row:
                        plot['width'] += additional_width
                rows.append(current_row)
                current_row = []
        return rows

    def _generate_plot(self, plot_name: str, ratios: List[int], center_element: str = None, target_elements: List[str] = None, size: tuple = None) -> go.Figure:
        figure = None
        width, height = size or (600, 500)

        if plot_name == "composition":
            figure = self.plot_configs[plot_name](raw_ratios=ratios)
        elif plot_name == "rdf":
            figure = self.plot_configs[plot_name](raw_ratios=ratios,
                                                  center_element=center_element,
                                                  target_elements=target_elements)
        elif plot_name == "model":
            figure = self.plot_configs[plot_name](raw_ratios=ratios)
        elif plot_name == "relativeposition":
            figure = self.plot_configs[plot_name](raw_ratios=ratios)
        elif plot_name == "energy":
            visualize_start=0
            visualize_width=3
            figure = self.plot_configs[plot_name](raw_ratios=ratios,
                                                  visualize_start=visualize_start,
                                                  visualize_width=visualize_width)
            
        # レイアウトの更新
        if figure:
            figure.update_layout(width=width,
                                 height=height,
                                 margin=dict(l=40, r=40, t=40, b=40),
                                 autosize=False)
        return figure

    def setup_callbacks(self):
        for element in self.elements:
            self.app.callback(Output(f'value-{element.lower()}', 'children'),
                              Input(f'slider-{element.lower()}', 'value'))(lambda value, element=element: str(value)+'%')

        @self.app.callback(Output('composition-display', 'children'),
                           [Input(f'slider-{element.lower()}', 'value') for element in self.elements])
        def update_composition_display(*ratios):
            if sum(ratios) == 0:
                return html.Div("Please adjust the ratios", style={'color': 'red', 'text-align': 'center'}) 
            compositions = self.viz_manager.calculate_composition(ratios)
            composition_text = ''.join(f"{element}{compositions[element]}" for element in self.elements)
    
            return html.Div([html.Div([html.Strong(composition_text, style={'font-size': '24px'})], style={'text-align': 'center'})])

        @self.app.callback(Output('plot-container', 'children'),
                           [Input('plot-selector', 'value'),
                           Input('center-element-dropdown', 'value'),
                           Input('target-elements-dropdown', 'value')] +
                           [Input(f'slider-{element.lower()}', 'value') for element in self.elements])
        
        def update_layout(selected_plots, center_element, target_elements, *ratios):
            if not selected_plots:
                return html.Div("図を選択してください", className="text-center mt-4")

            plot_rows = self._organize_plots(selected_plots,
                                             list(ratios),
                                             center_element=center_element,
                                             target_elements=target_elements)
            row_components = []
            for row in plot_rows:
                cols = []
                for plot in row:
                    cols.append(dbc.Col(dcc.Graph(figure=plot['figure'], id=plot['id']),
                                                  width=plot['width']))
                row_components.append(dbc.Row(cols))
            
            return row_components

        # Run the Dash application
        self.app.run_server(debug=True,port=8051)