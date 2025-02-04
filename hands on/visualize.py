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

class VisualizationManager:
    """Class of visualization figure creation"""
    def __init__(self, resultdir, elements, total_atoms):
        self.resultdir = resultdir
        self.elements = elements
        self.total_atoms = total_atoms
        self.candidates = self.get_composition_candidate()
        self.default_layout = {'template': 'plotly_white',
                               'showlegend': True,
                               'margin': dict(t=50, l=50, r=50, b=50),
                               'title_font': dict(size=16),
                               'plot_bgcolor': 'white',
                               'width': 600,
                               'height': 500,
                               'xaxis': dict(title=dict(font=dict(color='black')),
                                             linecolor='black', linewidth=1, mirror=True,
                                             ticks='inside', tickcolor='black', tickwidth=1, ticklen=5,
                                             tickformat=".1f", tickfont=dict(color='black')),
                               'yaxis': dict(title=dict(font=dict(color='black')),
                                             linecolor='black', linewidth=1, mirror=True,
                                             ticks='inside', tickcolor='black', tickwidth=1, ticklen=5,
                                             tickformat=".0f", tickfont=dict(color='black')),
                               'legend_font_color': 'black'}

    def get_available_plots(self) -> Dict[str, callable]:
        """Return available display function"""
        return {'model': self.view_model,
                'composition': self.view_composition,
                'rdf': self.view_rdf,
                'relativeposition': self.view_relativeposition,
                'energy': self.view_energy}

    def get_composition_candidate(self):
        modelpathlist = np.sort(glob.glob(self.resultdir+'*.cif'))
        candidates = []
        for modelpath in modelpathlist:
            filename = modelpath.split('/')[-1]
            candidate = []
            for element in self.elements:
                # 元素名の後の数字を検索
                pattern = f"{element}(\\d+)"
                match = re.search(pattern, filename)
                ratio = int(match.group(1)) if match else 0
                candidate.append(ratio)
            candidates.append(candidate)
        return np.array(candidates)

    def calculate_composition(self,ratios):
        total_percentage = sum(ratios)
        composition = np.array([round((ratio / total_percentage) * self.total_atoms) for ratio in ratios])
        remaining_atoms = self.total_atoms - sum(composition)
        for i in range(abs(remaining_atoms)):
            composition[i % len(ratios)] += 1 if remaining_atoms > 0 else -1
        distances = np.sqrt(np.sum((self.candidates-composition) ** 2, axis=1))
        min_idx = np.argmin(distances)
        compositions = {element: ratio for element, ratio in zip(self.elements, self.candidates[min_idx].tolist())}
        return compositions
    
    def view_model(self, raw_ratios,size=15):
        normalized_compositions = self.calculate_composition(raw_ratios)
        alloy_ratios = [normalized_compositions[element] for element in self.elements]
        modelname = ''.join(f"{element}{ratio:.0f}" for element, ratio in zip(self.elements, alloy_ratios))
        model = read(self.resultdir+f'{modelname}.cif')

        #set viewpoint of plotly figure
        camera = dict(up=dict(x=0, y=0, z=1),
                      center=dict(x=0, y=0, z=0),
                      eye=dict(x=0, y=-2, z=1))    
        #set color of atom        
        colors = []
        for symbol in model.symbols:
            index = np.where(np.array(self.elements) == symbol)[0][0]
            colors.append(plotly.colors.qualitative.Plotly[index])
        #get periodic boundary conditions
        cell = model.get_cell()
        vertices = np.array([[0, 0, 0],
                             [1, 0, 0],
                             [1, 1, 0],
                             [0, 1, 0],
                             [0, 0, 1],
                             [1, 0, 1],
                             [1, 1, 1],
                             [0, 1, 1]]) @ cell
        edges = [(0,1), (1,2), (2,3), (3,0),
                (4,5), (5,6), (6,7), (7,4),
                (0,4), (1,5), (2,6), (3,7)]
        fig = go.Figure()

        #display surface model
        fig.add_trace(go.Scatter3d(x=model.positions[:, 0],
                                   y=model.positions[:, 1],
                                   z=model.positions[:, 2],
                                   mode='markers',
                                   marker=dict(size=size,
                                               color=colors,
                                               symbol='circle',
                                               opacity=1.0,
                                               line=dict(width=2, color='#000000'))))
        #vew periodic boundary conditions
        for start_idx, end_idx in edges:
            start = vertices[start_idx]
            end = vertices[end_idx]
            fig.add_trace(go.Scatter3d(x=[start[0], end[0]],
                                       y=[start[1], end[1]],
                                       z=[start[2], end[2]],
                                       mode='lines',
                                       line=dict(color='black', width=2),
                                       hoverinfo='none'))
        # set layout
        fig.update_layout(title='<b>Surface model</b>',
                          scene=dict(xaxis=dict(visible=False),
                                     yaxis=dict(visible=False),
                                     zaxis=dict(visible=False),
                                     bgcolor="rgba(0,0,0,0)",
                                     camera=camera),
                          paper_bgcolor="rgba(0,0,0,0)",
                          margin=dict(l=0, r=0, b=0, t=0),
                          showlegend=False)
        return fig

    def view_composition(self, raw_ratios):
        normalized_compositions = self.calculate_composition(raw_ratios)
        alloy_ratios = [normalized_compositions[element] for element in self.elements]

        modelname = ''.join(f"{element}{ratio:.0f}" for element, ratio in zip(self.elements, alloy_ratios))
        model = read(self.resultdir+f'{modelname}.cif')
        sortindex = np.argsort(model.positions[:, 2])
        height = model.positions[sortindex][:, 2]
        layer_indices = np.where(np.diff(height) > 0.5)[0] + 1
        layers = np.split(height, layer_indices)

        atom_count = 0
        composition_list = []
        for layer in layers:
            average_height = np.average(layer)
            target_index = sortindex[atom_count:atom_count + len(layer)]
            data = [average_height]
            for element in self.elements:
                ratio = len([1 for i in model[target_index] if i.symbol == element]) *100 / len(layer)
                data.append(ratio)
            atom_count += len(layer)
            composition_list.append(data)
        df = pd.DataFrame(composition_list, columns=['height'] + self.elements)

        fig = go.Figure()
        for element in self.elements:
            fig.add_trace(go.Scatter(x=df['height'],
                                     y=df[element],
                                     mode='lines+markers',
                                     marker=dict(size=13),
                                     line=dict(dash='dot'),
                                     name=element))
        fig.update_layout(title='<b>Composition of Elements by Height</b>',
                          xaxis_title="Height of atoms [angstrom]",
                          yaxis_title="Composition ratio [at%]",
                          legend_title="Elements",
                          **self.default_layout)        
        return fig
    
    def view_rdf(self, raw_ratios, center_element, target_elements):
        normalized_compositions = self.calculate_composition(raw_ratios)
        alloy_ratios = [normalized_compositions[element] for element in self.elements]

        # Get information
        modelname = ''.join(f"{element}{alloy_ratio}" for element, alloy_ratio in zip(self.elements, alloy_ratios))
        model = read(self.resultdir+f'{modelname}.cif')
        sortindex = np.argsort(model.positions[:, 2])
        height = model.positions[sortindex][:, 2]
        layer_indices = np.where(np.diff(height) > 0.5)[0] + 1
        layers = np.split(height, layer_indices)
        layers_index = np.split(sortindex, layer_indices)

        # Create dataset for RDF
        rdf_data = {}
        for i in reversed(range(len(layers))):
            layer_model = model.copy()
            del layer_model[:]
            for atom in model[layers_index[i]]:
                layer_model += Atom(atom.symbol, atom.position)
            
            # Calculate RDF for all target elements
            all_rdfs = []
            for target_element in target_elements:
                bin_centers, rdf = self.calculate_rdf(layer_model * (3, 3, 1), center_element, target_element, num_bins=100, r_max=5)
                all_rdfs.append(rdf)

            # Average RDF for the current layer
            average_rdf = np.mean(all_rdfs, axis=0)
            rdf_data[f'Layer {i + 1}'] = (bin_centers, average_rdf)

        # Create figure
        fig = go.Figure()
        for layer_name, (bin_centers, average_rdf) in rdf_data.items():
            fig.add_trace(go.Scatter(x=bin_centers, 
                                     y=average_rdf, 
                                     mode='lines', 
                                     name=f'{layer_name}'))
        fig.update_layout(title=f'<b>RDF of {center_element} to {", ".join(target_elements)} by Layer</b>',
                          xaxis_title='Distance [angstrom]',
                          yaxis_title='RDF',
                          legend_title="Layers",
                          **self.default_layout)
        return fig

    def calculate_rdf(self, atoms, center_element, target_elements, num_bins=100, r_max=10):
        center_indices = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == center_element]
        target_indices = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == target_elements]

        distances = []
        
        for center in center_indices:
            for target in target_indices:
                distance = np.linalg.norm(atoms.positions[center] - atoms.positions[target])
                if 0 < distance < r_max:
                    distances.append(distance)

        hist, bin_edges = np.histogram(distances, bins=num_bins, range=(0, r_max))
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        rdf = hist / (4/3 * np.pi * (bin_centers**3) * (bin_edges[1] - bin_edges[0]))

        return bin_centers, rdf
    
    def view_relativeposition(self, raw_ratios):
        #Get the difference of median for each metal ad the relative position
        normalized_compositions = self.calculate_composition(raw_ratios)
        alloy_ratios = [normalized_compositions[element] for element in self.elements]

        csv_path = self.resultdir+'relativeposition_data.csv'
        if not os.path.exists(csv_path):
            modelpathlist = np.sort(glob.glob(self.resultdir+'*.cif'))
            data = []
            for modelpath in modelpathlist:
                model = read(modelpath)
                median_position_list = []
                composition_list = []
                for element in self.elements:
                    height = [i.position[2] for i in model if i.symbol == element]
                    if height:
                        median_position_list.append(np.median(height))
                        composition_list.append(len(height))
                    else:
                        median_position_list.append(np.nan)
                        composition_list.append(0)
                data.append(median_position_list + composition_list)
            df = pd.DataFrame(data, columns=[element + '_position' for element in self.elements] + [element + '_count' for element in self.elements])
            df.to_csv(csv_path)

        #Create figure
        fig = go.Figure()
        element_combination_list = list(itertools.combinations(range(len(self.elements)), 2))
        num_bins = 10
        df=pd.read_csv(csv_path)
        range_max=max([df['{}_position'.format(element)].max() for element in self.elements])

        targets = []
        for i in range(len(self.elements)):
            targets.append(f"{self.elements[i]}_count == {alloy_ratios[i]}")
        query_str = ' & '.join(targets)
        target_positions = np.array(df.query(query_str)[[col for col in df.columns if col.endswith('_position')]])[0]

        for index, element_combination in enumerate(element_combination_list):
            x = df[self.elements[element_combination[0]] + '_count'] / (df[self.elements[element_combination[0]] + '_count'] + df[self.elements[element_combination[1]] + '_count'])
            y = df[self.elements[element_combination[0]] + '_position'] - df[self.elements[element_combination[1]] + '_position']

            bins = np.linspace(x.min(), x.max(), num_bins + 1)
            x_bin_centers = []
            y_means = []
            y_maxs = []
            y_mins = []
            for i in range(num_bins):
                bin_mask = (x >= bins[i]) & (x < bins[i + 1])
                if bin_mask.sum() > 0:
                    x_bin_centers.append((bins[i] + bins[i + 1]) / 2)
                    y_means.append(y[bin_mask].mean())
                    y_maxs.append(y[bin_mask].max())
                    y_mins.append(y[bin_mask].min())

            default_colors = plotly.colors.qualitative.Plotly
            mean_color = default_colors[index % len(default_colors)]
            r = int(mean_color[1:3], 16)
            g = int(mean_color[3:5], 16)
            b = int(mean_color[5:7], 16)
            fill_color = f'rgba({r}, {g}, {b}, 0.3)' 

            fig.add_trace(go.Scatter(x=x_bin_centers + x_bin_centers[::-1],
                                     y=y_maxs + y_mins[::-1],
                                     fill="toself",
                                     fillcolor=fill_color,
                                     line=dict(color="rgba(255,255,255,0)"),
                                     name='{0}/{1} (range)'.format(self.elements[element_combination[0]], self.elements[element_combination[1]])))
            fig.add_trace(go.Scatter(x=x_bin_centers,
                                     y=y_means,
                                     mode='lines+markers',
                                     line=dict(color=mean_color, dash='dot'),
                                     marker=dict(color=mean_color, size=13),
                                     name='{0}/{1} (ave.)'.format(self.elements[element_combination[0]], self.elements[element_combination[1]])))
            if (alloy_ratios[element_combination[0]]+alloy_ratios[element_combination[1]]) not in [0, 1]:
                fig.add_trace(go.Scatter(x=[alloy_ratios[element_combination[0]]/(alloy_ratios[element_combination[0]]+alloy_ratios[element_combination[1]])],
                                         y=[target_positions[element_combination[0]]-target_positions[element_combination[1]]],
                                         mode='markers',
                                         marker=dict(color=mean_color, line=dict(color='black',width=2),size=15),
                                         name='{0}/{1} (now)'.format(self.elements[element_combination[0]], self.elements[element_combination[1]])))
        fig.update_layout(title="<b>Relative position</b>",
                          xaxis_title="Alloy ratio M<sub>1</sub>/(M<sub>1</sub>+M<sub>2</sub>) [-]",
                          yaxis_title="Relative height [angstrom]",
                          legend_title="Combinations",
                          **self.default_layout)
        return fig
    
    def view_energy(self, raw_ratios, visualize_start=12, visualize_width=4):
        normalized_compositions = self.calculate_composition(raw_ratios)
        alloy_ratios = [normalized_compositions[element] for element in self.elements]

        energy_csvpath = self.resultdir+'energy.csv'
        df = pd.read_csv(energy_csvpath)
        fig = make_subplots(rows=len(self.elements)-1, cols=len(self.elements)-1, horizontal_spacing=0.05, vertical_spacing=0.05)
        
        colorlist = plotly.colors.qualitative.Plotly[:len(range(visualize_start, self.total_atoms+1, visualize_width))] * (len(self.elements) * (len(self.elements) - 1) // 2)
        showlegend = True
        color_count = 0

        targets = []
        for i in range(len(self.elements)):
            targets.append(f"{self.elements[i]} == {alloy_ratios[i]}")
        query_str = ' & '.join(targets)
        target_energy = np.array(df.query(query_str)['energy'])

        for i in range(len(self.elements)-1):
            for j in range(i, len(self.elements)-1):
                x = np.array(df[self.elements[j+1]] / (df[self.elements[i]] + df[self.elements[j+1]]))
                size = (df[self.elements[i]] + df[self.elements[j+1]])
                for unit in range(visualize_start, self.total_atoms+1, visualize_width):
                    index = np.where((np.array(size) == unit))[0]
                    fig.add_trace(go.Scatter(x=np.sort(x[index]), 
                                             y=np.array(df['energy'][index])[np.argsort(x[index])],
                                             mode='lines',
                                             line=dict(width=3, color=colorlist[color_count]),
                                             name='N=' + str(unit),
                                             showlegend=showlegend, 
                                             legendgroup=unit), 
                                  row=j+1, col=i+1)
                    color_count += 1
                    fig.update_yaxes(range=(df['energy'].min(), df['energy'].max()), row=j+1, col=i+1)
                showlegend = False
                if (alloy_ratios[i]+alloy_ratios[j+1]) not in [0, 1]:
                    fig.add_trace(go.Scatter(x=[alloy_ratios[j+1]/(alloy_ratios[i]+alloy_ratios[j+1])],
                                                y=target_energy,
                                                mode='markers',
                                                marker=dict(color='grey', line=dict(color='black',width=2),size=15),
                                                name='{0}/{1} (now)'.format(self.elements[j+1],self.elements[i])),
                                    row=j+1, col=i+1)
                        
        
        for i, name in enumerate(self.elements[:-1]):
            fig.update_xaxes(title='<b>/({0}+M)'.format(name), row=len(self.elements)-1, col=i+1)
            
        for j, name in enumerate(self.elements[1:]):
            fig.update_yaxes(title='<b>M=' + name, row=j+1, col=1)

        #set layout
        for i in range(len(self.elements)-1):
            for j in range(len(self.elements)-1):
                fig.update_xaxes(title=dict(font=dict(color='black')),
                                 linecolor='black', linewidth=1, mirror=True,
                                 ticks='inside', tickcolor='black', tickwidth=1, ticklen=5,
                                 tickfont=dict(color='black'),
                                 range=(0, 1), row=j+1, col=i+1)

                fig.update_yaxes(title=dict(font=dict(color='black')),
                                 linecolor='black', linewidth=1, mirror=True,
                                 ticks='inside', tickcolor='black', tickwidth=1, ticklen=5,
                                 tickfont=dict(color='black'), row=j+1, col=i+1)

        fig.update_layout(title="<b>Surface energy</b>",
                          title_font=dict(size=16),
                          legend_title="size",
                          plot_bgcolor='white',
                          width=600,
                          height=500,
                          margin=dict(t=50, b=50, l=50, r=50),
                          legend=dict(x=1, y=1, xanchor='right', yanchor='top', bgcolor="rgba(255,255,255,0.8)",
                                      traceorder="normal",itemsizing="constant",title_font_size=12,font=dict(size=10),orientation="v", tracegroupgap=5),
                          legend_font_color='black')
        return fig