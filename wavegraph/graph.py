# (C) 2014-2018
# Contributed to by Eve Chase, Eric Chassande-Mottin, Eric Lebigot, Philippe Bacon, Quentin Bammey

import os
import sys
import time
import glob
import numpy
import getpass
import collections
import itertools
import logging
import subprocess
import yaml

# Information about a GridPoint in a Graph:
GraphInfo = collections.namedtuple('GraphNodeInfo', 'value_avg value_stdev')

class Graph(object):
    """
    Directed Acyclic Graph of GridPoints endowed with GraphInfo information.
    The graph is built by grouping multiple GridPaths: together, they define the graph.

    Main attributes:
    ---------------
    graph -- dictionary that maps each GridPoint to the set of its ancestor 
    GridPoints.
    point_info -- dictionary that maps each GridPoint to additional 
    information (GraphInfo object).
    head_nodes -- the set of 'head' nodes, i.e. nodes which are not an ancestor 
    in some of the original clusters.
    sorted_list -- a topologically sorted list of the GridPoints
    """
    def __init__(self, clusters):
        """
        Create a new graph from clusters.

        Attributes:
        -----------
        clusters [list] -- Cluster objects list that will make up the graph
        """
        self.graph = collections.defaultdict(set)
        self.point_info = {}
        self.head_nodes = set()
        
        # Values at each GridPoint are first stored: they will later
        # yield an average and a standard deviation:
        values = collections.defaultdict(list)
        
        for cluster in clusters:
            
            cluster_points = cluster.grid_points
            
            # Test if cluster is empty
            if not cluster_points:
                continue
            
            # Test if cluster has duplicated nodes
            if len(cluster_points) != len(set(cluster_points)):
                logging.warning("Cluster {} has duplicates".format(cluster.metadata))
                
            # The first node has no ancestors (by definition of
            # "ancestor"):
            self.graph[cluster_points[0]] = set()
            
            for (ancestor, node) in itertools.izip(cluster_points[:-1],
                                                   cluster_points[1:]):
                self.graph[node].add(ancestor)
                
            # The values contributed to by the cluster are registered
            # at each of its points:
            for (grid_point, value) in itertools.izip(cluster_points,
                                                      cluster.values):
                values[grid_point].append(value)
                
            # The last node in a cluster is by definition a head node:
            self.head_nodes.add(node)
            
        # Statistics on the values at all grid point:
        for (grid_point, values) in values.iteritems():
            logging.debug("Graph point: {}".format(grid_point))
            logging.debug("  Values: {}".format(values))
            
            self.point_info[grid_point] = GraphInfo(
                numpy.mean(values), numpy.std(values))
            
            logging.debug("  Point info: {}".format(
                self.point_info[grid_point]))
        
        self.sorted_list = self._topological_sorting()
        
    def __str__(self, grid):
        """
        Returns a string representing the graph.
        Each node is assigned an ID.
        
        For each node, print its ID, time index, frequency index,
        scale index, endnode boolean, and the node ID's of all
        ancestors, in that order. The endnode boolean is 0 if this
        node is an ancestor for another node in the same cluster or 1
        otherwise. In a comment line, it also adds the time and frequency 
        coordinates in the node in physical units.

        Input:
        ------
        grid [CoherentWaveBurstGrid object] -- cWB grid

        Output:
        -------
        string [str] -- text version of selected pixels in Cluster object
        """
        # Init with the header of the string description
        output_list = [
            "## nodeID time_idx freq_idx scale_idx" \
            " value_avg value_stdev" \
            " endnode ancestors"
        ]
        
        node_IDs = {}
        
        # Loop over the node in the graph sorted in the topological order
        for (node_ID, node) in enumerate(self.sorted_list):
        
            extra_info = self.point_info[node]
        
            # Collect the relevant informations about node in a list
            node_info = [
                node_ID,
                node.time_index, node.freq_index, node.scale_index,
                extra_info.value_avg, extra_info.value_stdev
            ]
        
            node_info.append(1 if node in self.head_nodes else 0)
        
            # Retrieve the node IDs of the ancestors for node
            for ancestor in self.graph[node]:
                node_info.append(node_IDs[ancestor])
            
            # Convert node infos to a string
            output_list.append(' '.join(map(str, node_info)))
        
            # Add a comment line with node coordinates in physical units
            output_list.append('## {ID}: a={scale}, t={time} s, f={frequency} Hz' \
                               .format(ID=node_ID, **node.phys_coords(grid)))
            
            node_IDs[node] = node_ID
        
        return '\n'.join(output_list)

    def _topological_sorting(self):
        """
        Return a topologically sorted list of nodes, representing the graph.
        The first nodes have no ancestors.
        The graph must have no cycles.

        Output:
        -------
        result [list] -- sorted list of nodes
        """
        # Initially add each node to a list of unvisited nodes.
        nodes_unvisited = set(self.graph)
    
        def sorted_ancestors(node):
            """
            Return a topological sorting of a node and of its
            ancestors at all levels (i.e., recursively).
            node -- a tuple of GridPoint objects
            """
            result = []
        
            #  If an ancestor has not been visited, sort that ancestor
            #  and all of its ancestors.
            for ancestor in self.graph[node]:
                try:
                    nodes_unvisited.remove(ancestor)
                except KeyError:
                    pass
                else:
                    result.extend(sorted_ancestors(ancestor))
                
            result.append(node)
        
            return result
    
    
        sorted_nodes = []
        while nodes_unvisited:
            sorted_nodes.extend(sorted_ancestors(nodes_unvisited.pop()))
        
        return sorted_nodes

    def span(self, grid):
        """
        Calculate the span index for a given graph in units of the largest
        timescales.
        
        Input:
        -------
        grid [CoherentWaveBurstGrid object] -- cWB grid
        
        Output:
        -------
        span [int] -- coverage of the graph
        """
        #  Express the maximum time in any timescale plane as its index in the
        #  largest timescale plane
        scale_index_max = max(node.scale_index for node in self.sorted_list)
        time_max = max(node.time_index * grid.timescales[node.scale_index]
                       for node in self.sorted_list)

        return int(numpy.floor(time_max/grid.timescales[scale_index_max]) + 1)

def read_graph():
    """
    Read graph from .txt file.
    """
    raise NotImplementedError

def write_graph(filename, graph, metadata, grid, format='custom'):
    """
    Write graph to .txt file.

    Input
    -----
    filename       [str]  -- name of output file
    graph [Graph object] -- graph to be stored
    metadata       [str] -- additional info about input data 

    Output file structure:
    ---------------------
    metadata in a header
    nodeID time_idx freq_idx scale_idx value_avg value_stdev endnode ancestors
    """
    # Write in output .txt file.
    with open(filename, 'w') as outfile:
        print >> outfile, metadata
        print >> outfile, graph.__str__(grid)
        
    logging.info('Wrote graph ({} nodes) in {}'.format(len(graph.sorted_list), \
                                                       filename))
