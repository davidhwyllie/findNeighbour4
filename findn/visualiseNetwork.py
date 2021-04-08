""" visualises networks of related samples using cytoscape.js

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
"""
import networkx as nx
import datetime
import logging

class snvNetwork():
    """ build and output Cytoscape compatible  networks """
    def __init__(self, snv_threshold=None):
        """ makes a undirected graph of samples;
        makes edges if pairwise SNV <= snv_threshold SNV
        assigns edge weight of 1/(1+SNV_distance)
                
        Note that in all the below documentation, 'guid' refers to a node identified by a guid,
        and 'guids' to multiple such nodes.
        """
        self.G= nx.Graph()
        self.snv_threshold = snv_threshold
            
                                
    def is_mixed(self,guid):
        """ returns True if the is_mixed attribute is set True """
        try:
            if self.G.nodes[guid]['is_mixed']==True:
                return True
        except KeyError:
            # there's no 'is_mixed' attribute
            pass
        return False     
    
    def snv2weight(self, x):
        """ returns 1/(1+x) """
        x = float(x)
        return 1/(1+x)
    
    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging """
        raise ZeroDivisionError(token)
    
    def to_dict(self):
        """ converts snv_clustering object to a dictionary.  """
        return nx.json_graph.node_link_data(self.G)

    def set_mixed(self,guid):
        """ marks guid as being mixed.
        """
           
        # set is_mixed attribute
        self.G.nodes[guid]['is_mixed']=True

    def add_sample(self, starting_guid, guids=None, neighbours=[], **kwargs):
        """ adds a sample, guid, linked to neighbours.
        - guid should be a string
        - guids are all the samples which are in the cluster.  edges outside these will not be displayed
        If None, then all edges will be displayed.
        - neighbours should be a list of tuples (guid, SNV)
        - additional arguments are added as node properties (e.g. surname = 'Smith')
        """

        # create a list of guid - neighbour links,
        # suitable for importing into networkx,
        # from the input data

        self.G.add_node(starting_guid, **kwargs)
        for item in neighbours:
            if not len(item)==2:
                raise TypeError("Neighbours must be a list of tuples (guid, snv) but it is {0}".format(item))
            add_edge = False
            if guids is None:
                add_edge = True
            else:
                if starting_guid in guids and item[0] in guids:
                    add_edge = True
            if add_edge is True:
                self.G.add_edge(starting_guid,item[0],weight=self.snv2weight(item[1]), snv=item[1])

    def guids(self):
        """ returns a set of all guids in the graph """
        return set(self.G.nodes)


    def network2cytoscapejs(self, max_edges = 1e5):
        """ this function is used to convert networkx to Cytoscape.js JSON format
            used by the elements property http://js.cytoscape.org/#notation/elements-json
            returns dictionary {'success': (0 or 1),  'message': description of result
            'elements':{data usable by cytoscape.js as elements property on construction}}
            
            will not return any network with > max_edges"""
        
        # lookup snv
        w = nx.get_edge_attributes(self.G, 'weight')
        snv = nx.get_edge_attributes(self.G, 'snv')
        # load all nodes into elements array
        elements = []
        for node in self.G.nodes(data = True):
           
            dat = {'id':node[0], **node[1]}
            elements.append( {
                  'group': 'nodes', 
                  'data': dat
                }
            )
        nNodes = len(elements)
        if nNodes == 0:
            return {'elements':elements, 'success':1, 'message':'No nodes found','nNodes':0, 'nEdges':0}

        #load all edges to edges array
        edge_id = 0
        for edge in self.G.edges():
            snv_dist = snv[(edge[0],edge[1])]
            if snv_dist <= self.snv_threshold:
                edge_id +=1
                
                if edge_id > max_edges:
                    return {'elements':{}, 'message':'Not rendered; > {0} edges present'.format(max_edges), 'success':0}
                elements.append(
                 { 'group':'edges',
                         'data':{
                            'id':edge_id,
                            'source':edge[0],
                            'target':edge[1],
                            'weight':w[(edge[0],edge[1])],
                            'snv':snv_dist
                            }
                        }
                )
        return {'elements':elements, 'success':1, 'message':'Graph with {0} nodes and {1} edges <= {2} SNV'.format(nNodes, edge_id, self.snv_threshold),
                'nNodes':nNodes, 'nEdges':edge_id}
