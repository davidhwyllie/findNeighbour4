import networkx as nx
import datetime
import unittest
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
            if self.G.node[guid]['is_mixed']==True:
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
        self.G.node[guid]['is_mixed']=True

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

# unittests
class test_nv_init(unittest.TestCase):
    """ tests init method of snv_clustering """
    def runTest(self):
        """ tests init """
        
        # no parameters except SNV threshold
        nv =snvNetwork(snv_threshold=12)
        
        self.assertEqual(type(nv.G), nx.classes.graph.Graph)

class test_nv_to_dict(unittest.TestCase):
    """ tests serialisation to dict """
    def runTest(self):
        nv =snvNetwork(snv_threshold=12)
        nv.G.add_node('a', is_mixed=True)
        nv.G.add_node('b', is_mixed=False)
        d = nv.to_dict()
        self.assertTrue(d['directed']==False)

class test_set_mixed(unittest.TestCase):
    """ tests _is_mixed function """
    def runTest(self):
        # set up
        nv = snvNetwork(snv_threshold=12)
        nv.add_sample('c')
        self.assertFalse(nv.is_mixed('c'))
        nv.set_mixed('c')
        self.assertTrue(nv.is_mixed('c'))
         
class test_is_mixed(unittest.TestCase):
    """ tests _is_mixed function """
    def runTest(self):
        # set up
        snvc = snvNetwork(snv_threshold=12)
        snvc.G.add_node('a', is_mixed=True)
        snvc.G.add_node('b', is_mixed=False)
        snvc.G.add_node('c')
        
        self.assertTrue(snvc.is_mixed('a'))
        self.assertFalse(snvc.is_mixed('b'))
        self.assertFalse(snvc.is_mixed('c'))
         
class test_guids(unittest.TestCase):
    """ tests recovery of list of guids """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
        self.assertEqual(snvc.guids(), set([]))
            
        # add two samples
        snvc.add_sample('n1')      
        snvc.add_sample('n2')      
        snvc.add_sample('n3')      
        self.assertEqual(snvc.guids(),set(['n1','n2','n3']))

class test_add_guids1(unittest.TestCase):
    """ tests recovery of list of guids """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
        self.assertEqual(snvc.guids(),set([]))
            
        # add two samples
        snvc.add_sample('n1')      
        snvc.add_sample('n1')      
        self.assertEqual(snvc.guids(),set(['n1']))

class test_add_guids_2(unittest.TestCase):
    """ tests addition of guids """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
        self.assertEqual(snvc.guids(),set([]))
        guid = "3a10b14a-218a-49f4-9397-c3dc7ef73818"
        neighbours = [["68f914ea-b2b0-493a-ba56-155235dcec30", 17], ["29a7737f-2298-4e57-bb0e-8678c84bedd8", 13], ["17922571-c82f-440f-bda4-3ebfb20fa554", 10]]

        snvc.add_sample(guid, neighbours = neighbours)      
        self.assertEqual(len(snvc.guids()),4)

class test_network2cytoscapejs_1(unittest.TestCase):
    """ tests addition of guids """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
        self.assertEqual(snvc.guids(),set([]))
        guid = "3a10b14a-218a-49f4-9397-c3dc7ef73818"
        neighbours = [["68f914ea-b2b0-493a-ba56-155235dcec30", 17], ["29a7737f-2298-4e57-bb0e-8678c84bedd8", 13], ["17922571-c82f-440f-bda4-3ebfb20fa554", 10]]

        snvc.add_sample(guid, neighbours = neighbours)      
        res = snvc.network2cytoscapejs()
        self.assertEqual(set(res.keys()), set(['message','success','elements', 'nNodes','nEdges']))
        self.assertEqual(res['nNodes'],4)
        self.assertEqual(res['nEdges'],1)
 
class test_network2cytoscapejs_2(unittest.TestCase):
    """ tests addition of guids """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
        self.assertEqual(snvc.guids(), set([]))
        guid = "3a10b14a-218a-49f4-9397-c3dc7ef73818"
        neighbours = []

        snvc.add_sample(guid, neighbours = neighbours)      
        res = snvc.network2cytoscapejs()
        self.assertEqual(set(res.keys()), set(['message','success','elements', 'nNodes','nEdges']))
        self.assertEqual(res['nNodes'],1)
        self.assertEqual(res['nEdges'],0)

class test_network2cytoscapejs_3(unittest.TestCase):
    """ tests addition of guids """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=20)
        self.assertEqual(snvc.guids(),set([]))
        guid = "3a10b14a-218a-49f4-9397-c3dc7ef73818"
        neighbours = [["68f914ea-b2b0-493a-ba56-155235dcec30", 17], ["29a7737f-2298-4e57-bb0e-8678c84bedd8", 13], ["17922571-c82f-440f-bda4-3ebfb20fa554", 10]]

        snvc.add_sample(guid, neighbours = neighbours)      
        res = snvc.network2cytoscapejs()
        self.assertEqual(set(res.keys()), set(['message','success','elements', 'nNodes','nEdges']))
        self.assertEqual(res['nNodes'],4)
        self.assertEqual(res['nEdges'],3)

class test_network2cytoscapejs_4(unittest.TestCase):
    """ tests addition of guid annotation """
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
        self.assertEqual(snvc.guids(),set([]))
        guid = "3a10b14a-218a-49f4-9397-c3dc7ef73818"
        neighbours = [["17922571-c82f-440f-bda4-3ebfb20fa554", 10]]

        snvc.add_sample(guid, neighbours = neighbours, surname='Smith')      
        res = snvc.network2cytoscapejs()
        self.assertEqual(set(res.keys()), set(['message','success','elements', 'nNodes','nEdges']))
        self.assertEqual(res['nNodes'],2)
        self.assertEqual(res['nEdges'],1)
        for item in res['elements']:
            if item['group'] == 'nodes':
                if item['data']['id']==guid:
                    self.assertEqual(item['data']['surname'], 'Smith')
class test_Raise_error(unittest.TestCase):
    """ tests raise_error"""
    def runTest(self):
        snvc = snvNetwork(snv_threshold=12)
                      
        with self.assertRaises(ZeroDivisionError):
            snvc.raise_error("token")
        
        