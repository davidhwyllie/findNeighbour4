""" tests nucleicacid.py

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

import unittest
import networkx as nx
from findn.visualiseNetwork import snvNetwork


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
        
        

