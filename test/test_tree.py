import unittest
from tree.manipulate_tree import ManipulateTree

class Test_ManipulateTree(unittest.TestCase):
    """ tests the MSAStore class"""
    def runTest(self):

        with open("testdata/tree/bigtree.nwk", 'rt') as f:
            treetxt = f.read()
            
            mt = ManipulateTree()
            mt.reroot(treetxt, '--Wuhan-Reference--')



