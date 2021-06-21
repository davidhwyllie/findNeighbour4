import unittest
from tree.manipulate_tree import ManipulateTree


class Test_ManipulateTree(unittest.TestCase):
    """tests the MSAStore class"""

    def runTest(self):

        with open("testdata/tree/smalltree.nwk", "rt") as f:
            treetxt = f.read()

            mt = ManipulateTree(treetxt)
            n1 = set(mt.nodes())
            print(n1)
            mt.reroot("1")
            n2 = set(mt.nodes())
            self.assertEqual(n1 - n2, set([]))
            print(len(n2))
            mt.remove_outgroup()
            self.assertEqual(n1 - n2, set(["1"]))
