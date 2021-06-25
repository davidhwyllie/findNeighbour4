import unittest
import os
from collections import Counter
import shutil
import pickle
from tree.tree_utils import IQTree, ManipulateTree, DepictTree

iqtree_requiring_test = unittest.skipIf(
    os.getenv('IQTREE_DIR') is None, 'Test requires IQTREE v 2; no IQTREE environment found'
)

@iqtree_requiring_test
class Test_iqTree_1(unittest.TestCase):
    """ test the iqTree module """

    def runTest(self):
        """ initialise """

        iqt = IQTree(genome_length = 30000)

        # read fasta file
        inputfile = '/data/software/fn4dev/testdata/fasta/iqtree_test1.fasta'
        with open(inputfile, 'rt') as f:
            inputstring = f.read()       
        cnt = Counter(['A', 'A', 'C', 'C', 'T', 'T', 'G'])

        with self.assertRaises(TypeError):
            iqt.build(inputstring, 2, '/tmp/iqtree_test')           # cnt must be a dictionary
        
        # empty the target directory
        targetdir= '/tmp/iqtree_test'
        shutil.rmtree(targetdir)
        iqt.build(inputstring, cnt, targetdir)           # cnt must be a dictionary

        targetdir= '/tmp/iqtree_test'
        shutil.rmtree(targetdir)

        # check there is not data in the directory
        testfile = os.path.join(targetdir, 'alignment.fa.iqtree')
        self.assertFalse(os.path.exists(testfile))

        iqt.build(inputstring, {'A': 200, 'C': 200, 'G': 200, 'T': 200}, targetdir)           # cnt must be a dictionary

        # check that files are in the directory
        testfile = os.path.join(targetdir, 'alignment.fa.iqtree')
        self.assertTrue(os.path.exists(testfile))

class Test_iqTree_2(unittest.TestCase):
    """ test the iqTree module when there is no IQTREE_DIR environment variable """

    def runTest(self):
        """ initialise """

        # remove if present
        try:
            del os.environ['IQTREE_DIR']
        except KeyError:
            pass

        iqt = IQTree(genome_length = 30000)

        # read fasta file
        inputfile = '/data/software/fn4dev/testdata/fasta/iqtree_test1.fasta'
        with open(inputfile, 'rt') as f:
            inputstring = f.read()       

        res = iqt.build(inputstring, 2, '/tmp/iqtree_test')           # cnt must be a dictionary
        self.assertIsNone(res)


class Test_ManipulateTree(unittest.TestCase):
    """tests the ManipulateTree class"""

    def runTest(self):

        with open("testdata/tree/smalltree.nwk", "rt") as f:
            treetxt = f.read()

        mt = ManipulateTree(treetxt)

        #d1 = mt.rtd()
        #print(mt.nodes())

        n1 = set([x for x in mt.nodes()])
        mt.reroot("1")
        n2 = set([x for x in mt.nodes()])
        self.assertEqual(n1 - n2, set([]))

        #d2 = mt.rtd()

        mt.remove_outgroup()
        n3 = set([x for x in mt.nodes()])
        self.assertEqual(n1 - n3, set(["1"]))

        #d3 = mt.rtd()
        #print(d1,d2,d3)

class Test_DepictTree_1(unittest.TestCase):
    """tests the DepictTree class"""

    def runTest(self):

        with open("testdata/ete3/test18.nwk", "rt") as f:
            treetxt = f.read()
        with open("testdata/ete3/test18.pickle", "rb") as f:
            metadata = pickle.load(f)
        
        mt = DepictTree(treetxt, metadata)
        mt.render('analysis/test.svg')

class Test_DepictTree_2(unittest.TestCase):
    """tests the DepictTree class"""

    def runTest(self):

        with open("testdata/ete3/test18.nwk", "rt") as f:
            treetxt = f.read()
        
        mt = ManipulateTree(treetxt)
        mt.zero_small_branch_lengths()
        mt.generate_simplified_tree()
