import unittest
import os
import pandas as pd
from collections import Counter
import shutil
from tree.tree_utils import IQTree, ManipulateTree, DepictTree

iqtree_requiring_test = unittest.skipIf(
    os.getenv("IQTREE_DIR") is None,
    "Test requires IQTREE v 2; no IQTREE environment found",
)


@iqtree_requiring_test
class Test_iqTree_1(unittest.TestCase):
    """test the iqTree module"""

    def runTest(self):
        """initialise"""

        iqt = IQTree(genome_length=30000)

        # read fasta file
        inputfile = "testdata/fasta/iqtree_test1.fasta"
        with open(inputfile, "rt") as f:
            inputstring = f.read()
        cnt = Counter(["A", "A", "C", "C", "T", "T", "G"])

        with self.assertRaises(TypeError):
            iqt.build(
                inputstring, 2, "unitTest_tmp/iqtree_test"
            )  # cnt must be a dictionary

        # empty the target directory
        targetdir = "/tmp/iqtree_test"
        shutil.rmtree(targetdir)
        iqt.build(inputstring, cnt, targetdir)  # cnt must be a dictionary

        targetdir = "/tmp/iqtree_test"
        shutil.rmtree(targetdir)

        # check there is not data in the directory
        testfile = os.path.join(targetdir, "alignment.fa.iqtree")
        self.assertFalse(os.path.exists(testfile))

        iqt.build(
            inputstring, {"A": 200, "C": 200, "G": 200, "T": 200}, targetdir
        )  # cnt must be a dictionary

        # check that files are in the directory
        testfile = os.path.join(targetdir, "alignment.fa.iqtree")
        self.assertTrue(os.path.exists(testfile))


class Test_iqTree_2(unittest.TestCase):
    """test the iqTree module when there is no IQTREE_DIR environment variable"""

    def runTest(self):
        """initialise"""

        # remove if present
        try:
            del os.environ["IQTREE_DIR"]
        except KeyError:
            pass

        iqt = IQTree(genome_length=30000)

        # read fasta file
        inputfile = "testdata/fasta/iqtree_test1.fasta"
        with open(inputfile, "rt") as f:
            inputstring = f.read()

        res = iqt.build(
            inputstring, 2, "unitTest_tmp/iqtree_test"
        )  # cnt must be a dictionary
        self.assertIsNone(res)


class Test_ManipulateTree(unittest.TestCase):
    """tests the ManipulateTree class"""

    def runTest(self):

        with open("testdata/tree/smalltree.nwk", "rt") as f:
            treetxt = f.read()

        mt = ManipulateTree(treetxt)

        # d1 = mt.rtd()
        # print(mt.nodes())

        n1 = set([x for x in mt.nodes()])
        mt.reroot("1")
        n2 = set([x for x in mt.nodes()])
        self.assertEqual(n1 - n2, set([]))

        # d2 = mt.rtd()

        mt.remove_outgroup()
        n3 = set([x for x in mt.nodes()])
        self.assertEqual(n1 - n3, set(["1"]))

        # d3 = mt.rtd()
        # print(d1,d2,d3)


@unittest.skip(
    """causes core dump in gitlab actions environment. runs in standard environment. reason requires investigation.
       probably to do with QT5 not being installed in github actions environment"""
)
class Test_DepictTree_1a(unittest.TestCase):
    """tests the DepictTree class"""

    def runTest(self):

        with open("testdata/ete3/test18.nwk", "rt") as f:
            treetxt = f.read()
        with open("testdata/ete3/test18.json", "rt") as f:
            json_str = f.read()
        metadata = pd.read_json(json_str)

        mt = DepictTree(treetxt, metadata, output_method = 'render_graphics')
        mt.render("unitTest_tmp/test.png")  # svg output causes github actions to fail


class Test_DepictTree_1b(unittest.TestCase):
    """tests the DepictTree class"""

    def runTest(self):

        with open("testdata/ete3/test18.nwk", "rt") as f:
            treetxt = f.read()
        with open("testdata/ete3/test18.json", "rt") as f:
            json_str = f.read()
        metadata = pd.read_json(json_str)

        outputfile = "unitTest_tmp/test.pickle"
        if os.path.exists(outputfile):
            os.unlink(outputfile)

        mt = DepictTree(treetxt, metadata)
        mt.render("unitTest_tmp/test")  # svg output causes github actions to fail

        self.assertTrue(os.path.exists(outputfile))

class Test_DepictTree_2(unittest.TestCase):
    """tests the DepictTree class"""

    def runTest(self):

        with open("testdata/ete3/test18.nwk", "rt") as f:
            treetxt = f.read()

        mt = ManipulateTree(treetxt)
        mt.zero_small_branch_lengths()
        mt.generate_simplified_tree()
