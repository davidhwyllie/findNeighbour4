""" tests malinkage.py

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

import os
import unittest
import json
import uuid
import random
import networkit as nk
from findn.common_utils import ConfigManager
from findn.persistence import Persistence
from findn.cw_seqComparer import cw_seqComparer
from catwalk.pycw_client import CatWalk
from snpclusters.ma_linkage import (
    MixtureChecker,
    MixPOREMixtureChecker,
    MixtureAwareLinkageResult,
    MixtureAwareLinkage,
)

# unittests
UNITTEST_MONGOCONN: str = "mongodb://localhost"
UNITTEST_RDBMSCONN: str = "sqlite://"


class MixtureCheckerTest(MixtureChecker):
    """a class which implements a MixtureChecker which sets guids as mixed if they begin with a number"""

    def __init__(self):
        """does nothing"""
        self.info = "Test check"

    def is_mixed(self, guid):
        """does not check for mixtures"""
        return {"mix_check_method": self.info, "is_mixed": guid[0] in ["0", "1", "2"]}


class MockPersistence:
    """simulates the fnPersistence class which provides access to findNeighbour stored data;
    the objective is to allow testing of the Clustering class, which requires a Persistence object to be provided to allow it to access findNeighbour SNP distances.
    therefore, only methods relevant to SNP distance provision are provided.
    The MockPersistence class generates data resembling a collection of samples which are related, their SNP distances,
    These methods are:
    guids()    lists the names of sequences present
    guid2neighbours  links between guids - returns type 1 output

    it also supports an isMixed() method, which is used for simulating whether a sample is mixed or not.
    """

    def cluster_delete_legacy(self, name):
        """delete any legacy data in the mock persistence store"""
        pass
        return

    def __init__(self, n_guids: int):
        """starts up a MockPersistence object; generates a set of pairwise links compatible with n_guids being part of a set of similar sequences."""
        self.latest_version = 0
        self.latest_version_behaviour = "increment"

        self.node2name = {}
        self.name2node = {}
        self.name2clusterid = {}
        self.node2clusterid = {}
        self._guid2neighbours = {}
        self.store = {}
        self.g = nk.generators.ClusteredRandomGraphGenerator(
            n_guids, int(n_guids / 2), 1, 0
        ).generate()
        for x in self.g.iterNodes():
            new_guid = str(uuid.uuid4())
            self.node2name[x] = new_guid
            self.name2node[new_guid] = x
            self._guid2neighbours[new_guid] = []

        # determine connected components ('single linkage clusters')
        cc = nk.components.ConnectedComponents(self.g)
        cc.run()
        for clusterid, component in enumerate(cc.getComponents()):
            for node in component:
                self.name2clusterid[self.node2name[node]] = clusterid
                self.node2clusterid[node] = clusterid

        # build a dictionary containing edges (SNV distances) in format 1
        for (x, y) in self.g.iterEdges():
            guid1 = self.node2name[x]
            guid2 = self.node2name[y]
            snv = random.sample(range(7), 1)[0]  # distances drawn randomly from 0-6
            self._guid2neighbours[guid1].append([guid2, snv])
            self._guid2neighbours[guid2].append([guid1, snv])

    def cluster_latest_version(self, clustering_version):
        """returns fake version information; increments if latest_version_behaviour is 'increment'"""
        if self.latest_version_behaviour == "increment":
            self.latest_version += 1
        return self.latest_version

    def guids(self):
        """returns all guids (sequence identifiers) in the network"""
        return set(self.name2node.keys())

    def guids_valid(self):
        """returns all guids (sequence identifiers) in the network , which are assumed to be valid"""
        return set(self.name2node.keys())

    def guid2neighbours(self, guid: str, returned_format: int = 1) -> dict:
        """returns neighbours of a guid in type 1 format [[guid1, snvdist1], [guid2, snvdist2], ...]

        note: attempts to change the returned_format result in a NotImplementedError

        Parameters:
        guid: the identifier of the sequence whose neighbours are sought
        returned_type: a placeholder, only a value of 1 is accepted
        """
        if not returned_format == 1:
            raise NotImplementedError(
                "the MockPersistence.guid2neighbours() method always returns type 1 output, but type {0} was requested".format(
                    returned_format
                )
            )
        return {"guid": guid, "neighbours": self._guid2neighbours[guid]}

    def cluster_store(self, key, serialisation):
        """stores serialisation in a dictionary using key"""
        self.store[key] = serialisation
        return None

    def cluster_read(self, key):
        try:
            return self.store[key]
        except KeyError:
            return None


class Test_MP(unittest.TestCase):
    """tests the MockPersistance object"""

    def runTest(self):

        p = MockPersistence(n_guids=20)
        self.assertEqual(20, len(p.node2name))
        self.assertEqual(20, len(p.name2node))

        # get guids
        guids = p.guids()
        self.assertEqual(len(guids), 20)
        self.assertTrue(isinstance(guids, set))

        guid = min(guids)
        with self.assertRaises(NotImplementedError):
            res = p.guid2neighbours(guid, returned_format=2)

        res = p.guid2neighbours(guid, returned_format=1)
        self.assertTrue(isinstance(res, dict))

        to_store = {"one": 1, "two": 2}
        res = p.cluster_store("myKey", to_store)
        self.assertIsNone(res)
        res = p.cluster_read("NoKey")
        self.assertIsNone(res)
        res = p.cluster_read("myKey")

        self.assertEqual(to_store, res)


class Test_MAL_1(unittest.TestCase):
    """tests the MixtureLinkage startup"""

    def runTest(self):
        p = MockPersistence(n_guids=20)
        m = MixtureAwareLinkage(PERSIST=p)

        # test there are no guids in the object on instatiation
        guids = m.guids()
        self.assertEqual(len(guids), 0)
        self.assertTrue(isinstance(guids, set))

        with self.assertRaises(TypeError):
            m = MixtureAwareLinkage(
                PERSIST=p, MIXCHECK=p
            )  # MIXCHECK is not a mixture checker


class Test_MAL_2(unittest.TestCase):
    """tests the MixtureLinkage _filter_guid2neighbours_by_cutoff method"""

    def runTest(self):
        p = MockPersistence(n_guids=20)
        m = MixtureAwareLinkage(PERSIST=p)

        # test there are no guids in the object on instatiation
        input = [["a", 6], ["b", 5], ["c", 4]]
        output = m._filter_guid2neighbours_by_snpcutoff(input, 5)
        self.assertEqual(output, [["b", 5], ["c", 4]])


class Test_MAL_3(unittest.TestCase):
    """tests the MixtureLinkage _filter_guid2neighbours_by_targets method"""

    def runTest(self):
        p = MockPersistence(n_guids=20)
        m = MixtureAwareLinkage(PERSIST=p)

        # test there are no guids in the object on instatiation
        input = [["a", 6], ["b", 5], ["c", 4]]
        output = m._filter_guid2neighbours_by_targets(input, ["c"])
        self.assertEqual(output, [["c", 4]])


class Test_MAL_4(unittest.TestCase):
    """tests the MixtureLinkage add method"""

    def runTest(self):
        p = MockPersistence(n_guids=20)
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        change_id = m.add(guids_to_add)
        self.assertEqual(m.guids(), p.guids())  # check all were added
        self.assertEqual(
            m.g.numberOfEdges(), p.g.numberOfEdges()
        )  # same number of edges
        self.assertEqual(
            m.g.numberOfNodes(), p.g.numberOfNodes()
        )  # same number of nodes
        self.assertEqual(
            m.g.numberOfNodes() - 1, change_id
        )  # change_id is the nodeid of the last item inserted, which is zero indexed
        self.assertEqual(
            len(m._name2meta), p.g.numberOfNodes()
        )  # one meta data entry per guid


class Test_MAL_4b(unittest.TestCase):
    """tests the MixtureLinkage add_sample method"""

    def runTest(self):
        p = MockPersistence(n_guids=20)
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        for guid in guids_to_add:
            change_id = m.add_sample(guid)
        self.assertEqual(m.guids(), p.guids())  # check all were added
        self.assertEqual(
            m.g.numberOfEdges(), p.g.numberOfEdges()
        )  # same number of edges
        self.assertEqual(
            m.g.numberOfNodes(), p.g.numberOfNodes()
        )  # same number of nodes
        self.assertEqual(
            m.g.numberOfNodes() - 1, change_id
        )  # change_id is the nodeid of the last item inserted, which is zero indexed
        self.assertEqual(
            len(m._name2meta), p.g.numberOfNodes()
        )  # one meta data entry per guid


class Test_MAL_5(unittest.TestCase):
    """tests the _connectedComponents() method"""

    def runTest(self):
        p = MockPersistence(n_guids=20)
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        m.add(guids_to_add)

        result = m._connectedComponents(what="node_id")
        for key in result.keys():  # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], int))

        result = m._connectedComponents(what="name")
        for key in result.keys():  # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], str))


class Test_MAL_7(unittest.TestCase):
    """tests the _connectedComponents() method, at scale, with check of whether it got the right answer"""

    def runTest(self):

        n_guids = 20
        p = MockPersistence(n_guids=n_guids)

        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        m.add(guids_to_add)

        # we compute the number of clusters in the original (as defined for node by p.node2clusterid[node])

        # for each cluster in the test (defined by the keys of result)
        resultkey2original = dict()

        result = m._connectedComponents(what="name")
        for key in result.keys():  # keys are hashes of contents
            resultkey2original[key] = set()

        for key in result.keys():  # keys are hashes of contents
            self.assertTrue(
                isinstance(result[key][0], str)
            )  # contents are lists of integers

            for i, name in enumerate(result[key]):  # for each of our clusters
                resultkey2original[key].add(
                    p.name2clusterid[name]
                )  # we expect exactly one original cluster: i.e. the clustering works

        for key in resultkey2original.keys():
            self.assertEqual(len(resultkey2original[key]), 1)

        # test guid conversion works
        result = m._connectedComponents(what="node_id")
        for key in result.keys():  # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], int))


class Test_MAL_8(unittest.TestCase):
    """tests the update() method, which adds any guids not in the MixtureAwareLinkage object to it."""

    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())
        first_five = set(guids_to_add[0:5])
        m.add(first_five)

        self.assertEqual(len(m.guids()), 5)
        self.assertEqual(m.guids(), first_five)

        m.update()
        self.assertEqual(len(m.guids()), n_guids)

        # check update adds all guids if none present
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        self.assertEqual(len(m.guids()), 0)
        m.update()
        self.assertEqual(len(m.guids()), n_guids)


class Test_MAL_9(unittest.TestCase):
    """tests the simplify() and ensure_edges() methods, which respectively
    rewires a graph to minimise edges while preserving connected components.
    restore all edges"""

    def runTest(self):

        n_guids = 2000
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())
        m.update()
        self.assertEqual(len(m.guids()), n_guids)

        # ensure the connected components remain the same after simplication
        pre = m._connectedComponents(what="node_id")
        edges_pre = m.g.numberOfEdges()

        self.assertFalse(m.is_simplified)
        m.simplify()
        self.assertTrue(m.is_simplified)

        post = m._connectedComponents(what="node_id")
        edges_post = m.g.numberOfEdges()

        self.assertEqual(pre, post)
        self.assertTrue(edges_pre > edges_post)

        m.ensure_edges(set(guids_to_add[0:10]))  # just a small number
        self.assertTrue(m.is_simplified)  # still simplified

        m.ensure_edges(set(guids_to_add))
        edges_restored = m.g.numberOfEdges()
        self.assertTrue(edges_restored > edges_post)
        self.assertEqual(edges_restored, edges_pre)
        self.assertFalse(m.is_simplified)


class Test_MAL_10(unittest.TestCase):
    """tests the cluster() method, which persists transformation of the data generated by _connectedComponents()"""

    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        # guids_to_add = list(p.guids())
        m.update()
        self.assertEqual(len(m.guids()), n_guids)

        self.assertTrue(m.clustered_at_change_id is None)
        self.assertEqual(0, len(m.name2cluster.keys()))
        self.assertEqual(0, len(m.cluster2names.keys()))

        # ensure the connected components remain the same after simplication
        m.cluster()

        self.assertEqual(n_guids - 1, m.clustered_at_change_id)
        self.assertEqual(n_guids, len(m.name2cluster.keys()))
        self.assertEqual(
            len(set([x for x in p.node2clusterid.values()])),
            len(m.cluster2names.keys()),
        )


class Test_MAL_11(unittest.TestCase):
    """tests the serialise() method, which returns a representation of the MixtureAwareLinkage object as a json string"""

    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        # guids_to_add = list(p.guids())
        m.update()

        retVal = m.serialise()
        self.assertEqual(
            set(retVal.keys()),
            set(
                [
                    "clusterid2clusterlabel",
                    "parameters",
                    "is_simplified",
                    "name",
                    "snv_threshold",
                    "mixed_sample_management",
                    "_edges",
                    "_node2name",
                    "_name2meta",
                ]
            ),
        )

        retVal2 = m.to_dict()  # synonym
        self.assertEqual(retVal, retVal2)

        # test jsonification round trip
        json_repr = json.dumps(retVal)
        retVal2 = json.loads(json_repr)

        # check we recover the graph
        m._deserialise_from_dict(retVal2)

        self.assertEqual(m.g.numberOfNodes(), p.g.numberOfNodes())
        self.assertEqual(m.g.numberOfEdges(), p.g.numberOfEdges())

        # check edges are all inplace
        for u, v in m.g.iterEdges():
            gu = m._node2name[u]
            gv = m._node2name[v]

            ou = p.name2node[gu]
            ov = p.name2node[gv]
            self.assertTrue(
                p.g.hasEdge(ou, ov)
            )  # check we wound up with what we started with


class Test_MAL_12(unittest.TestCase):
    """tests recreation of MixtureAwareLinkage from serialisation"""

    def runTest(self):

        n_guids = 50
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        # print("Creating",datetime.datetime.now())
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        # guids_to_add = list(p.guids())
        m.update()
        # print("Done",datetime.datetime.now())

        retVal = m.serialise()
        self.assertEqual(
            set(retVal.keys()),
            set(
                [
                    "clusterid2clusterlabel",
                    "snv_threshold",
                    "name",
                    "parameters",
                    "is_simplified",
                    "_edges",
                    "mixed_sample_management",
                    "_node2name",
                    "_name2meta",
                ]
            ),
        )

        # test jsonification round trip
        json_repr = json.dumps(retVal)

        retVal2 = json.loads(json_repr)

        # check we recover the graph
        # print("Recovering",datetime.datetime.now())
        m2 = MixtureAwareLinkage(PERSIST=p, serialisation=retVal2)
        # print("Done",datetime.datetime.now())

        self.assertEqual(m2.g.numberOfNodes(), p.g.numberOfNodes())
        self.assertEqual(m2.g.numberOfEdges(), p.g.numberOfEdges())

        # check edges are all inplace
        for u, v in m2.g.iterEdges():
            gu = m2._node2name[u]
            gv = m2._node2name[v]

            ou = p.name2node[gu]
            ov = p.name2node[gv]
            self.assertTrue(
                p.g.hasEdge(ou, ov)
            )  # check we wound up with what we started with


class Test_MAL_13(unittest.TestCase):
    """tests centrality measurement"""

    def runTest(self):

        n_guids = 50
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        # print("Creating",datetime.datetime.now())
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        # guids_to_add = list(p.guids())
        m.update()

        retVal = m.centrality()
        self.assertEqual(len(retVal.keys()), len(p.guids()))  # one for each guid


class Test_MAL_14a(unittest.TestCase):
    """tests the remove_edges() method, which
    removes the edges of one or more guids"""

    def runTest(self):

        n_guids = 2000
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p, snv_threshold=20, mixed_sample_management="include"
        )  # with this setting, mixtures are evaluated and edge removals occur
        guids_to_add = list(p.guids())
        m.update()
        self.assertEqual(len(m.guids()), n_guids)

        # the graph is not simplified.
        self.assertFalse(m.is_simplified)

        # remove the edges of  200 guids
        edges_pre = m.g.numberOfEdges()
        guids_to_remove = guids_to_add[0:200]

        # check it works
        m._connectedComponents(what="node_id")
        m.remove_edges(set(guids_to_remove))
        edges_post = m.g.numberOfEdges()
        self.assertTrue(edges_post < edges_pre)

        # simplify and check again
        self.assertFalse(m.is_simplified)
        m.simplify()
        self.assertTrue(m.is_simplified)

        # remove the edges of  200 guids
        edges_pre = m.g.numberOfEdges()
        guids_to_remove = guids_to_add[200:400]

        # check it works post simplification
        m._connectedComponents(what="node_id")
        m.remove_edges(set(guids_to_remove))
        edges_post = m.g.numberOfEdges()


class Test_MAL_14b(unittest.TestCase):
    """tests the add() method, which
    removes the edges of one or more guids if they are mixed and mixed_sample_management is set"""

    def runTest(self):

        n_guids = 2000
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m1 = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
        )  # with this setting, mixtures are  not evaluated and edge removals don't occur
        # guids_to_add = list(p.guids())
        m1.update()
        self.assertEqual(len(m1.guids()), n_guids)

        # the graph is not simplified.
        self.assertFalse(m1.is_simplified)

        # check update adds remaining guids
        m2 = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="ignore",
        )  # with this setting, mixtures are  not evaluated and edge removals don't occur
        # guids_to_add = list(p.guids())
        m2.update()
        self.assertEqual(len(m2.guids()), n_guids)

        # the graph is not simplified.
        self.assertFalse(m2.is_simplified)

        self.assertTrue(m1.g.numberOfEdges() < m2.g.numberOfEdges())


class Test_MAL_15(unittest.TestCase):
    """tests the name2meta method, which exports metadata about the samples"""

    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        # guids_to_add = list(p.guids())
        m.update()
        self.assertEqual(len(m.guids()), n_guids)
        self.assertEqual(len(m.name2meta()), n_guids)

        # check output after clustering
        m.cluster()

        self.assertEqual(len(m.name2meta()), n_guids)


class Test_TMT_1(unittest.TestCase):
    """tests the setting of mixed samples by the test mixture checker"""

    def runTest(self):

        MIXCHECK = MixtureCheckerTest()
        res = MIXCHECK.is_mixed("a12")
        self.assertFalse(res["is_mixed"])  # doesn't begin with a number
        res = MIXCHECK.is_mixed("1a2")
        self.assertTrue(res["is_mixed"])  # does begin with a number


class Test_MAL_16(unittest.TestCase):
    """tests the setting of mixed samples in 'include' mode"""

    def runTest(self):

        n_guids = 1000  # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
        # this is very much a worst case scenario:
        # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
        )
        # guids_to_add = list(p.guids())
        m.update()

        self.assertEqual(len(m.guids()), n_guids)
        self.assertEqual(len(m.name2meta()), n_guids)

        # check output after clustering
        m.cluster()

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta()

        for guid in res.index:
            begins_with_number = guid[0] in ["0", "1", "2"]
            has_neighbours = res.loc[guid, "nneighbours"] > 0
            if begins_with_number and has_neighbours:
                self.assertTrue(res.at[guid, "is_mixed"])
            elif not begins_with_number and has_neighbours:
                self.assertFalse(res.at[guid, "is_mixed"])
            else:
                self.assertIsNone(res.at[guid, "is_mixed"])

        cl = m.cluster()
        self.assertTrue(cl is not None)
        # now: the guids which begin with 0-2 should be mixed.
        # the others are not mixed.
        # the unmixed samples should all be in 1 cluster.
        in_clusters = set()
        for guid in res.index:
            begins_with_number = guid[0] in ["0", "1", "2"]
            has_neighbours = res.loc[guid, "nneighbours"] > 0
            if not begins_with_number:
                self.assertTrue(len(m.name2cluster[guid]["cluster_id"]), 1)  # not mixed
            else:
                in_clusters.add(len(m.name2cluster[guid]["cluster_id"]))

        #    all the samples should be clustered, either on their own or separately.
        self.assertTrue(min(in_clusters) > 0)


class Test_MAL_17(unittest.TestCase):
    """tests the setting of mixed samples in 'ignore' mode"""

    def runTest(self):

        n_guids = 1000  # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
        # this is very much a worst case scenario:
        # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="ignore",
        )
        # guids_to_add = list(p.guids())
        m.update()

        self.assertEqual(len(m.guids()), n_guids)
        self.assertEqual(len(m.name2meta()), n_guids)

        # check output after clustering
        m.cluster()

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta()

        # it does get scored in ignore mode
        for guid in res.index:
            begins_with_number = guid[0] in ["0", "1", "2"]
            has_neighbours = res.loc[guid, "nneighbours"] > 0
            if begins_with_number and has_neighbours:
                self.assertTrue(res.at[guid, "is_mixed"])
            elif not begins_with_number and has_neighbours:
                self.assertFalse(res.at[guid, "is_mixed"])
            else:
                self.assertIsNone(res.at[guid, "is_mixed"])

        cl = m.cluster()
        self.assertTrue(cl is not None)
        # now: the guids which begin with 0-2 should be ignored and treated like all the others.
        # the others are not mixed.
        # the unmixed samples should all be in 1 cluster.

        for guid in res.index:
            begins_with_number = guid[0] in ["0", "1", "2"]
            has_neighbours = res.loc[guid, "nneighbours"] > 0
            if not begins_with_number:
                self.assertTrue(len(m.name2cluster[guid]["cluster_id"]), 1)
            else:
                self.assertTrue(
                    len(m.name2cluster[guid]["cluster_id"]), 1
                )  # each sample is in its own cluster; no cross cluster samples.


class Test_MAL_18(unittest.TestCase):
    """tests the setting of mixed samples in 'exclude' mode"""

    def runTest(self):

        n_guids = 1000  # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
        # this is very much a worst case scenario:
        # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="exclude",
        )
        # guids_to_add = list(p.guids())
        m.update()

        self.assertEqual(len(m.guids()), n_guids)
        self.assertEqual(len(m.name2meta()), n_guids)

        # check output after clustering
        m.cluster()

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta()

        for guid in res.index:
            begins_with_number = guid[0] in ["0", "1", "2"]
            has_neighbours = res.loc[guid, "nneighbours"] > 0
            if begins_with_number and has_neighbours:
                self.assertTrue(res.at[guid, "is_mixed"])
                self.assertEqual(res.at[guid, "is_mixed"], m.is_mixed(guid))
            elif not begins_with_number and has_neighbours:
                self.assertFalse(res.at[guid, "is_mixed"])
            else:
                self.assertIsNone(res.at[guid, "is_mixed"])
        cl = m.cluster()
        self.assertTrue(cl is not None)

        # now: the guids which begin with 0-2 should be mixed.
        # the others are not mixed.
        # all samples should be in 1 cluster.

        for guid in res.index:
            begins_with_number = guid[0] in ["0", "1", "2"]
            has_neighbours = res.loc[guid, "nneighbours"] > 0
            if not begins_with_number:
                self.assertTrue(len(m.name2cluster[guid]["cluster_id"]), 1)
            else:
                self.assertTrue(len(m.name2cluster[guid]["cluster_id"]), 1)


class Test_MAL_19(unittest.TestCase):
    """tests output functions"""

    def runTest(self):

        n_guids = 10  # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
        # this is very much a worst case scenario:
        # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
            parameters={
                "Param1": 1,
                "Param2": 2,
                "snv_threshold": 12,
                "uncertain_base_type": "M",
            },
            name="MAL_19",
        )
        # guids_to_add = list(p.guids())
        m.update()

        self.assertEqual(len(m.guids()), n_guids)
        self.assertEqual(len(m.name2meta()), n_guids)

        # check output after clustering
        m.cluster()

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta()
        self.assertEqual(len(res.index), len(p.guids()))

        res = m.serialise_output()
        self.assertTrue(res["parameters"] is not None)
        self.assertTrue(res["guid2clustermeta"] is not None)


class Test_MAL_20(unittest.TestCase):
    """tests persistance and recovery functions"""

    def runTest(self):

        n_guids = 10  # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
        # this is very much a worst case scenario:
        # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids=n_guids)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
            parameters={
                "Param1": 1,
                "Param2": 2,
                "snv_threshold": 12,
                "uncertain_base_type": "M",
            },
            name="MAL_20",
        )
        # guids_to_add = list(p.guids())
        m.update()

        # check output after clustering
        m.cluster()

        m.persist(what="graph")
        m.persist(what="output")
        m.remove_legacy()

        # reload from persistence
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
            parameters={
                "Param1": 1,
                "Param2": 2,
                "snv_threshold": 12,
                "uncertain_base_type": "M",
            },
            name="MAL_20",
        )
        self.assertEqual(m.guids(), p.guids())

        self.assertEqual(m.guid2cluster_labels(), {})  # no labels assigned

        # apply labels; make some up
        cl2label = {}
        expected_labels = list()
        for cl in m.cluster2names:
            if len(m.cluster2names[cl]["guids"]) > 1:
                new_label = "LABEL-{0}".format(cl)
                cl2label[cl] = {"cluster_label": new_label}
                expected_labels.append(new_label)

        # apply the labels
        m._clusterid2clusterlabel = cl2label

        # check they are used correctly

        for guid in m.guid2cluster_labels():
            self.assertEqual(
                "LABEL-{0}".format(m.name2cluster[guid]["cluster_id"][0]),
                m.guid2cluster_labels()[guid][0],
            )
        self.assertEqual(set(m.existing_labels()), set(expected_labels))
        # check that the labels get persisted
        m.persist(what="graph")
        m.persist(what="output")
        m.remove_legacy()

        # reload from persistence
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
            parameters={
                "Param1": 1,
                "Param2": 2,
                "snv_threshold": 12,
                "uncertain_base_type": "M",
            },
            name="MAL_20",
        )

        # test that the guid2cluster maps are the same
        self.assertEqual(m._clusterid2clusterlabel, cl2label)

        # test apply_cluster_labels()
        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
            parameters={
                "Param1": 1,
                "Param2": 2,
                "snv_threshold": 12,
                "uncertain_base_type": "M",
            },
            name="MAL_20_v2",
        )
        # guids_to_add = list(p.guids())
        m.update()

        # check output after clustering
        m.cluster()

        # there are no labels
        self.assertEqual(m.guid2cluster_labels(), {})  # no labels assigned

        # apply labels; make some up
        cl2label = {}
        for cl in m.cluster2names:
            if len(m.cluster2names[cl]["guids"]) > 1:
                cl2label[cl] = {"cluster_label": "LABEL-{0}".format(cl)}

        # apply the labels
        m.apply_cluster_labels(cl2label)

        # check they are used correctly
        for guid in m.guid2cluster_labels():
            # print(guid, m.name2cluster[guid]['cluster_id'][0], m.guid2cluster_labels()[guid])
            self.assertEqual(
                "LABEL-{0}".format(m.name2cluster[guid]["cluster_id"][0]),
                m.guid2cluster_labels()[guid][0],
            )


class test_Raise_error(unittest.TestCase):
    """tests raise_error"""

    def runTest(self):
        p = MockPersistence(n_guids=50)

        # check update adds remaining guids
        m = MixtureAwareLinkage(
            PERSIST=p, MIXCHECK=MixtureCheckerTest(), snv_threshold=20
        )

        with self.assertRaises(ZeroDivisionError):
            m.raise_error("token")


class Test_MALR_1(unittest.TestCase):
    """tests MixtureAwareLinkageResults"""

    def runTest(self):

        # generate something to analyse
        n_guids = 10
        p = MockPersistence(n_guids=n_guids)
        m = MixtureAwareLinkage(
            PERSIST=p,
            MIXCHECK=MixtureCheckerTest(),
            snv_threshold=20,
            mixed_sample_management="include",
            parameters={
                "Param1": 1,
                "Param2": 2,
                "snv_threshold": 12,
                "uncertain_base_type": "M",
            },
            name="test",
        )
        # guids_to_add = list(p.guids())
        m.update()
        m.cluster()
        serialisation = m.serialise_output()

        # create and test all methods
        malr = MixtureAwareLinkageResult(serialisation)
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            for item in res:
                keys = set(item.keys())

                self.assertEqual(
                    keys,
                    set(
                        [
                            "clustering_algorithm",
                            "guid",
                            "cluster_id",
                            "cluster_label",
                            "cluster_loadtime",
                        ]
                    ),
                )
            self.assertIsNotNone(res)  # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))
        self.assertTrue(isinstance(malr.snv_threshold, int))

        m.persist(what="output")

        malr = MixtureAwareLinkageResult(PERSIST=p, name="test")
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            self.assertIsNotNone(res)  # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))

        # load an empty or missing data set
        malr = MixtureAwareLinkageResult(PERSIST=p, name="nodata")
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), set())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        self.assertTrue(malr.snv_threshold is None)

        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), 0)

        # load an empty or missing data set
        p.latest_version_behaviour = "nochange"  # don't increment version

        malr = MixtureAwareLinkageResult(PERSIST=p, name="test")
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            self.assertIsNotNone(res)  # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))
        t1 = malr.current_version_load_time
        malr.refresh()
        t2 = malr.current_version_load_time  # should not reload
        self.assertEqual(t1, t2)

        # load an empty or missing data set
        p.latest_version_behaviour = "increment"  # increment version

        malr = MixtureAwareLinkageResult(PERSIST=p, name="test")
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            self.assertIsNotNone(res)  # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))
        t1 = malr.current_version_load_time
        malr.refresh()
        t2 = malr.current_version_load_time
        self.assertNotEqual(t1, t2)  # updated


## database dependent tests

rdbms_test = unittest.skipIf(os.environ.get("NO_RDBMS_TESTS", False), "No rdbms tests")


@rdbms_test
class test_MIXCHECK_1_rdbms(unittest.TestCase):
    """tests mixpore mixture checker"""

    def runTest(self):

        rc = ConfigManager(os.path.join("config", "default_test_config.json"))
        CONFIG = rc.read_config()

        # get a clustering object's settings
        for clustering_name in CONFIG["CLUSTERING"].keys():
            clustering_setting = CONFIG["CLUSTERING"][clustering_name]

            pm = Persistence()
            PERSIST = pm.get_storage_object(
                dbname=CONFIG["SERVERNAME"],
                connString=UNITTEST_RDBMSCONN,
                debug=CONFIG["DEBUGMODE"],
                verbose=True,
            )
            PERSIST._delete_existing_data()

            # empty any test catwalk server
            cw = CatWalk(
                cw_binary_filepath=None,
                reference_name="H37RV",
                reference_filepath="reference/TB-ref.fasta",
                mask_filepath="reference/TB-exclude-adaptive.txt",
                max_n_positions=130000,
                bind_host="localhost",
                bind_port=5998,
            )

            # stop the server if it is running
            cw.stop()
            self.assertFalse(cw.server_is_running())

            hc = cw_seqComparer(
                reference=CONFIG["reference"],
                maxNs=CONFIG["MAXN_STORAGE"],
                snpCeiling=CONFIG["SNPCEILING"],
                excludePositions=CONFIG["excludePositions"],
                preComparer_parameters=CONFIG["PRECOMPARER_PARAMETERS"],
                PERSIST=PERSIST,
                unittesting=True,
            )

            mpmc = MixPOREMixtureChecker(hc, **clustering_setting)

            # check update adds remaining guids
            m = MixtureAwareLinkage(
                PERSIST=PERSIST,
                MIXCHECK=mpmc,
                mixed_sample_management=clustering_setting["mixed_sample_management"],
                snv_threshold=clustering_setting["snv_threshold"],
            )

            # add fake data
            guids_inserted = list()
            for i in range(1, 4):
                # print("Inserting",i)

                seq = list(str(CONFIG["reference"]))

                if i % 3 == 0:  # every third sample is mixed
                    is_mixed = True
                    guid_to_insert = "mixed_{0}".format(i)
                else:
                    is_mixed = False
                    guid_to_insert = "nomix_{0}".format(i)

                # make 5 mutations at position 500,000

                offset = 700000
                for j in range(i + 5):
                    mutbase = offset + j
                    ref = seq[mutbase]
                    if is_mixed is False:
                        if i % 2 == 0:
                            if not ref == "T":
                                seq[mutbase] = "T"
                            else:
                                seq[mutbase] = "A"
                        elif i % 2 == 1:
                            if not ref == "C":
                                seq[mutbase] = "C"
                            else:
                                seq[mutbase] = "G"
                    else:
                        seq[mutbase] = "M"
                seq = "".join(seq)
                # print(i,guid_to_insert, seq[699995:700020])
                guids_inserted.append(guid_to_insert)
                obj = hc.compress(seq)
                loginfo = hc.persist(obj, guid_to_insert)
                # print(loginfo)
                self.assertTrue(loginfo is not None)
            # test the mixporemixture checker.
            for guid in guids_inserted:
                res = mpmc.is_mixed(guid)
                # print(guid, res)
                self.assertEqual(
                    "mixed" in guid, res["is_mixed"]
                )  # should identify all mixed guids

            m.update()
            m.cluster()

            res = m.name2meta()
            # print(res)
            for guid in res.index:

                self.assertEqual(
                    "mixed" in guid, res.at[guid, "is_mixed"]
                )  # should identify all mixed guids

                # check whether it's what we expect
                if clustering_setting["mixed_sample_management"] in [
                    "ignore",
                    "exclude",
                ]:
                    self.assertEqual(
                        len(res.at[guid, "cluster_id"]), 1
                    )  # samples are in one cluster
                if clustering_setting["mixed_sample_management"] == "include":
                    if "mixed" in guid:
                        self.assertTrue(len(res.at[guid, "cluster_id"]) > 0)
                    else:
                        self.assertEqual(len(res.at[guid, "cluster_id"]), 1)

            PERSIST.closedown()


# skip these tests if the NO_MONGO_TESTS variable exists
mongo_test = unittest.skipIf(
    os.environ.get("NO_MONGO_TESTS", False), "no mongo tests performed"
)


@mongo_test
class test_MIXCHECK_1_mongo(unittest.TestCase):
    """tests mixpore mixture checker"""

    def runTest(self):

        rc = ConfigManager(os.path.join("config", "default_test_config.json"))
        CONFIG = rc.read_config()

        # get a clustering object's settings
        for clustering_name in CONFIG["CLUSTERING"].keys():
            clustering_setting = CONFIG["CLUSTERING"][clustering_name]

            pm = Persistence()
            PERSIST = pm.get_storage_object(
                dbname=CONFIG["SERVERNAME"],
                connString=UNITTEST_MONGOCONN,
                debug=CONFIG["DEBUGMODE"],
                verbose=True,
            )
            PERSIST._delete_existing_data()

            # empty any test catwalk server
            cw = CatWalk(
                cw_binary_filepath=None,
                reference_name="H37RV",
                reference_filepath="reference/TB-ref.fasta",
                mask_filepath="reference/TB-exclude-adaptive.txt",
                max_n_positions=130000,
                bind_host="localhost",
                bind_port=5998,
            )

            # stop the server if it is running
            cw.stop()
            self.assertFalse(cw.server_is_running())

            hc = cw_seqComparer(
                reference=CONFIG["reference"],
                maxNs=CONFIG["MAXN_STORAGE"],
                snpCeiling=CONFIG["SNPCEILING"],
                excludePositions=CONFIG["excludePositions"],
                preComparer_parameters=CONFIG["PRECOMPARER_PARAMETERS"],
                PERSIST=PERSIST,
                unittesting=True,
            )

            mpmc = MixPOREMixtureChecker(hc, **clustering_setting)

            # check update adds remaining guids
            m = MixtureAwareLinkage(
                PERSIST=PERSIST,
                MIXCHECK=mpmc,
                mixed_sample_management=clustering_setting["mixed_sample_management"],
                snv_threshold=clustering_setting["snv_threshold"],
            )

            # add fake data
            guids_inserted = list()
            for i in range(1, 4):
                # print("Inserting",i)

                seq = list(str(CONFIG["reference"]))

                if i % 3 == 0:  # every third sample is mixed
                    is_mixed = True
                    guid_to_insert = "mixed_{0}".format(i)
                else:
                    is_mixed = False
                    guid_to_insert = "nomix_{0}".format(i)

                # make 5 mutations at position 500,000

                offset = 700000
                for j in range(i + 5):
                    mutbase = offset + j
                    ref = seq[mutbase]
                    if is_mixed is False:
                        if i % 2 == 0:
                            if not ref == "T":
                                seq[mutbase] = "T"
                            else:
                                seq[mutbase] = "A"
                        elif i % 2 == 1:
                            if not ref == "C":
                                seq[mutbase] = "C"
                            else:
                                seq[mutbase] = "G"
                    else:
                        seq[mutbase] = "M"
                seq = "".join(seq)
                # print(i,guid_to_insert, seq[699995:700020])
                guids_inserted.append(guid_to_insert)
                obj = hc.compress(seq)
                loginfo = hc.persist(obj, guid_to_insert)
                # print(loginfo)
                self.assertTrue(loginfo is not None)
            # test the mixporemixture checker.
            for guid in guids_inserted:
                res = mpmc.is_mixed(guid)
                # print(guid, res)
                self.assertEqual(
                    "mixed" in guid, res["is_mixed"]
                )  # should identify all mixed guids

            m.update()
            m.cluster()

            res = m.name2meta()
            # print(res)
            for guid in res.index:

                self.assertEqual(
                    "mixed" in guid, res.at[guid, "is_mixed"]
                )  # should identify all mixed guids

                # check whether it's what we expect
                if clustering_setting["mixed_sample_management"] in [
                    "ignore",
                    "exclude",
                ]:
                    self.assertEqual(
                        len(res.at[guid, "cluster_id"]), 1
                    )  # samples are in one cluster
                if clustering_setting["mixed_sample_management"] == "include":
                    if "mixed" in guid:
                        self.assertTrue(len(res.at[guid, "cluster_id"]) > 0)
                    else:
                        self.assertEqual(len(res.at[guid, "cluster_id"]), 1)
            PERSIST.closedown()
