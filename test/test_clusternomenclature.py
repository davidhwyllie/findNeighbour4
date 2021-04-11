""" tests clusternomenclature.py

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
from snpclusters.clusternomenclature import ClusterNameAssigner, ClusterNomenclature

# unittests


class Test_ClusterNomenclature_1(unittest.TestCase):
    """ tests cluster nomenclature exposed functions """

    def runTest(self):
        # test whether input is checked properly
        with self.assertRaises(ValueError):
            n = ClusterNomenclature(cluster_nomenclature_method="NoMethod")
        with self.assertRaises(ValueError):
            n = ClusterNomenclature(deserialise_from=[], existing_labels=[])

        # make an integer ClusterNomenclature generating object
        n = ClusterNomenclature()
        self.assertEqual(n.cluster_nomenclature_method, "integer")
        self.assertEqual(n._labels, set([]))
        self.assertEqual(n.existing_labels(), [])
        self.assertEqual(n.serialise(), {"cluster_nomenclature_method": "integer", "labels": []})
        res = n.serialise()
        n = ClusterNomenclature(deserialise_from=res, existing_labels=None)
        res = n.new_label()
        self.assertEqual(res, 1)

        res = n.new_label()
        self.assertEqual(res, 2)

        # make an integer ClusterNomenclature generating object
        n = ClusterNomenclature(cluster_nomenclature_method="TB")
        self.assertEqual(n.cluster_nomenclature_method, "TB")
        self.assertEqual(n._labels, set([]))
        self.assertEqual(n.existing_labels(), [])
        self.assertEqual(n.serialise(), {"cluster_nomenclature_method": "TB", "labels": []})
        res = n.serialise()

        n = ClusterNomenclature(deserialise_from=res, existing_labels=None)
        self.assertEqual(n.cluster_nomenclature_method, "TB")
        self.assertEqual(n._labels, set([]))

        res = n.new_label()
        self.assertNotEqual(res, 1)
        self.assertEqual(res, "AA001")
        res = n.new_label()
        self.assertEqual(res, "AA002")
        res = n.new_label()
        self.assertEqual(res, "AA003")

        serialise_1 = n.serialise()
        n2 = ClusterNomenclature(deserialise_from=serialise_1, existing_labels=None)
        res = n2.new_label()
        self.assertEqual(res, "AA004")

        n = ClusterNomenclature(cluster_nomenclature_method="TB")
        res = n.new_label()
        self.assertEqual(res, "AA001")
        res = n.new_label(from_existing="AA001")
        self.assertEqual(res, "AA001v1")
        res = n.new_label(from_existing="AA001v1")
        self.assertEqual(res, "AA001v2")

        serialise_1 = n.serialise()
        n2 = ClusterNomenclature(deserialise_from=serialise_1, existing_labels=None)
        res = n2.new_label()
        self.assertEqual(res, "AA002")
        res = n2.new_label(from_existing="AA001")
        self.assertEqual(res, "AA001v3")


class Test_ClusterNomenclature_2(unittest.TestCase):
    """ tests internal cluster naming function """

    def runTest(self):

        # make an integer generating object
        n = ClusterNomenclature()
        res = n._encode_integer(1)
        self.assertEqual(res, "AA001")
        res = n._encode_integer(2)
        self.assertEqual(res, "AA002")

        with self.assertRaises(TypeError):
            n._decode_TB(0)
        with self.assertRaises(ValueError):
            n._decode_TB("AA")
        with self.assertRaises(ValueError):
            n._decode_TB("--123")
        with self.assertRaises(ValueError):
            n._decode_TB("AA---")

        for i in range(25000):
            res1 = n._encode_integer(i)
            res2 = n._decode_TB(res1)
            self.assertEqual(i, res2)


class Test_ClusterNameAssigner_1(unittest.TestCase):
    """ test cluster name assigner """

    def runTest(self):
        n = ClusterNomenclature(cluster_nomenclature_method="TB")
        cna = ClusterNameAssigner(n)
        previous_guid2cluster_label = {"a": ["AA001"], "b": ["AA001"], "e": ["AB001"], "f": ["AB001"]}
        clusterid2guid = {1: {"guids": ["a", "b", "c", "d"]}, 2: {"guids": ["e", "f"]}}
        retVal = cna.assign_new_clusternames(
            clusterid2guid=clusterid2guid, previous_guid2cluster_label=previous_guid2cluster_label
        )

        self.assertEqual(retVal, {1: {"cluster_label": "AA001"}, 2: {"cluster_label": "AB001"}})  # use previous names
        cna = ClusterNameAssigner(n)

        retVal = cna.assign_new_clusternames(clusterid2guid=clusterid2guid)
        self.assertEqual(retVal, {1: {"cluster_label": "AA001"}, 2: {"cluster_label": "AA002"}})  # generate new names

        cna = ClusterNameAssigner(n)
        previous_guid2cluster_label = {
            "a": ["AA001"],
            "b": ["AA001"],
            "e": ["AA001"],
            "f": ["AA001"],
            "g": ["AA002"],
            "h": ["AA002"],
        }
        clusterid2guid = {1: {"guids": ["a", "b", "c", "d"]}, 2: {"guids": ["e", "f"]}, 3: {"guids": ["g", "h"]}}

        retVal = cna.assign_new_clusternames(
            clusterid2guid=clusterid2guid, previous_guid2cluster_label=previous_guid2cluster_label
        )
        expected_result = {1: {"cluster_label": "AA001"}, 2: {"cluster_label": "AA001v1"}, 3: {"cluster_label": "AA002"}}
        self.assertEqual(
            retVal, expected_result
        )  # generate new names from previous ones, keeping one unchanged and without altering irrelevant clusters

        # on repetition, it should not change
        previous_guid2cluster_label = {
            "a": ["AA001"],
            "b": ["AA001"],
            "e": ["AA001v1"],
            "f": ["AA001v1"],
            "g": ["AA002"],
            "h": ["AA002"],
        }
        retVal = cna.assign_new_clusternames(
            clusterid2guid=clusterid2guid, previous_guid2cluster_label=previous_guid2cluster_label
        )
        expected_result = {1: {"cluster_label": "AA001"}, 2: {"cluster_label": "AA001v1"}, 3: {"cluster_label": "AA002"}}
        self.assertEqual(
            retVal, expected_result
        )  # generate new names from previous ones, keeping one unchanged and without altering irrelevant clusters
        retVal = cna.assign_new_clusternames(
            clusterid2guid=clusterid2guid, previous_guid2cluster_label=previous_guid2cluster_label
        )
        expected_result = {1: {"cluster_label": "AA001"}, 2: {"cluster_label": "AA001v1"}, 3: {"cluster_label": "AA002"}}
        self.assertEqual(
            retVal, expected_result
        )  # generate new names from previous ones, keeping one unchanged and without altering irrelevant clusters
