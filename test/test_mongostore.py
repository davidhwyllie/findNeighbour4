""" test mongo db connection object """
import unittest
import time
import json
import pymongo  # type: ignore
import os
import pandas as pd  # type: ignore
import pickle
import datetime
from findn.NucleicAcid import NucleicAcid

from findn.mongoStore import fn3persistence

## persistence unit tests
UNITTEST_MONGOCONN: str = "mongodb://localhost"

# skip these tests if the NO_MONGO_TESTS variable exists
mongo_test = unittest.skipIf(
    os.environ.get("NO_MONGO_TESTS", False), "no mongo tests performed"
)


@mongo_test
class Test_Server_Monitoring_0(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.server_monitoring_store(message="one")

        res = p.recent_server_monitoring(100)

        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))


@mongo_test
class Test_Server_Monitoring_1(unittest.TestCase):
    """tests recovery of database monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        res1 = p.recent_database_monitoring(100)

        db_summary = p.summarise_stored_items()
        p.server_monitoring_store(
            what="dbManager", message="Repacking", guid="-", content=db_summary
        )
        res2 = p.recent_database_monitoring(100)
        self.assertEqual(
            res1, {"latest_stats": {"storage_ratio": 1}, "recompression_data": False}
        )
        self.assertEqual(res2["recompression_data"], True)
        self.assertEqual(res2["latest_stats"]["storage_ratio"], 1)

        json.dumps(res2)  # should succeed


@mongo_test
class Test_Server_Monitoring_2(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.server_monitoring_store(message="one")
        p.server_monitoring_store(message="two")
        p.server_monitoring_store(message="three")

        res = p.recent_server_monitoring(0)
        self.assertEqual(len(res), 0)
        self.assertTrue(isinstance(res, list))

        res = p.recent_server_monitoring(1)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))

        res = p.recent_server_monitoring(3)
        self.assertEqual(len(res), 3)
        self.assertTrue(isinstance(res, list))

        res = p.recent_server_monitoring(5)
        self.assertEqual(len(res), 3)
        self.assertTrue(isinstance(res, list))

        with self.assertRaises(ValueError):
            res = p.recent_server_monitoring(-1)

        with self.assertRaises(TypeError):
            res = p.recent_server_monitoring("thing")  # type: ignore


@mongo_test
class Test_Server_Monitoring_3(unittest.TestCase):
    """checks whether server_monitoring_min_interval_msec control works"""

    def runTest(self):
        p = fn3persistence(
            connString=UNITTEST_MONGOCONN,
            debug=2,
            server_monitoring_min_interval_msec=1000,
        )  # no logging for within 1 secs of another event

        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 0)

        retVal = p.server_monitoring_store(message="one")  # should insert
        self.assertEqual(retVal, True)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))

        retVal = p.server_monitoring_store(message="two")  # should not insert
        self.assertEqual(retVal, False)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))

        time.sleep(3)  # seconds
        retVal = p.server_monitoring_store(message="three")  # should insert
        self.assertEqual(retVal, True)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 2)
        self.assertTrue(isinstance(res, list))


@mongo_test
class Test_Server_Monitoring_4(unittest.TestCase):
    """checks whether delete_server_monitoring_entries"""

    def runTest(self):
        p = fn3persistence(
            connString=UNITTEST_MONGOCONN,
            debug=2,
            server_monitoring_min_interval_msec=0,
        )
        retVal = p.server_monitoring_store(message="one")  # should insert
        self.assertEqual(retVal, True)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))
        p.delete_server_monitoring_entries(1)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))

        time.sleep(2)  # seconds

        p.delete_server_monitoring_entries(1)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 0)
        self.assertTrue(isinstance(res, list))


@mongo_test
class Test_SeqMeta_singleton(unittest.TestCase):
    """tests guid2neighboursOf"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 12},
                "guid2": {"dist": 0},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
                "guid5": {"dist": 5},
            },
        )
        res1 = p.guid2neighbours("srcguid", returned_format=1)

        self.assertEqual(5, len(res1["neighbours"]))
        singletons = p.singletons(method="exact")
        self.assertEqual(len(singletons.index), 5)
        singletons = p.singletons(method="approximate")
        self.assertEqual(len(singletons.index), 5)


@mongo_test
class Test_SeqMeta_guid2neighbour_8(unittest.TestCase):
    """tests guid2neighboursOf"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 12},
                "guid2": {"dist": 0},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
                "guid5": {"dist": 5},
            },
        )

        res1 = p.guid2neighbours("srcguid", returned_format=1)
        self.assertEqual(5, len(res1["neighbours"]))
        with self.assertRaises(NotImplementedError):
            res2 = p.guid2neighbours("srcguid", returned_format=2)
            self.assertTrue(res2 is not None)

        res3 = p.guid2neighbours("srcguid", returned_format=3)
        self.assertEqual(5, len(res3["neighbours"]))
        res4 = p.guid2neighbours("srcguid", returned_format=4)
        self.assertEqual(5, len(res4["neighbours"]))


@mongo_test
class Test_SeqMeta_guid2neighbour_7(unittest.TestCase):
    """tests guid2neighboursOf"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 12},
                "guid2": {"dist": 0},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
                "guid5": {"dist": 5},
            },
        )
        with self.assertRaises(NotImplementedError):
            res1 = p.guid2neighbours("srcguid", returned_format=2)
        res1 = p.guid2neighbours("srcguid", returned_format=1)
        self.assertEqual(5, len(res1["neighbours"]))
        p.guid2neighbour_repack(guid="srcguid")
        res2 = p.guid2neighbours("srcguid", returned_format=1)
        self.assertEqual(5, len(res2["neighbours"]))


@mongo_test
class Test_SeqMeta_guid2neighbour_5(unittest.TestCase):
    """tests repack"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid", {"guid1": {"dist": 12}, "guid2": {"dist": 0}}
        )

        # check the insert worked
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 1)

        p.guid2neighbour_repack(guid="srcguid")  # will make no difference

        # should make no change for 'srcguid'
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 1)

        # the src guid entry should contain two keys, 'guid1' and 'guid2' in its neighbours: section
        res = p.db.guid2neighbour.find_one({"guid": "srcguid"})
        self.assertEqual(set(res["neighbours"].keys()), set(["guid1", "guid2"]))

        # add some more links, creating 2 entries for guid1
        res = p.guid2neighbour_add_links("srcguid2", {"guid1": {"dist": 12}})

        # check the insert worked
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 2)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid2"})
        self.assertEqual(res, 1)

        p.guid2neighbour_repack(guid="srcguid2")  # will make no difference

        # should make no change for 'srcguid2'
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 2)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid2"})
        self.assertEqual(res, 1)

        p.guid2neighbour_repack(
            guid="guid1"
        )  # will compress guid1's two entries into one.

        # should make no change for 'srcguid1'
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid2"})
        self.assertEqual(res, 1)


@mongo_test
class Test_SeqMeta_guid2neighbour_4a(unittest.TestCase):
    """tests creation of new guid2neighbour entries"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid", {"guid1": {"dist": 12}, "guid2": {"dist": 0}}
        )
        self.assertEqual(p.max_neighbours_per_document, 3)  # debug setting
        a, b, c = p._audit_storage("srcguid")

        self.assertEqual(len(b.index), 2)  # 2 guids
        self.assertEqual(set(b["rstat"]), set(["m"]))  # in a multirecord

        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 1)


@mongo_test
class Test_SeqMeta_guid2neighbour_4b(unittest.TestCase):
    """tests creation of new guid2neighbour entries"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid",
            {"guid1": {"dist": 12}, "guid2": {"dist": 0}, "guid3": {"dist": 6}},
        )
        self.assertEqual(p.max_neighbours_per_document, 3)
        a, b, c = p._audit_storage("srcguid")

        self.assertEqual(len(b.index), 3)
        self.assertEqual(set(b["rstat"]), set(["f"]))  # one full, one singleton

        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid2"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid3"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 1)


@mongo_test
class Test_SeqMeta_guid2neighbour_4c(unittest.TestCase):
    """tests creation of new guid2neighbour entries"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 1},
                "guid2": {"dist": 2},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
            },
        )
        self.assertEqual(p.max_neighbours_per_document, 3)
        a, b, c = p._audit_storage("srcguid")

        self.assertEqual(len(b.index), 4)
        self.assertEqual(set(b["rstat"]), set(["f", "s"]))  # one full, one singleton

        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid2"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid3"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid4"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 2)


@mongo_test
class Test_SeqMeta_guid2neighbour_4d(unittest.TestCase):
    """tests creation of new guid2neighbour entries"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 1},
                "guid2": {"dist": 2},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
                "guid5": {"dist": 5},
            },
        )
        self.assertEqual(p.max_neighbours_per_document, 3)
        a, b, c = p._audit_storage("srcguid")

        self.assertEqual(len(b.index), 5)
        self.assertEqual(set(b["rstat"]), set(["f", "m"]))  # one full, one mixed

        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid2"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid3"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid4"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid5"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 2)


@mongo_test
class Test_SeqMeta_guid2neighbour_4e(unittest.TestCase):
    """tests creation of new guid2neighbour entries"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 1},
                "guid2": {"dist": 2},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
                "guid5": {"dist": 5},
                "guid6": {"dist": 6},
            },
        )
        self.assertEqual(p.max_neighbours_per_document, 3)
        a, b, c = p._audit_storage("srcguid")

        self.assertEqual(len(b.index), 6)
        self.assertEqual(set(b["rstat"]), set(["f"]))  # one full, one mixed

        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid2"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid3"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid4"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid5"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid6"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 2)


@mongo_test
class Test_SeqMeta_guid2neighbour_4f(unittest.TestCase):
    """tests creation of new guid2neighbour entries"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links(
            "srcguid",
            {
                "guid1": {"dist": 1},
                "guid2": {"dist": 2},
                "guid3": {"dist": 3},
                "guid4": {"dist": 4},
                "guid5": {"dist": 5},
                "guid6": {"dist": 6},
                "guid7": {"dist": 7},
            },
        )
        self.assertEqual(p.max_neighbours_per_document, 3)
        a, b, c = p._audit_storage("srcguid")

        self.assertEqual(len(b.index), 7)
        self.assertEqual(set(b["rstat"]), set(["f", "s"]))  # two full, one single

        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid2"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid3"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid4"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid5"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid6"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "guid7"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 3)


@mongo_test
class Test_SeqMeta_guid2neighbour_3(unittest.TestCase):
    """tests creation of a new guid2neighbour entry"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links("srcguid", {"guid1": {"dist": 12}})
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 1)
        res = p.db.guid2neighbour.count_documents({"guid": "srcguid"})
        self.assertEqual(res, 1)
        a, b, c = p._audit_storage("srcguid")
        self.assertEqual(len(b.index), 1)
        self.assertEqual(set(b["rstat"]), set(["s"]))  # singleton


@mongo_test
class Test_SeqMeta_audit_storage_2(unittest.TestCase):
    """tests audit of the situation where there are no links"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links("srcguid", {})
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 0)
        a, b, c = p._audit_storage("srcguid")
        self.assertEqual(a, set([]))
        self.assertIsInstance(b, pd.DataFrame)
        self.assertEqual(
            c,
            {
                "singletons": 0,
                "stored_in_m_records": 0,
                "stored_in_f_records": 0,
                "f_records_with_duplicates": 0,
                "neighbouring_guids": 0,
            },
        )


@mongo_test
class Test_SeqMeta_guid2neighbour_2(unittest.TestCase):
    """tests creation of a new guid2neighbour entry"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links("srcguid", {})
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 0)


@mongo_test
class Test_SeqMeta_version(unittest.TestCase):
    """tests version of library.  only tested with > v3.0"""

    def runTest(self):
        self.assertTrue(pymongo.__version__ >= "3.0")


@mongo_test
class Test_SeqMeta_file_store1(unittest.TestCase):
    """tests storage of sequence objects database"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        obj1 = {1, 2, 3}
        guid = "guid1"
        p.rcs.delete({"filename": guid})  # delete if present
        p.refcompressedseq_store(guid, obj1)
        res = p.rcs.find_one({"filename": guid}).read()
        obj2 = pickle.loads(res)
        self.assertEqual(obj1, obj2)


@mongo_test
class Test_SeqMeta_file_store2(unittest.TestCase):
    """tests storage of pickle files in database"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        obj1 = {1, 2, 3}
        guid = "guid1"
        p.rcs.delete({"filename": guid})  # delete if present
        pickled_obj = pickle.dumps(obj1, protocol=2)
        self.assertTrue(pickled_obj is not None)

        p.refcompressedseq_store(guid, obj1)
        with self.assertRaises(FileExistsError):
            p.refcompressedseq_store(guid, obj1)


@mongo_test
class Test_SeqMeta_file_store3(unittest.TestCase):
    """tests storage of pickle files in database"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        obj1 = {1, 2, 3}
        guid = "guid1"
        p.rcs.delete({"filename": "guid1"})  # delete if present
        p.rcs.delete({"filename": "guid2"})  # delete if present
        res1 = p.refcompressedsequence_guids()
        pickled_obj = pickle.dumps(obj1, protocol=2)
        p.refcompressedseq_store(guid, pickled_obj)
        guid = "guid2"
        p.refcompressedseq_store(guid, pickled_obj)
        res2 = p.refcompressedsequence_guids()
        self.assertEqual(res2 - res1, set(["guid1", "guid2"]))


@mongo_test
class Test_SeqMeta_guid_annotate_1(unittest.TestCase):
    """tests insert of new data item"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        #  add the dictionary 'payload' to the namespace 'ns'
        guid = 1
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload)


@mongo_test
class Test_SeqMeta_guid_annotate_2(unittest.TestCase):
    """tests addition to existing data item"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = 1

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload)

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"three": 3}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], {"one": 1, "two": 2, "three": 3})


@mongo_test
class Test_SeqMeta_guid_annotate_3(unittest.TestCase):
    """tests addition to existing data item"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = 1

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload)

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"two": 3}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], {"one": 1, "two": 3})


@mongo_test
class Test_SeqMeta_guid_valid_1(unittest.TestCase):
    """tests insert of new data item and validity check"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = "valid"
        namespace = "DNAQuality"
        payload = {"invalid": 0}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "invalid"
        namespace = "DNAQuality"
        payload = {"invalid": 1}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "missing"
        namespace = "DNAQuality"
        payload = {"N": 1}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

        res = p.guid_valid("valid")
        self.assertEqual(res, 0)
        res = p.guid_valid("invalid")
        self.assertEqual(res, 1)
        res = p.guid_valid("missing")
        self.assertEqual(res, -2)
        res = p.guid_valid("noguid")
        self.assertEqual(res, -1)

        valid_guids = p.guids_valid()
        self.assertEqual(valid_guids, set(["valid"]))
        invalid_guids = p.guids_invalid()
        self.assertEqual(invalid_guids, set(["invalid"]))


@mongo_test
class Test_SeqMeta_guid_valid_2(unittest.TestCase):
    """tests insert of new data item and validity check"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = "valid1"
        namespace = "DNAQuality"
        payload = {"invalid": 0}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "valid2"
        namespace = "DNAQuality"
        payload = {"invalid": 0}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

        guid = "invalid"
        namespace = "DNAQuality"
        payload = {"invalid": 1}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "missing"
        namespace = "DNAQuality"
        payload = {"N": 1}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

        res = p.guids_valid()
        self.assertEqual(res, set(["valid1", "valid2"]))
        res = p.guids_invalid()
        self.assertEqual(res, set(["invalid"]))


@mongo_test
class Test_SeqMeta_guid_annotate_5(unittest.TestCase):
    """tests update of existing data item with same namespace"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace = "ns"
        payload1 = {"one": 1, "two": 2}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload1)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload1)
        payload2 = {"one": 1, "two": 2}
        p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload2)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload2)


@mongo_test
class Test_SeqMeta_guid_annotate_6(unittest.TestCase):
    """tests update of existing data item with different namespace"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        payload1 = {"one": 1, "two": 2}

        p.guid_annotate(guid=guid, nameSpace="ns1", annotDict=payload1)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns1"], payload1)

        payload2 = {"one": 1, "two": 2}
        p.guid_annotate(guid=guid, nameSpace="ns2", annotDict=payload2)
        res = p.db.guid2meta.find_one({"_id": 1})

        payloads = {"ns1": payload1, "ns2": payload2}
        self.assertEqual(res["sequence_meta"], payloads)


@mongo_test
class Test_SeqMeta_init(unittest.TestCase):
    """tests database creation"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        self.assertTrue(p.first_run() is True)
        res = p.config_read("preComparer")
        self.assertEqual(res, None)

        p.config_store("config", {"item": 1})
        self.assertTrue(p.first_run() is False)
        res = p.config_read("config")
        self.assertEqual(res, {"_id": "config", "item": 1})

        p.config_store("preComparer", {"item": 2})
        res = p.config_read("config")
        self.assertEqual(res, {"_id": "config", "item": 1})
        res = p.config_read("preComparer")
        self.assertEqual(res, {"_id": "preComparer", "item": 2})

        p.config_store("preComparer", {"item": 3})
        res = p.config_read("config")
        self.assertEqual(res, {"_id": "config", "item": 1})
        res = p.config_read("preComparer")
        self.assertEqual(res, {"_id": "preComparer", "item": 3})


@mongo_test
class Test_SeqMeta_guids(unittest.TestCase):
    """tests recovery of sequence guids"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        startup = {"_id": 1}
        res = p.db.guid2meta.insert_one(startup)
        startup = {"_id": 2}
        res = p.db.guid2meta.insert_one(startup)
        res = p.guids()
        self.assertEqual(res, set([1, 2]))


@mongo_test
class Test_SeqMeta_Base(unittest.TestCase):
    """sets up a connection for unit testing"""

    def setUp(self):
        self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        self.assertTrue(self.t.first_run() is True)


@mongo_test
class Test_SeqMeta_guid_quality_check_1(Test_SeqMeta_Base):
    def runTest(self):
        """tests return of sequences and their qualities"""
        # set up nucleic acid object
        na = NucleicAcid()
        na.examine("ACGTACGTNN")  # 20% bad

        self.t.guid_annotate(
            guid="g1", nameSpace="DNAQuality", annotDict=na.composition
        )
        na.examine("ACGTACNNNN")  # 40% bad
        self.t.guid_annotate(
            guid="g2", nameSpace="DNAQuality", annotDict=na.composition
        )
        na.examine("ACGTNNNNNN")  # 60% bad
        self.t.guid_annotate(
            guid="g3", nameSpace="DNAQuality", annotDict=na.composition
        )

        r1 = self.t.guid_quality_check("g1", 0.80)  # valid
        r2 = self.t.guid_quality_check("g2", 0.80)  # invalid
        r3 = self.t.guid_quality_check("g3", 0.80)  # invalid
        r4 = self.t.guid_quality_check("g4", 0.80)  # invalid; does not exist.

        self.assertEqual(r1, True)
        self.assertEqual(r2, False)
        self.assertEqual(r3, False)
        self.assertEqual(r4, None)


@mongo_test
class Test_SeqMeta_guid2quality1(Test_SeqMeta_Base):
    def runTest(self):
        """tests return of sequences and their qualities"""
        # set up nucleic acid object
        na = NucleicAcid()
        na.examine("ACGTACGTNN")  # 20% bad

        self.t.guid_annotate(
            guid="g1", nameSpace="DNAQuality", annotDict=na.composition
        )
        na.examine("ACGTACNNNN")  # 40% bad
        self.t.guid_annotate(
            guid="g2", nameSpace="DNAQuality", annotDict=na.composition
        )
        na.examine("ACGTNNNNNN")  # 60% bad
        self.t.guid_annotate(
            guid="g3", nameSpace="DNAQuality", annotDict=na.composition
        )

        r1 = self.t.guid_quality_check("g1", 0.80)  # valid
        r2 = self.t.guid_quality_check("g2", 0.80)  # invalid
        r3 = self.t.guid_quality_check("g3", 0.80)  # invalid
        r4 = self.t.guid_quality_check("g4", 0.80)  # invalid; does not exist.

        self.assertEqual(r1, True)
        self.assertEqual(r2, False)
        self.assertEqual(r3, False)
        self.assertEqual(r4, None)

        resDict = self.t.guid2quality(None)  # restrict to nothing - return all
        self.assertTrue(resDict is not None)
        assert resDict is not None  # for typing purposes
        self.assertEqual(resDict["g1"], 0.80)
        self.assertEqual(resDict["g2"], 0.60)
        self.assertEqual(resDict["g3"], 0.40)


@mongo_test
class Test_SeqMeta_guid2quality2(Test_SeqMeta_Base):
    def runTest(self):
        """tests return of sequences and their qualities"""
        # set up nucleic acid object

        na = NucleicAcid()
        na.examine("ACGTACGTNN")  # 20% bad

        self.t.guid_annotate(
            guid="g1", nameSpace="DNAQuality", annotDict=na.composition
        )
        na.examine("ACGTACNNNN")  # 40% bad
        self.t.guid_annotate(
            guid="g2", nameSpace="DNAQuality", annotDict=na.composition
        )
        na.examine("ACGTNNNNNN")  # 60% bad
        self.t.guid_annotate(
            guid="g3", nameSpace="DNAQuality", annotDict=na.composition
        )

        r1 = self.t.guid_quality_check("g1", 0.80)  # valid
        r2 = self.t.guid_quality_check("g2", 0.80)  # invalid
        r3 = self.t.guid_quality_check("g3", 0.80)  # invalid

        self.assertEqual(r1, True)
        self.assertEqual(r2, False)
        self.assertEqual(r3, False)  # check the db insert works

        resDict = self.t.guid2quality(["g1", "g2", "g3"])
        self.assertTrue(resDict is not None)
        assert resDict is not None  # for typing purposes
        self.assertEqual(resDict["g1"], 0.80)
        self.assertEqual(resDict["g2"], 0.60)
        self.assertEqual(resDict["g3"], 0.40)


@mongo_test
class Test_SeqMeta_Base1(unittest.TestCase):
    """initialise FN persistence and adds data"""

    def setUp(self):
        self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        dna = NucleicAcid()

        # add some sequences
        seqs = {"guid1": "ACGT", "guid2": "NACT", "guid3": "TTTT", "guid4": "NNNN"}
        for guid in seqs.keys():
            seq = seqs[guid]
            dna.examine(seq)
            self.t.guid_annotate(
                guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
            )


@mongo_test
class Test_SeqMeta_Base1t(unittest.TestCase):
    """initialise FN persistence and adds data, 0.1 secs apart.
    Used for testing queries examining order of recovery of samples."""

    def setUp(self):
        self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        dna = NucleicAcid()

        # add some sequences
        seqs = {"guid1": "ACGT", "guid2": "NACT", "guid3": "TTTT", "guid4": "NNNN"}
        for guid in seqs.keys():
            time.sleep(0.1)
            seq = seqs[guid]
            dna.examine(seq)
            self.t.guid_annotate(
                guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
            )
        self.seqs = seqs


@mongo_test
class Test_SeqMeta_guid2ExaminationDateTime(Test_SeqMeta_Base1):
    """recovering guids and examination times;"""

    def runTest(self):
        res = self.t.guid2ExaminationDateTime()
        self.assertIsNotNone(res)
        assert res is not None  # for typing purposes
        expected = 4
        self.assertEqual(len(res.keys()), expected)


@mongo_test
class Test_SeqMeta_guid2ExaminationDateTime_order(Test_SeqMeta_Base1t):
    """tests guid2ExaminationDateTime"""

    def runTest(self):
        res = self.t.guid2ExaminationDateTime()
        self.assertIsNotNone(res)
        assert res is not None  # for typing purposes
        expected = 4
        self.assertEqual(len(res.keys()), expected)

        # check that the sample were added in order, with increasing examination times.
        previous_addition_time = None
        for i, guid in enumerate(sorted(self.seqs.keys())):  # the order added
            if i > 0:
                self.assertGreater(res[guid], previous_addition_time)
            previous_addition_time = res[guid]


@mongo_test
class Test_SeqMeta_guid_examination_time(Test_SeqMeta_Base1t):
    """tests guid_examination_time()"""

    def runTest(self):
        res = self.t.guid2ExaminationDateTime()
        self.assertIsNotNone(res)
        assert res is not None  # for typing purposes
        expected = 4
        self.assertEqual(len(res.keys()), expected)

        # check that the sample were added in order, with increasing examination times.
        previous_addition_time = None
        for i, guid in enumerate(sorted(self.seqs.keys())):  # the order added
            this_examination_time = self.t.guid_examination_time(guid)
            if i > 0:
                self.assertGreater(this_examination_time, previous_addition_time)
            previous_addition_time = this_examination_time

        this_examination_time = self.t.guid_examination_time("missing-guid")
        self.assertIsNone(this_examination_time)


@mongo_test
class Test_SeqMeta_guid_considered_after(Test_SeqMeta_Base1t):
    """recovering guids and examination times;"""

    def runTest(self):
        res = self.t.guid2ExaminationDateTime()
        self.assertIsNotNone(res)
        assert res is not None  # for typing purposes
        expected = 4
        self.assertEqual(len(res.keys()), expected)

        # check that the sample were added in order, with increasing examination times.
        for i, guid in enumerate(sorted(self.seqs.keys())):  # the order added
            this_examination_time = self.t.guid_examination_time(guid)
            self.assertIsNotNone(this_examination_time)
            assert this_examination_time is not None  # for typing purposes
            res2 = self.t.guids_considered_after(this_examination_time)
            self.assertEqual(
                len(res2), 3 - i
            )  # with guid1, we expect three; with guid2, we expect 2; etc

            res2 = self.t.guids_considered_after_guid(guid)
            self.assertEqual(
                len(res2), 3 - i
            )  # with guid1, we expect three; with guid2, we expect 2; etc


@mongo_test
class Test_SeqMeta_guid_added_after(Test_SeqMeta_Base1t):
    """recovering guids added after another guid;"""

    def runTest(self):
        res = self.t.guids_added_after_sample("noguid")  # should be None
        self.assertIsNone(res)

        res = self.t.guids_added_after_sample("guid4")  # should be empty set
        self.assertEqual(set([]), res)

        res = self.t.guids_added_after_sample("guid1")  # should be empty set
        self.assertEqual(set(["guid2", "guid3", "guid4"]), res)

        res = self.t.guids_added_after_sample("guid2")  # should be empty set
        self.assertEqual(set(["guid3", "guid4"]), res)

        res = self.t.guids_added_after_sample("guid3")  # should be empty set
        self.assertEqual(set(["guid4"]), res)


@mongo_test
class Test_SeqMeta_propACTG_filteredSequenceGuids(Test_SeqMeta_Base1):
    """recovered guids filtered by the propACTG criterion"""

    def runTest(self):
        n = 0
        for guid in self.t.guid2propACTG_filtered(cutoff=0.85):
            n += 1
        expected = 2
        self.assertEqual(n, expected)


@mongo_test
class Test_SeqMeta_allAnnotations(Test_SeqMeta_Base1):
    """tests recovery of all annoations"""

    def runTest(self):
        df = self.t.guid_annotations()
        self.assertIsNotNone(df)
        assert df is not None  # for typing purposes
        self.assertEqual(len(df.keys()), 4)


@mongo_test
class Test_SeqMeta_oneAnnotation(Test_SeqMeta_Base1):
    """tests recovery of one annotations"""

    def runTest(self):
        df = self.t.guid_annotation("guid3")
        self.assertIsNotNone(df)
        assert df is not None  # for typing purposes
        self.assertEqual(len(df.keys()), 1)
        df = self.t.guid_annotation("missing")
        self.assertIsNotNone(df)
        assert df is not None  # for typing purposes
        self.assertEqual(len(df.keys()), 0)


@mongo_test
class Test_Clusters(unittest.TestCase):
    """tests saving and recovery of dictionaries to Clusters"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        self.assertIsNone(p.cluster_latest_version("cl1"))

        self.assertEqual(0, len(p.cluster_versions("cl1")))

        payload1 = {"one": 1, "two": 2}
        x = p.cluster_store("cl1", payload1)

        payload1b = {"2_one": 1, "2_two": 2}
        y = p.cluster_store("cl2", payload1b)
        self.assertTrue(y is not None)

        self.assertIsNotNone(p.cluster_latest_version("cl1"))
        clv = p.cluster_latest_version("cl1")
        self.assertEqual(x, clv)
        self.assertEqual(1, len(p.cluster_versions("cl1")))
        self.assertIsNone(p.cluster_read_update("cl1", clv))

        payload2 = p.cluster_read("cl1")
        self.assertEqual(payload1, payload2)

        payload3 = {"one": 10, "two": 20}
        p.cluster_store("cl1", payload3)  # this is now the latest version

        payload4 = {"2-one": 10, "2-two": 20}
        p.cluster_store("cl2", payload4)  # this is now the latest version

        self.assertEqual(p.cluster_keys(), ["cl1", "cl2"])
        self.assertEqual(p.cluster_keys(clustering_name="cl1"), ["cl1"])

        self.assertEqual(2, len(p.cluster_versions("cl1")))
        self.assertNotEqual(clv, p.cluster_latest_version("cl1"))
        self.assertIsNotNone(p.cluster_read_update("cl1", clv))

        payload5 = p.cluster_read("cl1")
        self.assertEqual(payload5, payload3)

        p.cluster_delete_legacy_by_key("cl1")

        self.assertEqual(1, len(p.cluster_versions("cl1")))

        p.cluster_delete_all("cl1")

        self.assertEqual(0, len(p.cluster_versions("cl1")))
        self.assertEqual(2, len(p.cluster_versions("cl2")))
        p.cluster_delete_legacy_by_key("cl2")
        self.assertEqual(1, len(p.cluster_versions("cl2")))


@mongo_test
class Test_Tree(unittest.TestCase):
    """tests saving and recovery of dictionaries to Tree"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        payload1 = {"one": 1, "two": 2}
        p.tree_store(tree_token="tree1", tree=payload1)
        payload2 = p.tree_read(tree_token="tree1")
        self.assertEqual(payload1, payload2)
        p.tree_delete(tree_token="tree1")
        payload3 = p.tree_read(tree_token="tree1")
        self.assertIsNone(payload3)

        payload1 = {"one": 1, "two": 2}
        p.tree_store(tree_token="tree1", tree=payload1)
        payload2 = {"one": 3, "two": 4}
        p.tree_store(tree_token="tree2", tree=payload2)
        self.assertEqual(2, len(p.tree_stored_ids()))
        p.tree_delete_unless_whitelisted(whitelist=["tree1", "tree2"])
        self.assertEqual(2, len(p.tree_stored_ids()))
        p.tree_delete_unless_whitelisted(whitelist=["tree1"])
        self.assertEqual(1, len(p.tree_stored_ids()))


@mongo_test
class Test_MSA(unittest.TestCase):
    """tests saving and recovery of dictionaries to MSA"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        payload1 = {"one": 1, "two": 2}
        p.msa_store(msa_token="msa1", msa=payload1)
        payload2 = p.msa_read(msa_token="msa1")
        self.assertEqual(payload1, payload2)
        p.msa_delete(msa_token="msa1")
        payload3 = p.msa_read(msa_token="msa1")
        self.assertIsNone(payload3)

        payload1 = {"one": 1, "two": 2}
        p.msa_store(msa_token="msa1", msa=payload1)
        payload2 = {"one": 3, "two": 4}
        p.msa_store(msa_token="msa2", msa=payload2)
        self.assertEqual(2, len(p.msa_stored_ids()))
        p.msa_delete_unless_whitelisted(whitelist=["msa1", "msa2"])
        self.assertEqual(2, len(p.msa_stored_ids()))
        p.msa_delete_unless_whitelisted(whitelist=["msa1"])
        self.assertEqual(1, len(p.msa_stored_ids()))


@mongo_test
class Test_Monitor(unittest.TestCase):
    """tests saving and recovery of strings to monitor"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        payload1 = "line1"
        p.monitor_store("r1", payload1)
        payload2 = p.monitor_read("r1")
        self.assertEqual(payload1, payload2)
        payload3 = p.monitor_read("nil")
        self.assertIsNone(payload3)


@mongo_test
class test_Raise_error(unittest.TestCase):
    """tests raise_error"""

    def runTest(self):
        # generate compressed sequences
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        with self.assertRaises(ZeroDivisionError):
            p.raise_error("token")


@mongo_test
class Test_summarise_stored_items(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.summarise_stored_items()
        self.assertTrue(res is not None)


@mongo_test
class Test_rotate_log(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.rotate_log()


@mongo_test
class Test_lock_1(unittest.TestCase):
    """tests locking.

    Note: does not test concurrent operations"""

    def runTest(self):

        pdm = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        self.assertTrue(pdm.unlock(1, force=True))
        self.assertEqual(0, pdm.lock_status(1)["lock_status"])
        self.assertEqual("-NotSpecified-", pdm.lock_status(1)["sequence_id"])

        res = pdm.lock_details(1)
        self.assertIsNone(res)

        self.assertTrue(pdm.unlock(0, force=True))
        self.assertEqual(0, pdm.lock_status(0)["lock_status"])
        self.assertEqual("-NotSpecified-", pdm.lock_status(0)["sequence_id"])
        res = pdm.lock_details(0)
        self.assertIsNone(res)

        self.assertTrue(pdm.lock(1, "guid1"))  # lock open; should succeed
        self.assertEqual(1, pdm.lock_status(1)["lock_status"])
        self.assertEqual("guid1", pdm.lock_status(1)["sequence_id"])
        res = pdm.lock_details(1)
        self.assertEqual(res["sequence_id"], "guid1")
        self.assertIsInstance(res["uuid"], str)
        self.assertIsInstance(res["lock_set_date"], datetime.datetime)

        self.assertTrue(pdm.lock(0, "guid0"))  # lock open; should succeed

        self.assertFalse(pdm.lock(1, "guid3"))  # lock closed; should fail
        self.assertEqual(1, pdm.lock_status(1)["lock_status"])
        self.assertEqual("guid1", pdm.lock_status(1)["sequence_id"])
        res = pdm.lock_details(1)
        self.assertEqual(res["sequence_id"], "guid1")
        self.assertIsInstance(res["uuid"], str)
        self.assertIsInstance(res["lock_set_date"], datetime.datetime)

        self.assertTrue(pdm.unlock(1))  # lock closed should succeed
        self.assertEqual(0, pdm.lock_status(1)["lock_status"])
        self.assertEqual("-NotSpecified-", pdm.lock_status(1)["sequence_id"])
        res = pdm.lock_details(1)
        self.assertIsNone(res)

        self.assertTrue(pdm.unlock(2, force=True))
        self.assertEqual(0, pdm.lock_status(2)["lock_status"])
        self.assertEqual("-NotSpecified-", pdm.lock_status(2)["sequence_id"])
        self.assertTrue(pdm.lock(2, "start_catwalk"))
        self.assertTrue(pdm.unlock(2))


class test_lockmanager_2(unittest.TestCase):
    """tests whether the lockmanager runs and releases a lock running lock"""

    def runTest(self):

        pdm = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        # unlock
        pdm.unlock(1, force=True)

        # there is no lock
        lock_details = pdm.lock_details(1)
        self.assertIsNone(lock_details)

        # take a lock
        pdm.lock(1, "guid1")

        # check the lock is held
        lock_details = pdm.lock_details(1)
        self.assertIsNotNone(lock_details)

        os.system(
            "pipenv run python3 findNeighbour4_lockmanager.py --max_run_time 10 --run_once_only"
        )

        # check the lock is still held
        lock_details = pdm.lock_details(1)
        self.assertIsNotNone(lock_details)

        # unlock
        pdm.unlock(1)

        # there is no lock
        lock_details = pdm.lock_details(1)
        self.assertIsNone(lock_details)

        # take a lock
        pdm.lock(1, "guid1")
        lock_details = pdm.lock_details(1)
        self.assertIsNotNone(lock_details)

        # wait 15 seconds
        time.sleep(15)
        os.system(
            "pipenv run python3 findNeighbour4_lockmanager.py --max_run_time 10 --run_once_only"
        )

        # there is no lock
        lock_details = pdm.lock_details(1)
        self.assertIsNone(lock_details)


class Test_refcompressedseq_read_many(unittest.TestCase):
    """tests read_many method for referencecompressed sequences"""

    def runTest(self):

        # define a sequence object for testing
        seqobj1 = {
            "A": set([1]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 0,
        }
        seqobj2 = {
            "A": set([1, 2]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 0,
        }
        # define a sequence object for testing
        seqobj3 = {
            "A": set([1, 2, 3]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 0,
        }

        payloads = {"guid1": seqobj1, "guid2": seqobj2, "guid3": seqobj3}

        pdm = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        for guid in payloads.keys():
            pdm.refcompressedseq_store(guid, payloads[guid])

        i = 0
        for guid, rcs in pdm.refcompressedsequence_read_many(["guid1", "guid3"]):
            self.assertEqual(rcs, payloads[guid])
            self.assertTrue(guid in ["guid1", "guid3"])
            i += 1

        self.assertEqual(i, 2)


class Test_refcompressedseq_read_all(unittest.TestCase):
    """tests read_all method for referencecompressed sequences"""

    def runTest(self):

        # define a sequence object for testing
        seqobj1 = {
            "A": set([1]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 0,
        }
        seqobj2 = {
            "A": set([1, 2]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 0,
        }
        # define a sequence object for testing
        seqobj3 = {
            "A": set([1, 2, 3]),
            "C": set([6]),
            "T": set([4]),
            "G": set([5]),
            "M": {11: "Y", 12: "k"},
            "invalid": 0,
        }

        payloads = {"guid1": seqobj1, "guid2": seqobj2, "guid3": seqobj3}
        pdm = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        for guid in payloads.keys():
            print("STORING", guid)
            pdm.refcompressedseq_store(guid, payloads[guid])

        i = 0
        for guid, rcs in pdm.refcompressedsequence_read_all(internal_batch_size=1):
            self.assertEqual(rcs, payloads[guid])
            self.assertTrue(guid in ["guid1", "guid2", "guid3"])
            i += 1
        self.assertEqual(i, 3)
        i = 0
        for guid, rcs in pdm.refcompressedsequence_read_all(internal_batch_size=2):
            self.assertEqual(rcs, payloads[guid])
            self.assertTrue(guid in ["guid1", "guid2", "guid3"])
            i += 1
        self.assertEqual(i, 3)
        i = 0
        for guid, rcs in pdm.refcompressedsequence_read_all(internal_batch_size=3):
            self.assertEqual(rcs, payloads[guid])
            self.assertTrue(guid in ["guid1", "guid2", "guid3"])
            i += 1
        self.assertEqual(i, 3)
        i = 0
        for guid, rcs in pdm.refcompressedsequence_read_all(internal_batch_size=4):
            self.assertEqual(rcs, payloads[guid])
            self.assertTrue(guid in ["guid1", "guid2", "guid3"])
            i += 1
        self.assertEqual(i, 3)
