""" test mongo db connection object """
import unittest
import time
import json
import pymongo
import pandas as pd
import pickle
from findn.NucleicAcid import NucleicAcid


from findn.mongoStore import fn3persistence

## persistence unit tests
UNITTEST_MONGOCONN = "mongodb://localhost"


class Test_Server_Monitoring_0(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.server_monitoring_store(message="one")

        res = p.recent_server_monitoring(100)

        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))


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
            res = p.recent_server_monitoring("thing")


class Test_Server_Monitoring_3(unittest.TestCase):
    """checks whether server_monitoring_min_interval_msec control works"""

    def runTest(self):
        p = fn3persistence(
            connString=UNITTEST_MONGOCONN,
            debug=2,
            server_monitoring_min_interval_msec=2000,
        )
        retVal = p.server_monitoring_store(message="one")  # should insert
        self.assertEqual(retVal, True)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))

        retVal = p.server_monitoring_store(message="two")  # should not inserted
        self.assertEqual(retVal, False)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 1)
        self.assertTrue(isinstance(res, list))

        time.sleep(2)  # seconds
        retVal = p.server_monitoring_store(message="three")  # should insert
        self.assertEqual(retVal, True)
        res = p.recent_server_monitoring(100)
        self.assertEqual(len(res), 2)
        self.assertTrue(isinstance(res, list))


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


class Test_SeqMeta_guid2neighbour_2(unittest.TestCase):
    """tests creation of a new guid2neighbour entry"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.guid2neighbour_add_links("srcguid", {})
        res = p.db.guid2neighbour.count_documents({"guid": "guid1"})
        self.assertEqual(res, 0)


class Test_SeqMeta_version(unittest.TestCase):
    """tests version of library.  only tested with > v3.0"""

    def runTest(self):
        self.assertTrue(pymongo.__version__ >= "3.0")


class Test_SeqMeta_file_store1(unittest.TestCase):
    """tests storage of pickle files in database"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        obj1 = {1, 2, 3}
        guid = "guid1"
        p.rcs.delete({"filename": guid})  # delete if present
        p.refcompressedseq_store(guid, obj1)
        res = p.rcs.find_one({"filename": guid}).read()
        obj2 = pickle.loads(res)
        self.assertEqual(obj1, obj2)


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


class Test_SeqMeta_guid_annotate_1(unittest.TestCase):
    """tests insert of new data item"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        #  add the dictionary 'payload' to the namespace 'ns'
        guid = 1
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload)


class Test_SeqMeta_guid_annotate_2(unittest.TestCase):
    """tests addition to existing data item"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = 1

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload)

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"three": 3}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], {"one": 1, "two": 2, "three": 3})


class Test_SeqMeta_guid_annotate_3(unittest.TestCase):
    """tests addition to existing data item"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = 1

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload)

        # add the dictionary 'payload' to the namespace 'ns'
        namespace = "ns"
        payload = {"two": 3}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], {"one": 1, "two": 3})


class Test_SeqMeta_guid_exists_1(unittest.TestCase):
    """tests insert of new data item and existence check"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace = "ns"
        payload = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        res = p.guid_exists(guid)
        self.assertEqual(res, True)
        res = p.guid_exists(-1)
        self.assertEqual(res, False)


class Test_SeqMeta_guid_valid_1(unittest.TestCase):
    """tests insert of new data item and validity check"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = "valid"
        namespace = "DNAQuality"
        payload = {"invalid": 0}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "invalid"
        namespace = "DNAQuality"
        payload = {"invalid": 1}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "missing"
        namespace = "DNAQuality"
        payload = {"N": 1}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

        res = p.guid_valid("valid")
        self.assertEqual(res, 0)
        res = p.guid_valid("invalid")
        self.assertEqual(res, 1)
        res = p.guid_valid("missing")
        self.assertEqual(res, -2)
        res = p.guid_valid("noguid")
        self.assertEqual(res, -1)


class Test_SeqMeta_guid_valid_2(unittest.TestCase):
    """tests insert of new data item and validity check"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        guid = "valid1"
        namespace = "DNAQuality"
        payload = {"invalid": 0}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "valid2"
        namespace = "DNAQuality"
        payload = {"invalid": 0}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

        guid = "invalid"
        namespace = "DNAQuality"
        payload = {"invalid": 1}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
        guid = "missing"
        namespace = "DNAQuality"
        payload = {"N": 1}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

        res = p.guids_valid()
        self.assertEqual(res, set(["valid1", "valid2"]))
        res = p.guids_invalid()
        self.assertEqual(res, set(["invalid"]))


class Test_SeqMeta_guid_annotate_5(unittest.TestCase):
    """tests update of existing data item with same namespace"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        namespace = "ns"
        payload1 = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload1)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload1)
        payload2 = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload2)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns"], payload2)


class Test_SeqMeta_guid_annotate_6(unittest.TestCase):
    """tests update of existing data item with different namespace"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)

        # test there is no 'test' item; insert, and confirm insert
        guid = 1
        payload1 = {"one": 1, "two": 2}

        res = p.guid_annotate(guid=guid, nameSpace="ns1", annotDict=payload1)
        res = p.db.guid2meta.find_one({"_id": 1})
        self.assertEqual(res["sequence_meta"]["ns1"], payload1)

        payload2 = {"one": 1, "two": 2}
        res = p.guid_annotate(guid=guid, nameSpace="ns2", annotDict=payload2)
        res = p.db.guid2meta.find_one({"_id": 1})

        payloads = {"ns1": payload1, "ns2": payload2}
        self.assertEqual(res["sequence_meta"], payloads)


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


class Test_SeqMeta_Base(unittest.TestCase):
    """sets up a connection for unit testing"""

    def setUp(self):
        self.t = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        self.assertTrue(self.t.first_run() is True)


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
        self.assertEqual(resDict["g1"], 0.80)
        self.assertEqual(resDict["g2"], 0.60)
        self.assertEqual(resDict["g3"], 0.40)


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
        self.assertEqual(resDict["g1"], 0.80)
        self.assertEqual(resDict["g2"], 0.60)
        self.assertEqual(resDict["g3"], 0.40)


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


class Test_SeqMeta_guid2ExaminationDateTime(Test_SeqMeta_Base1):
    """recovering guids and examination times;"""

    def runTest(self):
        res = self.t.guid2ExaminationDateTime()
        expected = 4
        self.assertEqual(len(res.keys()), expected)


class Test_SeqMeta_propACTG_filteredSequenceGuids(Test_SeqMeta_Base1):
    """recovered guids filtered by the propACTG criterion"""

    def runTest(self):
        n = 0
        for guid in self.t.guid2propACTG_filtered(cutoff=0.85):
            n += 1
        expected = 2
        self.assertEqual(n, expected)


class Test_SeqMeta_allAnnotations(Test_SeqMeta_Base1):
    """tests recovery of all annoations"""

    def runTest(self):
        df = self.t.guid_annotations()
        self.assertEqual(len(df.keys()), 4)


class Test_SeqMeta_oneAnnotation(Test_SeqMeta_Base1):
    """tests recovery of one annotations"""

    def runTest(self):
        df = self.t.guid_annotation("guid3")
        self.assertEqual(len(df.keys()), 1)
        df = self.t.guid_annotation("missing")
        self.assertEqual(len(df.keys()), 0)


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

        payload4 = p.cluster_read("cl1")
        self.assertEqual(payload4, payload3)

        p.cluster_delete_legacy_by_key("cl1")

        self.assertEqual(1, len(p.cluster_versions("cl1")))

        p.cluster_delete_all("cl1")

        self.assertEqual(0, len(p.cluster_versions("cl1")))
        self.assertEqual(2, len(p.cluster_versions("cl2")))
        p.cluster_delete_legacy_by_key("cl2")
        self.assertEqual(1, len(p.cluster_versions("cl2")))


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


class test_Raise_error(unittest.TestCase):
    """tests raise_error"""

    def runTest(self):
        # generate compressed sequences
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        with self.assertRaises(ZeroDivisionError):
            p.raise_error("token")


class Test_summarise_stored_items(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        res = p.summarise_stored_items()
        self.assertTrue(res is not None)


class Test_rotate_log(unittest.TestCase):
    """adds server monitoring info"""

    def runTest(self):
        p = fn3persistence(connString=UNITTEST_MONGOCONN, debug=2)
        p.rotate_log()
