""" tests rdbms backend to findneighbour4

    to invoke, either run automatically with pytest or do
    pipenv run python3 -m unittest test/test_rdbmsstore.py

"""
import os
import json
import time
import unittest
import pandas as pd
from sqlalchemy import func
from findn.NucleicAcid import NucleicAcid
from findn.rdbmsstore import fn3persistence, BulkLoadTest, RDBMSError, RefCompressedSeq


class Test_Database(unittest.TestCase):
    """establishes database connection strings for cross-database testing.
    Currently tests OCI (if relevant environment variables are set) and Sqlite

    To test other databases, such as MySql, add the relevant connection string &
    database name to the dictionary self.engines"""

    def setUp(self):
        self.engines = {}
        self.engines["Sqlite"] = "sqlite://"  # in memory sqlite

        conn_detail_file = None

        # try to read the environment variable 'DB_CONNECTION_CONFIG_FILE'
        try:
            conn_detail_file = os.environ["DB_CONNECTION_CONFIG_FILE"]
        except KeyError:
            # doesn't exist; we just run with sqlite, which is the default if engine is None.
            print(
                "No environment variable DB_CONNECTION_CONFIG_FILE found.  Testing with sqlite only."
            )
        if conn_detail_file is None:
            print(
                "No environment variable DB_CONNECTION_CONFIG_FILE found.  Testing with sqlite only."
            )
        else:
            if not os.path.exists(conn_detail_file):
                raise FileNotFoundError(
                    "Connection file specified but not found: {0}".format(
                        conn_detail_file
                    )
                )
            # try to read the config file
            with open(conn_detail_file, "rt") as f:
                conn_detail = json.load(f)
                for key in conn_detail.keys():
                    if key.startswith("unittest_ora"):
                        self.engines[key] = key
                        pass

    def pdms(self, **kwargs):
        """yields fn3persistence objects, one for each database server being tested."""
        for engine, config in self.engines.items():
            print(engine, type(self).__name__)
            pdm = fn3persistence(connection_config=config, debug=2, **kwargs)
            yield pdm

            # explicitly close connection (required for unittesting)
            pdm.explicitly_close_connections()


class Test_to_string(Test_Database):
    """tests conversion of bytes to string transparently"""

    def runTest(self):
        for pdm in self.pdms():
            x = "david"
            r1 = pdm._to_string(x)
            y = b"david"
            r2 = pdm._to_string(y)
            self.assertIsInstance(r1, str)
            self.assertIsInstance(r2, str)


#@unittest.skip("Fails on Oracle with ORA-00054 error, related to schema modification")
class Test_create_database_1(Test_Database):
    """tests creating the database and internal functions dropping tables"""

    def runTest(self):
        expected_tables = set(
            [
                "fn4_bulk_load_test",
                "refcompressedseq",
                "edge",
                "sequence_cluster",
                "monitor",
                "server_monitoring",
                "msa",
                "tree",
            ]
        )
        for pdm in self.pdms():
            n1 = pdm._table_names()
            pdm._drop_existing_tables()
            n2 = pdm._table_names()
            self.assertEqual(len(expected_tables.intersection(set(n1))), 8)
            self.assertEqual(len(expected_tables.intersection(set(n2))), 0)


class Test_oracle_bulk_upload_2(Test_Database):
    """tests bulk upload with small number of samples"""

    def runTest(self):
        for pdm in self.pdms():
            initial_rows = 11
            upload_data = []

            for i in range(initial_rows):
                upload_data.append({"bulk1": i, "bulk2": i * 1000})

            upload_df = pd.DataFrame.from_records(upload_data)
            pdm._bulk_load(upload_df, "fn4_bulk_load_test", max_batch=10)

            (final_rows,) = pdm.session.query(func.count(BulkLoadTest.blt_int_id)).one()
            self.assertEqual(initial_rows, final_rows)

            with self.assertRaises(ValueError):
                pdm._bulk_load(upload_df, "nonexistant_table_name", max_batch=10)

            with self.assertRaises(TypeError):
                pdm._bulk_load(None, "fn4_bulk_load_test", max_batch=10)

            with self.assertRaises(RDBMSError):
                empty_df = pd.DataFrame.from_records([])
                pdm._bulk_load(empty_df, "fn4_bulk_load_test", max_batch=10)


class Test_Server_Monitoring_0(Test_Database):
    """adds server monitoring info"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.server_monitoring_store(message="one")

            res = pdm.recent_server_monitoring(100)

            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))


class Test_Server_Monitoring_1(Test_Database):
    """tests recovery of database monitoring info"""

    def runTest(self):
        for pdm in self.pdms():
            res1 = pdm.recent_database_monitoring(100)

            db_summary = pdm.summarise_stored_items()
            pdm.server_monitoring_store(
                what="dbManager", message="Repacking", guid="-", content=db_summary
            )
            res2 = pdm.recent_database_monitoring(100)
            self.assertEqual(
                res1,
                {"latest_stats": {"storage_ratio": 1}, "recompression_data": False},
            )
            self.assertFalse(res2["recompression_data"])
            self.assertEqual(res2["latest_stats"]["storage_ratio"], 1)

            json.dumps(res2)  # should succeed


class Test_Server_Monitoring_2(Test_Database):
    """adds server monitoring info"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.server_monitoring_store(message="one")
            pdm.server_monitoring_store(message="two")
            pdm.server_monitoring_store(message="three")

            res = pdm.recent_server_monitoring(0)
            self.assertEqual(len(res), 0)
            self.assertTrue(isinstance(res, list))

            res = pdm.recent_server_monitoring(1)
            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))

            res = pdm.recent_server_monitoring(3)
            self.assertEqual(len(res), 3)
            self.assertTrue(isinstance(res, list))

            res = pdm.recent_server_monitoring(5)
            self.assertEqual(len(res), 3)
            self.assertTrue(isinstance(res, list))

            with self.assertRaises(TypeError):
                res = pdm.server_monitoring_store(content=42)  # type: ignore

            with self.assertRaises(ValueError):
                res = pdm.recent_server_monitoring(-1)

            with self.assertRaises(TypeError):
                res = pdm.recent_server_monitoring("thing")  # type: ignore


class Test_Server_Monitoring_3(Test_Database):
    """checks whether server_monitoring_min_interval_msec control works"""

    def runTest(self):
        # no logging for within 1 secs of another event
        for pdm in self.pdms(server_monitoring_min_interval_msec=1000):
            retVal = pdm.server_monitoring_store(message="one")  # should insert
            self.assertEqual(retVal, True)
            res = pdm.recent_server_monitoring(100)
            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))

            retVal = pdm.server_monitoring_store(message="two")  # should not insert
            self.assertEqual(retVal, False)
            res = pdm.recent_server_monitoring(100)
            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))

            time.sleep(2)  # seconds
            retVal = pdm.server_monitoring_store(message="three")  # should insert
            self.assertEqual(retVal, True)
            res = pdm.recent_server_monitoring(100)
            self.assertEqual(len(res), 2)
            self.assertTrue(isinstance(res, list))


class Test_Server_Monitoring_4(Test_Database):
    """checks whether delete_server_monitoring_entries"""

    def runTest(self):
        for pdm in self.pdms():
            retVal = pdm.server_monitoring_store(message="one")  # should insert
            self.assertEqual(retVal, True)
            res = pdm.recent_server_monitoring(100)
            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))
            pdm.delete_server_monitoring_entries(1)
            res = pdm.recent_server_monitoring(100)
            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))

            time.sleep(2)  # seconds

            pdm.delete_server_monitoring_entries(1)
            res = pdm.recent_server_monitoring(100)
            self.assertEqual(len(res), 0)
            self.assertTrue(isinstance(res, list))


class Test_Server_Monitoring_5(Test_Database):
    """adds server monitoring info"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.server_monitoring_store(content={"x": "true"})
            pdm.server_monitoring_store(content={"x": "false"})

            res = pdm.recent_server_monitoring()
            self.assertEqual(len(res), 2)
            self.assertTrue(isinstance(res, list))

            res = pdm.recent_server_monitoring(
                selection_field="x", selection_string="true"
            )
            self.assertEqual(len(res), 1)
            self.assertTrue(isinstance(res, list))


class Test_SeqMeta_singleton(Test_Database):
    """tests guid2neighboursOf"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.guid2neighbour_add_links(
                "srcguid",
                {
                    "guid1": {"dist": 12},
                    "guid2": {"dist": 0},
                    "guid3": {"dist": 3},
                    "guid4": {"dist": 4},
                    "guid5": {"dist": 5},
                },
            )
            res1 = pdm.guid2neighbours("srcguid", returned_format=1)

            self.assertEqual(5, len(res1["neighbours"]))
            singletons = pdm.singletons(method="exact")
            self.assertEqual(len(singletons.index), 0)
            singletons = pdm.singletons(method="approximate")
            self.assertEqual(len(singletons.index), 0)


class Test_SeqMeta_guid2neighbour_8(Test_Database):
    """tests guid2neighboursOf"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.guid2neighbour_add_links(
                "srcguid",
                {
                    "guid1": {"dist": 12},
                    "guid2": {"dist": 0},
                    "guid3": {"dist": 3},
                    "guid4": {"dist": 4},
                    "guid5": {"dist": 5},
                },
            )

            res1 = pdm.guid2neighbours("srcguid", returned_format=1)
            self.assertEqual(5, len(res1["neighbours"]))
            with self.assertRaises(NotImplementedError):
                res2 = pdm.guid2neighbours("srcguid", returned_format=2)
                self.assertTrue(res2 is not None)

            res3 = pdm.guid2neighbours("srcguid", returned_format=3)
            self.assertEqual(5, len(res3["neighbours"]))
            res4 = pdm.guid2neighbours("srcguid", returned_format=4)
            self.assertEqual(5, len(res4["neighbours"]))

            with self.assertRaises(ValueError):
                pdm.guid2neighbours("srcguid", returned_format=5)


class Test_SeqMeta_guid2neighbour_7(Test_Database):
    """tests guid2neighboursOf"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.guid2neighbour_add_links(
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
                res1 = pdm.guid2neighbours("srcguid", returned_format=2)
            res1 = pdm.guid2neighbours("srcguid", returned_format=1)
            self.assertEqual(5, len(res1["neighbours"]))
            pdm.guid2neighbour_repack(guid="srcguid")
            res2 = pdm.guid2neighbours("srcguid", returned_format=1)
            self.assertEqual(5, len(res2["neighbours"]))


class Test_SeqMeta_guid_exists_1(Test_Database):
    """tests insert of new data item and existence check"""

    def runTest(self):
        for pdm in self.pdms():
            # test there is no 'test' item; insert, and confirm insert
            guid = "sequence1"
            namespace = "ns"
            payload = {"one": 1, "two": 2}
            pdm.refcompressedseq_store(guid, {"seq": "ACTG"})
            pdm.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
            res = pdm.guid_exists(guid)
            self.assertEqual(res, True)
            res = pdm.guid_exists("missing")
            self.assertEqual(res, False)

            with self.assertRaises(TypeError):
                pdm.guid_annotate(guid=guid, nameSpace=namespace, annotDict=42)


class Test_SeqMeta_guid_valid_1(Test_Database):
    """tests insert of new data item and validity check"""

    def runTest(self):
        for pdm in self.pdms():
            guid = "valid"
            namespace = "DNAQuality"
            payload = {"invalid": 0}
            pdm.refcompressedseq_store(guid, {"seq": "ACTG"})

            pdm.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)
            guid = "invalid"
            pdm.refcompressedseq_store(guid, {"seq": "ACTG"})

            namespace = "DNAQuality"
            payload = {"invalid": 1}
            pdm.guid_annotate(guid=guid, nameSpace=namespace, annotDict=payload)

            guid = "missing"
            namespace = "DNAQuality"
            payload = {"N": 1}
            with self.assertRaises(RDBMSError):
                pdm.guid_annotate(
                    guid=guid, nameSpace=namespace, annotDict=payload
                )  # does not exist

            res = pdm.guid_valid("valid")
            self.assertEqual(res, 0)
            res = pdm.guid_valid("invalid")
            self.assertEqual(res, 1)
            res = pdm.guid_valid("missing")
            self.assertEqual(res, -1)

            res = pdm.guids_valid()
            self.assertEqual(
                res,
                set(
                    [
                        "valid",
                    ]
                ),
            )
            res = pdm.guids_invalid()
            self.assertEqual(res, set(["invalid"]))


class Test_SeqMeta_init(Test_Database):
    """tests database creation"""

    def runTest(self):
        for pdm in self.pdms():
            self.assertTrue(pdm.first_run() is True)
            res = pdm.config_read("preComparer")
            self.assertEqual(res, None)

            pdm.config_store("config", {"item": 1})
            self.assertTrue(pdm.first_run() is False)
            res = pdm.config_read("config")
            self.assertEqual(res, {"_id": "config", "item": 1})

            pdm.config_store("preComparer", {"item": 2})
            res = pdm.config_read("config")
            self.assertEqual(res, {"_id": "config", "item": 1})
            res = pdm.config_read("preComparer")
            self.assertEqual(res, {"_id": "preComparer", "item": 2})

            pdm.config_store("preComparer", {"item": 3})
            res = pdm.config_read("config")
            self.assertEqual(res, {"_id": "config", "item": 1})
            res = pdm.config_read("preComparer")
            self.assertEqual(res, {"_id": "preComparer", "item": 3})

            with self.assertRaises(TypeError):
                pdm.config_store("", 42)


class Test_SeqMeta_guid_quality_check_1(Test_Database):
    def runTest(self):
        """tests return of sequences and their qualities"""
        # set up nucleic acid object
        na = NucleicAcid()

        for pdm in self.pdms():
            na.examine("ACGTACGTNN")  # 20% bad
            pdm.refcompressedseq_store("g1", {"seq": "ACGTACGTNN"})
            pdm.guid_annotate(
                guid="g1", nameSpace="DNAQuality", annotDict=na.composition
            )
            na.examine("ACGTACNNNN")  # 40% bad
            pdm.refcompressedseq_store("g2", {"seq": "ACGTACNNNN"})
            pdm.guid_annotate(
                guid="g2", nameSpace="DNAQuality", annotDict=na.composition
            )
            na.examine("ACGTNNNNNN")  # 60% bad
            pdm.refcompressedseq_store("g3", {"seq": "ACGNNNNNNN"})
            pdm.guid_annotate(
                guid="g3", nameSpace="DNAQuality", annotDict=na.composition
            )

            r1 = pdm.guid_quality_check("g1", 0.80)  # valid
            r2 = pdm.guid_quality_check("g2", 0.80)  # invalid
            r3 = pdm.guid_quality_check("g3", 0.80)  # invalid
            r4 = pdm.guid_quality_check("g4", 0.80)  # invalid; does not exist.

            self.assertEqual(r1, True)
            self.assertEqual(r2, False)
            self.assertEqual(r3, False)
            self.assertEqual(r4, None)

            with self.assertRaises(TypeError):
                pdm.guid_quality_check("g1", "0.80")

            with self.assertRaises(TypeError):
                pdm.guid_quality_check(1, 0.80)

            resDict = pdm.guid2quality(None)  # restrict to nothing - return all
            self.assertIsNotNone(resDict)
            self.assertEqual(resDict["g1"], 0.80)
            self.assertEqual(resDict["g2"], 0.60)
            self.assertEqual(resDict["g3"], 0.40)


class Test_SeqMeta_guid2quality2(Test_Database):
    def runTest(self):
        """tests return of sequences and their qualities"""
        # set up nucleic acid object

        na = NucleicAcid()

        for pdm in self.pdms():
            na.examine("ACGTACGTNN")  # 20% bad
            pdm.refcompressedseq_store("g1", {"seq": "ACGTACGTNN"})
            pdm.guid_annotate(
                guid="g1", nameSpace="DNAQuality", annotDict=na.composition
            )
            na.examine("ACGTACNNNN")  # 40% bad
            pdm.refcompressedseq_store("g2", {"seq": "ACGTACNNNN"})
            pdm.guid_annotate(
                guid="g2", nameSpace="DNAQuality", annotDict=na.composition
            )
            na.examine("ACGTNNNNNN")  # 60% bad
            pdm.refcompressedseq_store("g3", {"seq": "ACGNNNNNNN"})
            pdm.guid_annotate(
                guid="g3", nameSpace="DNAQuality", annotDict=na.composition
            )
            r1 = pdm.guid_quality_check("g1", 0.80)  # valid
            r2 = pdm.guid_quality_check("g2", 0.80)  # invalid
            r3 = pdm.guid_quality_check("g3", 0.80)  # invalid

            self.assertEqual(r1, True)
            self.assertEqual(r2, False)
            self.assertEqual(r3, False)  # check the db insert works

            resDict = pdm.guid2quality(["g1", "g2", "g3"])
            self.assertTrue(resDict is not None)
            assert resDict is not None  # for typing purposes
            self.assertEqual(resDict["g1"], 0.80)
            self.assertEqual(resDict["g2"], 0.60)
            self.assertEqual(resDict["g3"], 0.40)


class Test_SeqMeta_Base1(Test_Database):
    """initialise FN persistence and adds data"""

    def pdms(self):
        dna = NucleicAcid()
        seqs = {"guid1": "ACGT", "guid2": "NACT", "guid3": "TTTT", "guid4": "NNNN"}

        for pdm in super().pdms():
            for guid, seq in seqs.items():
                dna.examine(seq)
                pdm.refcompressedseq_store(
                    guid, {"seq": seq}
                )  # in real application, seq would be a dictionary of reference compressed data
                pdm.guid_annotate(
                    guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
                )

            res = set()
            for (x,) in pdm.session.query(RefCompressedSeq.sequence_id).all():
                res.add(x)
            self.assertTrue(x, set(seqs.keys()))
            yield pdm


class Test_SeqMeta_Base1t(Test_Database):
    """initialise FN persistence and adds data, 0.1 secs apart.
    Used for testing queries examining order of recovery of samples."""

    def pdms(self):
        dna = NucleicAcid()
        self.seqs = {"guid1": "ACGT", "guid2": "NACT", "guid3": "TTTT", "guid4": "NNNN"}

        for pdm in super().pdms():
            for guid, seq in self.seqs.items():
                time.sleep(0.1)
                dna.examine(seq)
                pdm.refcompressedseq_store(
                    guid, {"seq": seq}
                )  # in real application, seq would be a dictionary of reference compressed data
                pdm.guid_annotate(
                    guid=guid, nameSpace="DNAQuality", annotDict=dna.composition
                )
            res = set()
            for (x,) in pdm.session.query(RefCompressedSeq.sequence_id).all():
                res.add(x)
            self.assertTrue(x, set(self.seqs.keys()))
            yield pdm


class Test_SeqMeta_guid2ExaminationDateTime(Test_SeqMeta_Base1):
    """recovering guids and examination times;"""

    def runTest(self):
        for pdm in self.pdms():
            res = pdm.guid2ExaminationDateTime()
            self.assertIsNotNone(res)
            expected = 4
            self.assertEqual(len(res.keys()), expected)


class Test_SeqMeta_guid2ExaminationDateTime_order(Test_SeqMeta_Base1t):
    """tests guid2ExaminationDateTime"""

    def runTest(self):
        for pdm in self.pdms():
            res = pdm.guid2ExaminationDateTime()
            self.assertIsNotNone(res)
            expected = 4
            self.assertEqual(len(res.keys()), expected)

            # check that the sample were added in order, with increasing examination times.
            previous_addition_time = None
            for i, guid in enumerate(sorted(self.seqs.keys())):  # the order added
                if i > 0:
                    self.assertGreater(res[guid], previous_addition_time)
                previous_addition_time = res[guid]


class Test_SeqMeta_guid_examination_time(Test_SeqMeta_Base1t):
    """tests guid_examination_time()"""

    def runTest(self):
        for pdm in self.pdms():
            res = pdm.guid2ExaminationDateTime()
            self.assertIsNotNone(res)
            assert res is not None  # for typing purposes
            expected = 4
            self.assertEqual(len(res.keys()), expected)

            # check that the sample were added in order, with increasing examination times.
            previous_addition_time = None
            for i, guid in enumerate(sorted(self.seqs.keys())):  # the order added
                this_examination_time = pdm.guid_examination_time(guid)
                if i > 0:
                    self.assertGreater(this_examination_time, previous_addition_time)
                previous_addition_time = this_examination_time

            this_examination_time = pdm.guid_examination_time("missing-guid")
            self.assertIsNone(this_examination_time)


class Test_SeqMeta_guid_considered_after(Test_SeqMeta_Base1t):
    """recovering guids and examination times;"""

    def runTest(self):
        for pdm in self.pdms():
            res = pdm.guid2ExaminationDateTime()
            self.assertIsNotNone(res)
            assert res is not None  # for typing purposes
            expected = 4
            self.assertEqual(len(res.keys()), expected)

            # check that the sample were added in order, with increasing examination times.
            for i, guid in enumerate(sorted(self.seqs.keys())):  # the order added
                this_examination_time = pdm.guid_examination_time(guid)
                self.assertIsNotNone(this_examination_time)
                assert this_examination_time is not None  # for typing purposes
                res2 = pdm.guids_considered_after(this_examination_time)
                self.assertEqual(
                    len(res2), 3 - i
                )  # with guid1, we expect three; with guid2, we expect 2; etc

                res2 = pdm.guids_considered_after_guid(guid)
                self.assertEqual(
                    len(res2), 3 - i
                )  # with guid1, we expect three; with guid2, we expect 2; etc

                with self.assertRaises(ValueError):
                    pdm.guids_considered_after_guid("not_a_guid")

            with self.assertRaises(TypeError):
                pdm.guids_considered_after(42)


class Test_SeqMeta_propACTG_filteredSequenceGuids(Test_SeqMeta_Base1):
    """recovered guids filtered by the propACTG criterion"""

    def runTest(self):
        for pdm in self.pdms():
            n = 0
            for guid in pdm.guid2propACTG_filtered(cutoff=0.85):
                n += 1
            expected = 2
            self.assertEqual(n, expected)


class Test_SeqMeta_allAnnotations(Test_SeqMeta_Base1):
    """tests recovery of all annoations"""

    def runTest(self):
        for pdm in self.pdms():
            df = pdm.guid_annotations()
            self.assertIsNotNone(df)
            assert df is not None  # for typing purposes
            self.assertEqual(len(df.keys()), 4)


class Test_SeqMeta_oneAnnotation(Test_SeqMeta_Base1):
    """tests recovery of one annotations"""

    def runTest(self):
        for pdm in self.pdms():
            df = pdm.guid_annotation("guid3")
            self.assertIsNotNone(df)
            assert df is not None  # for typing purposes
            self.assertEqual(len(df.keys()), 1)
            df = pdm.guid_annotation("missing")
            self.assertIsNotNone(df)
            assert df is not None  # for typing purposes
            self.assertEqual(len(df.keys()), 0)


class Test_Clusters(Test_Database):
    """tests saving and recovery of dictionaries to Clusters"""

    def runTest(self):
        for pdm in self.pdms():

            # there aren't any clusters initially
            self.assertIsNone(pdm.cluster_latest_version("cl1"))
            self.assertEqual(0, len(pdm.cluster_versions("cl1")))

            # we store a dictionary as cluster 1 (cl1)
            payload1 = {"one": 1, "two": 2}
            x = pdm.cluster_store("cl1", payload1)  # returns the cluster number

            # and another s cluster2
            payload1b = {"2_one": 1, "2_two": 2}
            y = pdm.cluster_store("cl2", payload1b)  # returns the cluster number
            self.assertIsNotNone(y)
            self.assertIsNotNone(pdm.cluster_latest_version("cl1"))

            # we measure the id of the latest version of cluster 1.
            clv = pdm.cluster_latest_version("cl1")
            self.assertEqual(x, clv)  # it's the same number we were given initially
            self.assertEqual(
                1, len(pdm.cluster_versions("cl1"))
            )  # there is only one version

            # this should not find any newer versions than clv; it should return none
            self.assertIsNone(pdm.cluster_read_update("cl1", clv))

            # should get back what we put in
            payload2 = pdm.cluster_read("cl1")
            self.assertEqual(payload1, payload2)

            # now we load an update.
            payload3 = {"one": 10, "two": 20}
            x2 = pdm.cluster_store(
                "cl1", payload3
            )  # this is now the latest version of cl1
            self.assertTrue(x2 > x)  # will have a higher number

            payload4 = {"2-one": 10, "2-two": 20}
            y2 = pdm.cluster_store(
                "cl2", payload4
            )  # this is now the latest version of cl2
            self.assertTrue(y2 > y)  # will have a higher number

            self.assertEqual(
                pdm.cluster_keys(), ["cl1", "cl2"]
            )  # there should be two clusters in the dictionary
            self.assertEqual(
                pdm.cluster_keys(clustering_name="cl1"), ["cl1"]
            )  # selecting one should yield only one

            self.assertEqual(
                2, len(pdm.cluster_versions("cl1"))
            )  # there are two different cluster 1 versions
            self.assertEqual(
                x2, pdm.cluster_latest_version("cl1")
            )  # should return the latest version
            self.assertNotEqual(x2, x)

            self.assertNotEqual(clv, pdm.cluster_latest_version("cl1"))  # fails
            self.assertIsNotNone(pdm.cluster_read_update("cl1", clv))

            payload5 = pdm.cluster_read("cl1")
            self.assertEqual(payload5, payload3)

            pdm.cluster_delete_legacy_by_key("cl1")

            self.assertEqual(1, len(pdm.cluster_versions("cl1")))

            pdm.cluster_delete_all("cl1")

            self.assertEqual(0, len(pdm.cluster_versions("cl1")))
            self.assertEqual(2, len(pdm.cluster_versions("cl2")))
            pdm.cluster_delete_legacy_by_key("cl2")
            self.assertEqual(1, len(pdm.cluster_versions("cl2")))

            pdm.cluster_store("cl1", {"datum": 42})
            pdm.cluster_store("cl1", {"datum": 42})
            pdm.cluster_store("cl2", {"datum": 42})
            pdm.cluster_store("other", {"datum": 42})
            pdm.cluster_store("other", {"datum": 42})

            self.assertEqual(2, len(pdm.cluster_versions("cl1")))
            self.assertEqual(2, len(pdm.cluster_versions("cl2")))
            self.assertEqual(2, len(pdm.cluster_versions("other")))
            pdm.cluster_delete_legacy("cl")
            self.assertEqual(1, len(pdm.cluster_versions("cl1")))
            self.assertEqual(1, len(pdm.cluster_versions("cl2")))
            self.assertEqual(2, len(pdm.cluster_versions("other")))

            with self.assertRaises(TypeError):
                pdm.cluster_store("", 42)


class Test_Tree(Test_Database):
    """tests saving and recovery of dictionaries to Tree"""

    def runTest(self):
        for pdm in self.pdms():
            payload1 = {"one": 1, "two": 2}
            pdm.tree_store(tree_token="tree1", tree=payload1)
            payload2 = pdm.tree_read(tree_token="tree1")
            self.assertEqual(payload1, payload2)
            pdm.tree_delete(tree_token="tree1")
            payload3 = pdm.tree_read(tree_token="tree1")
            self.assertIsNone(payload3)

            payload1 = {"one": 1, "two": 2}
            pdm.tree_store(tree_token="tree1", tree=payload1)
            payload2 = {"one": 3, "two": 4}
            pdm.tree_store(tree_token="tree2", tree=payload2)
            self.assertEqual(2, len(pdm.tree_stored_ids()))
            pdm.tree_delete_unless_whitelisted(whitelist=["tree1", "tree2"])
            self.assertEqual(2, len(pdm.tree_stored_ids()))
            pdm.tree_delete_unless_whitelisted(whitelist=["tree1"])
            self.assertEqual(1, len(pdm.tree_stored_ids()))

            with self.assertRaises(TypeError):
                pdm.tree_store("", 42)


class Test_MSA(Test_Database):
    """tests saving and recovery of dictionaries to MSA"""

    def runTest(self):
        for pdm in self.pdms():
            payload1 = {"one": 1, "two": 2}
            pdm.msa_store(msa_token="msa1", msa=payload1)
            payload2 = pdm.msa_read(msa_token="msa1")
            self.assertEqual(payload1, payload2)
            pdm.msa_delete(msa_token="msa1")
            payload3 = pdm.msa_read(msa_token="msa1")
            self.assertIsNone(payload3)

            payload1 = {"one": 1, "two": 2}
            pdm.msa_store(msa_token="msa1", msa=payload1)
            payload2 = {"one": 3, "two": 4}
            pdm.msa_store(msa_token="msa2", msa=payload2)
            self.assertEqual(2, len(pdm.msa_stored_ids()))
            pdm.msa_delete_unless_whitelisted(whitelist=["msa1", "msa2"])
            self.assertEqual(2, len(pdm.msa_stored_ids()))
            pdm.msa_delete_unless_whitelisted(whitelist=["msa1"])
            self.assertEqual(1, len(pdm.msa_stored_ids()))

            with self.assertRaises(TypeError):
                pdm.msa_store("", 42)


class Test_Monitor(Test_Database):
    """tests saving and recovery of strings to monitor"""

    def runTest(self):
        for pdm in self.pdms():
            payload1 = "line1"
            pdm.monitor_store("r1", payload1)
            payload2 = pdm.monitor_read("r1")
            self.assertEqual(payload1, payload2)
            payload3 = pdm.monitor_read("nil")
            self.assertIsNone(payload3)

            with self.assertRaises(TypeError):
                pdm.monitor_store("", 42)


class test_Raise_error(Test_Database):
    """tests raise_error"""

    def runTest(self):
        for pdm in self.pdms():
            with self.assertRaises(ZeroDivisionError):
                pdm.raise_error("token")


class Test_summarise_stored_items(Test_Database):
    """adds server monitoring info"""

    def runTest(self):
        for pdm in self.pdms():
            self.assertIsNotNone(pdm.summarise_stored_items())


class Test_rotate_log(Test_Database):
    """dummy test for coverage, this is a NOP"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.rotate_log()


class Test_no_progressbar(Test_Database):
    """dummy test for coverage, this is not easily testable"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.no_progressbar()


class Test_connect(Test_Database):
    """dummy test for coverage, this is a NOP"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.connect()


# @unittest.skip("On oracle, connection drops, reason unknown")
class Test_delete_existing_data(Test_Database):
    """check that all data is deleted"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.refcompressedseq_store("guid", {"datum": 42})
            pdm.config_store("config", {"datum": 42})
            pdm.server_monitoring_store()
            pdm.monitor_store("monitor", "<html></html>")
            pdm.msa_store("msa", {"datum": 42})
            pdm.tree_store("tree", {"datum": 42})
            pdm.cluster_store("cluster", {"datum": 42})

            pdm._delete_existing_data()

            self.assertIsNone(pdm.config_read("config"))
            self.assertEqual(len(pdm.recent_server_monitoring()), 0)
            self.assertIsNone(pdm.monitor_read("monitor"))
            self.assertIsNone(pdm.msa_read("msa"))
            self.assertIsNone(pdm.tree_read("tree"))
            self.assertIsNone(pdm.cluster_read("cluster"))
            self.assertIsNone(pdm.refcompressedseq_store("guid", {"datum": 42}))


class Test_delete_existing_clustering_data(Test_Database):
    """check that all clustering data is deleted"""

    def runTest(self):
        for pdm in self.pdms():
            pdm.msa_store("msa", {"datum": 42})
            pdm.tree_store("tree", {"datum": 42})
            pdm.cluster_store("cluster", {"datum": 42})

            pdm._delete_existing_clustering_data()

            self.assertIsNone(pdm.msa_read("msa"))
            self.assertIsNone(pdm.tree_read("tree"))
            self.assertIsNone(pdm.cluster_read("cluster"))


class Test_memory_usage(Test_Database):
    """get memory usage"""

    def runTest(self):
        for pdm in self.pdms():
            res = pdm.memory_usage()
            self.assertIsNotNone(res["server|mstat|total"])
            self.assertIsNotNone(res["server|mstat|available"])
            self.assertIsNotNone(res["server|mstat|used"])
            self.assertIsNotNone(res["server|mstat|free"])
            self.assertTrue(0 <= res["server|mstat|percent"] <= 100)


class Test_refcompressedseq_store(Test_Database):
    """this is mostly tested by the edges, just test the edge cases"""

    def runTest(self):
        for pdm in self.pdms():
            payload = {"datum": 42}
            pdm.refcompressedseq_store("guid", payload)
            self.assertEqual(pdm.refcompressedsequence_read("guid"), payload)
            pdm.refcompressedseq_store("guid2", {})
            self.assertEqual(pdm.guids(), {"guid", "guid2"})

            self.assertIsNone(pdm.refcompressedsequence_read("not_a_guid"))

            with self.assertRaises(TypeError):
                pdm.refcompressedseq_store("", 42)


class Test_guid2items(Test_Database):
    """this is mostly tested by other test cases, just test the edge cases"""

    def runTest(self):

        for pdm in self.pdms():
            pdm.refcompressedseq_store("guid1", {"seq": "ACTG"})
            pdm.refcompressedseq_store("guid2", {"seq": "ACTG"})
            pdm.guid_annotate("guid1", "ns1", {"datum": 1})
            pdm.guid_annotate("guid1", "ns2", {"datum": 2})
            pdm.guid_annotate("guid2", "ns2", {"datum": 3})
            self.assertEqual(
                pdm.guid2items(None, None),
                {
                    "guid1": {"ns1": {"datum": 1}, "ns2": {"datum": 2}},
                    "guid2": {"ns2": {"datum": 3}},
                },
            )
            self.assertEqual(
                pdm.guid2items(["guid1"], None),
                {"guid1": {"ns1": {"datum": 1}, "ns2": {"datum": 2}}},
            )
            self.assertEqual(
                pdm.guid2items(None, ["ns1"]),
                {
                    "guid1": {"ns1": {"datum": 1}},
                    "guid2": {},
                },
            )
            self.assertEqual(
                pdm.guid2items(None, ["ns2"]),
                {
                    "guid1": {"ns2": {"datum": 2}},
                    "guid2": {"ns2": {"datum": 3}},
                },
            )

            self.assertEqual(
                pdm.guid2item(None, "ns2", "datum"), {"guid1": 2, "guid2": 3}
            )
            self.assertEqual(pdm.guid2item(["guid1"], "ns1", "datum"), {"guid1": 1})
