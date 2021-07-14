""" tests pcadb.py - software to do store output from PCA

    to invoke, either run automatically with pytest or do
    pipenv run python3 -m unittest test/test_rdbmsstore.py

"""
import os
import json
import unittest
import pandas as pd
from sqlalchemy import func
from findn.rdbmsstore import fn3persistence, BulkLoadTest


class Test_Database(unittest.TestCase):
    """establishes database connection strings for cross-database testing.
    Currently tests OCI (if relevant environment variables are set) and Sqlite

    To test other databases, such as MySql, add the relevant connection string &
    database name to the dictionary self.engines"""

    def setUp(self):
        self.engines = {}
        self.engines["Sqlite"] = "sqlite://"  # in memory sqlite

        conn_detail_file = None

        # try to read the environment variable 'PCA_CONNECTION_CONFIG_FILE'
        try:
            conn_detail_file = os.environ["PCA_CONNECTION_CONFIG_FILE"]
        except KeyError:
            # doesn't exist; we just run with sqlite, which is the default if engine is None.
            print(
                "No environment variable PCA_CONNECTION_CONFIG_FILE found.  Testing with sqlite only."
            )
        if conn_detail_file is None:
            print(
                "No environment variable PCA_CONNECTION_CONFIG_FILE found.  Testing with sqlite only."
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


class Test_create_database_1(Test_Database):
    """tests creating the database and internal functions dropping tables"""

    def runTest(self):
        expected_tables = set(
            [
                "fn4_bulk_load_test",
                "refcompressedseq",
                "edge",
                "cluster",
                "monitor",
                "msa",
                "tree",
            ]
        )
        for engine in self.engines.keys():
            print(engine, "#1")
            pdm = fn3persistence(connection_config=self.engines[engine], debug=2)
            n1 = pdm._table_names()
            pdm._drop_existing_tables()
            n2 = pdm._table_names()
            self.assertEqual(len(expected_tables.intersection(set(n1))), 7)
            self.assertEqual(len(expected_tables.intersection(set(n2))), 0)
           
            pdm = fn3persistence(connection_config=self.engines[engine], debug=0)
            pdm = fn3persistence(connection_config=self.engines[engine], debug=2)


class Test_oracle_bulk_upload_2(Test_Database):
    """tests bulk upload with small number of samples"""

    def runTest(self):
        for engine in self.engines.keys():
            print(engine, "#2")

            pdm = fn3persistence(connection_config=self.engines[engine], debug=True)
            initial_rows = 11
            upload_data = []

            for i in range(initial_rows):
                upload_data.append({"bulk1": i, "bulk2": i * 1000})

            upload_df = pd.DataFrame.from_records(upload_data)
            pdm._bulk_load(upload_df, "fn4_bulk_load_test", max_batch=10)

            (final_rows,) = pdm.session.query(func.count(BulkLoadTest.blt_int_id)).one()
            self.assertEqual(initial_rows, final_rows)
