""" tests pcadb.py - software to do store output from PCA

"""
import os
import json
import unittest
import pandas as pd
import pickle
from sqlalchemy import func
from pca.pca import VariationModel
from pca.pcadb import ContingencyTable, BulkLoadTest
from pca.pcadb import PCADatabaseManager, FeatureAssociation


class Test_PCA_Database(unittest.TestCase):
    """establishes database connection strings for cross-database testing.
    Currently tests OCI (if relevant environment variables are set) and Sqlite

    To test other databases, such as MySql, add the relevant connection string &
    database name to the dictionary self.engines"""

    def setUp(self):
        self.engines = {} # {'None_specified':None}
        
        # try to read the environment variable 'PCA_CONNECTION_CONFIG_FILE'
        try:
            conn_detail_file = os.environ['PCA_CONNECTION_CONFIG_FILE']
        except KeyError:
            # doesn't exist; we just run with sqlite, which is the default if engine is None.
            print("No environment variable PCA_CONNECTION_CONFIG_FILE found.  Testing with sqlite only.")
        if conn_detail_file is None:
            print("No environment variable PCA_CONNECTION_CONFIG_FILE found.  Testing with sqlite only.")
            return
        if not os.path.exists(conn_detail_file):
            raise FileNotFoundError("Connection file specified but not found: {0}".format(conn_detail_file))
        # try to read the config file
        with open(conn_detail_file,'rt') as f:
            conn_detail = json.load(f)
            for key in conn_detail.keys():
                if key.startswith('unittest_ora'):
                    self.engines[key] = key

        return
        oracle_credential_location = os.getenv("TNS_ADMIN")
        oracle_engine_name = os.getenv("OCI_ENGINE_NAME")
        if len(oracle_credential_location) > 0 and len(oracle_engine_name) > 0:
            print("Using Oracle credentials at {0}".format(oracle_credential_location))
            self.engines["Oracle"] = oracle_engine_name  # Debug
        self.engines["Sqlite"] = "sqlite:///test.sqlite"


##@unittest.skip(reason="too slow")
class Test_create_database_1(Test_PCA_Database):
    """tests creating the database and internal functions dropping tables"""

    def runTest(self):
        for engine in self.engines.keys():
            print(engine, "#1")
            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            n1 = len(pdm._table_names())
            pdm._drop_existing_tables()
            n2 = len(pdm._table_names())
            self.assertEqual(n1, 18)
            self.assertEqual(n2, 0)

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=False)
            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)

            # no builds
            res = pdm.latest_build_int_id()
            self.assertIsNone(res)


@unittest.skip(reason="too slow")
class Test_create_database_2(Test_PCA_Database):
    """tests storing the results of a VariantModel"""

    def runTest(self):

        # load a variation model for testiong
        inputfile = "testdata/pca/vm.pickle"
        with open(inputfile, "rb") as f:
            vm = pickle.load(f)
        self.assertIsInstance(vm, VariationModel)

        for engine in self.engines.keys():
            print(engine, "#2")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            b1 = pdm.latest_build_int_id()
            self.assertIsNone(b1)

            pdm.store_variation_model(vm)

            b2 = pdm.latest_build_int_id()
            self.assertIsNotNone(b2)
            self.assertEqual(b2, 1)


@unittest.skip(reason="too slow")
class Test_load_clinical_data_3(Test_PCA_Database):
    """tests storing clinical data from COG-UK in the data model"""

    def runTest(self):

        # load a variation model for testiong

        for engine in self.engines.keys():
            print(engine, "#3")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)

            x1 = pdm.existing_sample_ids_in_clinical_metadata()
            x3 = pdm.store_cog_metadata(cogfile="testdata/pca/cog_metadata_test.csv")
            b3 = pdm.latest_build_int_id()
            self.assertIsNone(b3)

            x2 = pdm.existing_sample_ids_in_clinical_metadata()
            self.assertEqual(len(x1), 0)
            self.assertEqual(len(x2), 99)
            self.assertEqual(len(x2), x3)

            # try to add again
            x4 = pdm.store_cog_metadata(cogfile="testdata/pca/cog_metadata_test.csv")
            self.assertEqual(x4, 0)
            b4 = pdm.latest_build_int_id()
            self.assertIsNone(b4)

            # try to add again,  with replace_all switched on.
            x5 = pdm.store_cog_metadata(
                cogfile="testdata/pca/cog_metadata_test.csv", replace_all=True
            )
            self.assertEqual(x5, 99)
            b5 = pdm.latest_build_int_id()
            self.assertIsNone(b5)

            sf = pdm.sequence_features(is_lineage=True)
            l1 = len(sf)
            sf = pdm.sequence_features(is_lineage=False)
            l2 = len(sf)
            sf = pdm.sequence_features(is_lineage=None)
            l3 = len(sf)
            self.assertEqual(l1 + l2, l3)


@unittest.skip(reason="not necessary")
class Test_Contingency_table(unittest.TestCase):
    """tests the ContigencyTable class, which analyses 2x2 tables"""

    def runTest(self):
        """
        parameters are :(build_int_id, pc_cat, sequencefeature, a, b, c, d):
        expects as parameters count data for a 2x2 contingency table

        Feature    Present       Absent
        PC-CAT Y      a            c          a+c
               N      b            d          b+d
                    a+b           c+d       a+b+c+d

        tc_int_id and sf_int_id are integer tokens used to identify the results

        computes a G-test (LR test) and an estimated odds ratio.
        """
        ct = ContingencyTable(1, "1_0", "lineage:1", 10, 10, 10, 10)
        res = ct.featureassociation()
        self.assertIsInstance(res, FeatureAssociation)

        ct = ContingencyTable(1, "1_0", "lineage:1", 0, 0, 0, 0)
        res = ct.featureassociation()
        self.assertIsInstance(res, FeatureAssociation)

        ct = ContingencyTable(1, "1_0", "lineage:1", 1000, 0, 0, 1000)
        res = ct.featureassociation()
        self.assertIsInstance(res, FeatureAssociation)


@unittest.skip(reason ='too slow')
class Test_create_contingency_tables_4(Test_PCA_Database):
    """tests creation & storage of 2x2 contingency tables"""

    def runTest(self):

        # load a variation model for testiong
        inputfile = "testdata/pca/vm.pickle"
        with open(inputfile, "rb") as f:
            vm = pickle.load(f)

        # write cog metadata from a large set of  samples in the vm into a file
        # modelled = set(vm.model['sample_id'])
        # inputfile = '/data/data/inputfasta/cog_metadata.csv'
        # cog = pd.read_csv(inputfile)
        # cog['sample_id'] = [x[1] for x in cog['sequence_name'].str.split('/')]
        # cog = cog[cog['sample_id'].isin(modelled)]
        # cog.drop(['sample_id'], axis =1, inplace =True)
        # cog.drop(['mutations'], axis =1, inplace =True)
        # cog.to_csv('testdata/pca/cog_metadata_vm.csv', index=False)

        for engine in self.engines.keys():
            print(engine, "#4")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            pdm.store_variation_model(vm)
            pdm.store_cog_metadata(cogfile="testdata/pca/cog_metadata_vm.csv")
            pdm.make_contingency_tables()


@unittest.skip(reason="not routinely necessary ")
class Test_oracle_bulk_upload_5a(Test_PCA_Database):
    """tests bulk upload"""

    def runTest(self):
        for engine in self.engines.keys():
            print(engine, "#5a")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            initial_rows = 10
            upload_data = []

            for i in range(initial_rows):
                upload_data.append({"bulk1": i, "bulk2": i * 1000})

            upload_df = pd.DataFrame.from_records(upload_data)
            pdm._bulk_load(upload_df, "test", max_batch=1)

            (final_rows,) = pdm.session.query(func.count(BulkLoadTest.blt_int_id)).one()
            self.assertEqual(initial_rows, final_rows)


@unittest.skip(reason="not routinely necessary ")
class Test_oracle_bulk_upload_5b(Test_PCA_Database):
    """tests bulk upload with small number of samples"""

    def runTest(self):
        for engine in self.engines.keys():
            print(engine, "#5b")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            initial_rows = 11
            upload_data = []

            for i in range(initial_rows):
                upload_data.append({"bulk1": i, "bulk2": i * 1000})

            upload_df = pd.DataFrame.from_records(upload_data)
            pdm._bulk_load(upload_df, "test", max_batch=10)

            (final_rows,) = pdm.session.query(func.count(BulkLoadTest.blt_int_id)).one()
            self.assertEqual(initial_rows, final_rows)


@unittest.skip(reason="not routinely necessary and slow")
class Test_oracle_bulk_upload_5c(Test_PCA_Database):
    """tests bulk upload with large numbers of samples """

    def runTest(self):
        for engine in self.engines.keys():
            print(engine, "#5c")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            initial_rows = 1000000
            upload_data = []

            for i in range(initial_rows):
                upload_data.append({"bulk1": i, "bulk2": i * 1000})

            upload_df = pd.DataFrame.from_records(upload_data)
            pdm._bulk_load(upload_df, "test")

            (final_rows,) = pdm.session.query(func.count(BulkLoadTest.blt_int_id)).one()
            self.assertEqual(initial_rows, final_rows)

@unittest.skip(reason="not routinely necessary ")
class Test_count_per_day_6(Test_PCA_Database):
    """tests bulk upload with large numbers of samples """

    def runTest(self):
        for engine in self.engines.keys():
            print(engine, "#6")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            initial_rows = 1000000
            upload_data = []

            for i in range(initial_rows):
                upload_data.append({"bulk1": i, "bulk2": i * 1000})

            upload_df = pd.DataFrame.from_records(upload_data)
            pdm._bulk_load(upload_df, "test")

            (final_rows,) = pdm.session.query(func.count(BulkLoadTest.blt_int_id)).one()
            self.assertEqual(initial_rows, final_rows)

@unittest.skip(reason="not routinely necessary ")
class Test_count_per_day_7(Test_PCA_Database):
    """tests computation of count data frames"""

    def runTest(self):

        # load a variation model for testiong
        inputfile = "testdata/pca/vm.pickle"
        with open(inputfile, "rb") as f:
            vm = pickle.load(f)

            for engine in self.engines.keys():
                print(engine, "#7")

                pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
                pdm.store_variation_model(vm)
                pdm.store_cog_metadata(cogfile="testdata/pca/cog_metadata_vm.csv")
                pdm.make_contingency_tables()

class Test_dumpneighbours_8(Test_PCA_Database):
    """ tests bulk loading of fn4 snp data into an RDBMS.
        Specifically, it tests the output of findNeighbour4_dumpneighbours.py, which is design for one-off dumping and loading into an empty database. """

    def runTest(self):

        # load a variation model for testiong
        dumpdir = 'testdata/fn4snp/dumpneighbours'
        for engine in self.engines.keys():
            print(engine, "#7")

            pdm = PCADatabaseManager(connection_config=self.engines[engine], debug=True)
            pdm.fn4_bulk_upload(dumpdir)
              

        # check all links loaded
        loaded_links, = self.session.query(func.count(Neighbour.neighbour_int_id)).one_or_none()
        self.assertEqual(loaded_links, 25)
        loaded_samples, = self.session.query(func.count(FN4Sample.fn4s_int_id)).one_or_none()
        self.assertEqual(loaded_samples, 44)
        

        

