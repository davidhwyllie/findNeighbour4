""" tests pcadb.py - software to do store output from PCA

"""
import os
import json
import unittest
import pandas as pd
import pickle
import datetime
from pca.fittrend import PoissonModel
from pca.pcadb import PCADatabaseManager


class Test_pois_1(unittest.TestCase):
    """test Poisson 1"""

    def runTest(self):
        # load data for testing - the output from PCADatabaseManager.count_table()
        inputfile = "testdata/pca/cntdata.pickle"
        with open(inputfile, "rb") as f:
            cntdata = pickle.load(f)

        # test constructor and data munging processes
        this_latest_date = datetime.date(2021, 6, 1)
        nb = PoissonModel(**cntdata, latest_date=this_latest_date)
        self.assertEqual(max(nb.date_range["sample_date"]), this_latest_date)
        self.assertIsInstance(nb.counts, pd.DataFrame)
        self.assertIsInstance(nb.denominators, pd.DataFrame)

        # test  select_data_to_model methods
        res = nb.data_to_model(1)
        self.assertIsInstance(res, pd.DataFrame)
        self.assertEqual(len(res.index), 87)


class Test_pois_2(unittest.TestCase):
    """test Poisson model 2"""

    def runTest(self):
        # load data for testing - the output from PCADatabaseManager.count_table()
        inputfile = "testdata/pca/cntdata.pickle"
        with open(inputfile, "rb") as f:
            cntdata = pickle.load(f)

        # test constructor and data munging processes
        this_latest_date = datetime.date(2021, 6, 1)

        nb = PoissonModel(**cntdata, latest_date=this_latest_date)

        retVal = nb.fit()
        self.assertEqual(
            set(["statistical_model", "modelled_data", "coefficients"]),
            set(retVal.keys()),
        )
        self.assertEqual("completed", retVal["statistical_model"]["analysis_status"])
        self.assertIsNone(retVal["statistical_model"]["errors_returned"])
        self.assertEqual(
            set(retVal["modelled_data"].columns.to_list()),
            set(["n_total", "sample_date", "day_of_week", "n", "t", "pred"]),
        )

        # inputfile = "testdata/pca/model_success.pickle"
        # with open(inputfile, "wb") as f:
        #    pickle.dump(retVal, f)

        retVal = nb.fit(raise_error_for_unittest=True)
        self.assertEqual(
            set(["statistical_model", "modelled_data", "coefficients"]),
            set(retVal.keys()),
        )
        self.assertEqual("failed", retVal["statistical_model"]["analysis_status"])
        self.assertIsNotNone(retVal["statistical_model"]["errors_returned"])

        self.assertEqual(
            set(retVal["modelled_data"].columns.to_list()),
            set(["n_total", "sample_date", "day_of_week", "n", "t", "pred"]),
        )
        # inputfile = "testdata/pca/model_failed.pickle"
        # with open(inputfile, "wb") as f:
        #    pickle.dump(retVal, f)


class Test_PCA_Database(unittest.TestCase):
    """establishes database connection strings for cross-database testing.
    Currently tests OCI (if relevant environment variables are set) and Sqlite

    To test other databases, such as MySql, add the relevant connection string &
    database name to the dictionary self.engines"""

    def setUp(self):
        self.engines = {}  # {'None_specified':None}

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
            return
        if not os.path.exists(conn_detail_file):
            raise FileNotFoundError(
                "Connection file specified but not found: {0}".format(conn_detail_file)
            )
        # try to read the config file
        with open(conn_detail_file, "rt") as f:
            conn_detail = json.load(f)
            for key in conn_detail.keys():
                if key.startswith("unittest_ora"):
                    self.engines[
                        key
                    ] = key  # if there is a configured oracle connection string, then test using this
                    pass

        self.engines["Sqlite"] = "sqlite://"  # in memory sqlite
        return  # DEBUG - just test oracle for now


class Test_create_pc_summary_13(Test_PCA_Database):
    """tests addition of statistical models"""

    def runTest(self):

        # load a variation model for testiong
        inputfile = "testdata/pca/vm.pickle"
        with open(inputfile, "rb") as f:
            vm = pickle.load(f)

        for engine in self.engines.keys():

            print(engine, "#13")
            pdm = PCADatabaseManager(
                connection_config=self.engines[engine], debug=True, show_bar=False
            )
            pdm.store_variation_model(vm)
            pdm.store_cog_metadata(cogfile="testdata/pca/cog_metadata_vm.csv")
            pdm.store_pca_summary()  # store a summary

            pcas_df = pdm.pca_summary()

            n_modelled = 0
            for pcas_int_id in pcas_df.index:
                n_modelled += 1

                if n_modelled == 20:
                    break

                pcas_obj = pdm.single_pcas_summary(pcas_int_id)
                cntdata = pdm.count_table(pcas_obj)

                # override earliest date as required
                this_latest_date = datetime.date(2021, 6, 1)
                nb = PoissonModel(**cntdata, latest_date=this_latest_date)
                res = nb.fit()

                pdm.store_pcas_model_output(res)
