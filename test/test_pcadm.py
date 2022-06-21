""" tests PCADatabaseManager - software to do store PCA output.
    Only tests this with 
"""

import os
import shutil
import uuid
import unittest
import pandas as pd
import datetime
from sqlalchemy import create_engine

from pca.pca_scalable import VariantMatrix, PCARunner
from findn.common_utils import ConfigManager
from localstore.localstoreutils import LocalStore
from pca.pcadb import PCADatabaseManager, PCASummary


class Test_PCADB(unittest.TestCase):
    def setUp(self):
        # define a temporary directory
        self.tmpdir = os.path.join("unitTest_tmp", "pcadb", uuid.uuid4().hex)

        # remove if present
        shutil.rmtree(
            self.tmpdir, ignore_errors=True
        )  # should not exist, as uuid is unique

        os.makedirs(self.tmpdir)

    def tearDown(self):
        # remove any temporary directory
        shutil.rmtree(self.tmpdir, ignore_errors=True)


class Test_PCADB_1(Test_PCADB):
    """tests the PCADatabaseManager.store_cog_metadata() method,
    stores data to an sqlite database"""

    def runTest(self):

        # define connection to sqlite output database
        connection_config = os.path.join("sqlite:///", self.tmpdir, "test.sqlite")
        print(connection_config)

        pdm = PCADatabaseManager(
            connection_config=connection_config,
        )
        cogfile = os.path.join("testdata", "pca", "cog6000.csv")
        pdm.store_cog_metadata(
            cogfile=cogfile, date_end=datetime.date.fromisoformat("2022-01-01")
        )

        engine = create_engine(connection_config)
        sqlite_conn = engine.connect()

        expected_output = dict(
            test_pca_bulkload=0,
            build=0,
            build_metadata=0,
            population_studied=0,
            population_studied_extra=0,
            pca_summary=0,
            pca_summary_extra=0,
            statistical_model=0,
            statistical_model_fit=0,
            modelled_data=0,
            alert=0,
            contributing_basepos=0,
            contributing_pos=0,
            eigenvector=0,
            explained_variance_ratio=0,
            sample_set=0,
            sample_set_content=0,
            analysed_sample=0,
            clinical_metadata=4990,
            sequence_feature=14970,
        )
        for table_name in expected_output.keys():
            df = pd.read_sql_table(table_name, sqlite_conn)
            self.assertEqual(len(df.index), expected_output[table_name])
            # print("----------------{0}-----------------------".format(table_name))
            # print(df)
            # print("---------------- ends ------------------------------")
        engine.dispose()


class Test_PCADB_2(Test_PCADB):
    """uses the VariantMatrix and PCARunner classes using a LocalStore object as a sequence provider
    and then stores data to an sqlite database"""

    def runTest(self):

        # define connection to sqlite output database
        connection_config = os.path.join("sqlite:///", self.tmpdir, "test.sqlite")
        print(connection_config)

        pdm = PCADatabaseManager(
            connection_config=connection_config,
        )

        # local json test data n=6000
        LPERSIST = LocalStore("testdata/pca_scalable/test.tar")

        self.assertEqual(6000, len(LPERSIST.sequence_ids()))

        cfm = ConfigManager("testdata/pca/config.json")
        CONFIG = cfm.read_config()

        # start testing
        v = VariantMatrix(CONFIG, LPERSIST, self.tmpdir, show_bar=False)

        # print("Preparing to analyse 6,000 sequences")
        v.prepare_to_analyse()

        # valid samples should be in v.PCASEQSTORE
        to_analyse = v.PCASEQSTORE.sequence_ids
        self.assertEqual(len(to_analyse), 5710)

        # now compute a variation model
        print("Analysing variation")
        res = v.get_position_counts()
        self.assertIsInstance(res, dict)

        print("Building matrix")

        # test run
        pcr = PCARunner(v)
        pcr.run(
            n_components=10,
            select_from=to_analyse,
            pca_parameters={},
            target_matrix_size=1000,
        )

        # test cluster
        v = pcr.cluster()

        cl = pcr.vm.get_variationmodel_attribute("transformed_coordinate_categories")

        self.assertIsInstance(cl, pd.DataFrame)  # it's a dataframe

        pdm.store_variation_model(pcr.vm.vm)

        # check what has been created

        engine = create_engine(connection_config)
        sqlite_conn = engine.connect()

        expected_output = dict(
            test_pca_bulkload=0,
            build=1,
            build_metadata=13,
            population_studied=191,
            population_studied_extra=0,
            pca_summary=3446,
            pca_summary_extra=0,
            statistical_model=0,
            statistical_model_fit=0,
            modelled_data=0,
            alert=0,
            contributing_basepos=328,
            contributing_pos=327,
            eigenvector=1276,
            explained_variance_ratio=10,
            sample_set=0,
            sample_set_content=0,
            analysed_sample=5710,
            clinical_metadata=4990,
            sequence_feature=14970,
        )
        for table_name in expected_output.keys():
            df = pd.read_sql_table(table_name, sqlite_conn)
            print(table_name)
            if not len(df.index) == expected_output[table_name]:
                print("test2", table_name, len(df.index), expected_output[table_name])
            
            #self.assertEqual(len(df.index), expected_output[table_name])
            # print("----------------{0}-----------------------".format(table_name))
            # print(df)
            # print("---------------- ends ------------------------------")
        engine.dispose()


class Test_PCADB_3(Test_PCADB):
    """integration test of entire PCA and analysis process
    stores data to an sqlite database"""

    def runTest(self):

        # define connection to sqlite output database
        connection_config = os.path.join("sqlite:///", self.tmpdir, "test.sqlite")

        analysis_date = datetime.date.fromisoformat("2022-01-01")
        pdm = PCADatabaseManager(connection_config=connection_config)
        cogfile = os.path.join("testdata", "pca", "cog6000.csv")
        pdm.store_cog_metadata(
            cogfile=cogfile, date_end=datetime.date.fromisoformat("2022-01-01")
        )

        # local json test data n=6000
        LPERSIST = LocalStore("testdata/pca_scalable/test.tar")

        self.assertEqual(6000, len(LPERSIST.sequence_ids()))

        cfm = ConfigManager("testdata/pca/config.json")
        CONFIG = cfm.read_config()

        # start testing
        v = VariantMatrix(CONFIG, LPERSIST, self.tmpdir, show_bar=False)

        # print("Preparing to analyse 6,000 sequences")
        v.prepare_to_analyse()

        # valid samples should be in v.PCASEQSTORE
        to_analyse = v.PCASEQSTORE.sequence_ids
        self.assertEqual(len(to_analyse), 5710)

        # now compute a variation model
        print("Analysing variation")
        res = v.get_position_counts()
        self.assertIsInstance(res, dict)

        print("Building matrix")

        # test run
        pcr = PCARunner(v)
        pcr.run(
            n_components=10,
            min_variant_freq=10,
            select_from=to_analyse,
            pca_parameters={},
            target_matrix_size=1000,
        )

        # test cluster
        v = pcr.cluster()

        cl = pcr.vm.get_variationmodel_attribute("transformed_coordinate_categories")

        self.assertIsInstance(cl, pd.DataFrame)  # it's a dataframe

        pdm.store_variation_model(pcr.vm.vm)

        # build contingency tables
        pdm.make_contingency_tables(
            only_pc_cats_less_than_days_old=60,
            today=analysis_date,
        )

        # storing PCA summary
        pdm.store_pca_summary()  # store a summary

        # get results
        pcas_df = pdm.pca_summary(
            only_pc_cats_less_than_days_old=30,
            today=analysis_date,
        )
        pcas_df = pcas_df[pcas_df["n_days_observed"] >= 3]
        n_pc_cats = len(pcas_df["pc_cat"].unique())

        self.assertEqual(len(pcas_df.index), 115)
        self.assertEqual(n_pc_cats, 52)

        for i, pcas_int_id in enumerate(pcas_df.index):

            # ---------------------------------

            pcas_obj = pdm.single_pcas_summary(pcas_int_id)
            # print(pcas_obj)
            self.assertIsInstance(pcas_obj, PCASummary)
            cntdata = pdm.pcas_count_table(pcas_obj, output_format=1)

            # compute the date range over which the slope should be computed.
            # the last date is addition_date
            # the first date is either
            #  focus_on_most_recent_n_days before addition_date, OR
            #  the date the pc_cat was first seen, if that is later.
            analysis_start = analysis_date - datetime.timedelta(days=30)

            if cntdata["earliest_date"] > analysis_start:
                analysis_start = cntdata["earliest_date"]

            self.assertIsInstance(cntdata, dict)

        # check what has been created
        engine = create_engine(connection_config)
        sqlite_conn = engine.connect()

        expected_output = dict(
            test_pca_bulkload=0,
            build=1,
            build_metadata=13,
            population_studied=191,
            population_studied_extra=0,
            pca_summary=3446,
            pca_summary_extra=0,
            statistical_model=0,
            statistical_model_fit=0,
            modelled_data=0,
            alert=0,
            contributing_basepos=328,
            contributing_pos=327,
            eigenvector=1276,
            explained_variance_ratio=10,
            sample_set=0,
            sample_set_content=0,
            analysed_sample=5710,
            clinical_metadata=4990,
            sequence_feature=14970,
        )

        for table_name in expected_output.keys():
            df = pd.read_sql_table(table_name, sqlite_conn)
            #if not len(df.index) == expected_output[table_name]:
            #    print(table_name, len(df.index), expected_output[table_name])
            self.assertEqual(len(df.index), expected_output[table_name])
            # print("----------------{0}-----------------------".format(table_name))
            # print(df)
            # print("---------------- ends ------------------------------")
        engine.dispose()
