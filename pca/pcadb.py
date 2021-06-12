#!/usr/bin/env python
""" 
A component of a findNeighbour4 server which provides relatedness information for bacterial genomes.
Provides a persistence layer for the output of PCA, including modelling of category frequencies over time.

Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without tcen the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
"""

import os
import json
import pandas as pd
import logging
import numpy as np
import datetime
import warnings
import time
import glob
import progressbar
import statsmodels.api as sm
from sqlalchemy import (
    Integer,
    Boolean,
    Column,
    Float,
    Index,
    MetaData,
    Table,
    Text,
    String,
    DateTime,
    Date,
    Identity,
    ForeignKey
)
from sqlalchemy import create_engine, inspect, func
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declarative_base

import cx_Oracle
from pca.pca import VariationModel 


# global: definition of database structure
# classes mapping to persistence database inherit from this

db_pc = declarative_base()  # global

metadata = MetaData()


class PCADBManagerError(Exception):
    "a general purpose error used by the pcadb module."
    pass


class FN4Sample(db_pc):
    """ samples present in fn4 """

    __tablename__ = "fn4_sample"
    fn4s_int_id = Column(Integer, Identity(start=1), primary_key=True)
    sample_id = Column(String(38), index=True, unique=True)
    is_invalid = Column(Boolean)
    neighbours = relationship("Neighbour", backref="fn4sample")

class Neighbour(db_pc):
    """edges of samples present in fn4 """

    __tablename__ = "neighbour"
    neighbour_int_id = Column(Integer, Identity(start=1), primary_key=True)
    fn4s_int_id_1 = Column(Integer, ForeignKey(FN4Sample.fn4s_int_id), index=True)
    fn4s_int_id_2 = Column(Integer)
    dist = Column(Integer)
    
class BulkLoadTest(db_pc):
    """used for testing bulk uploads"""

    __tablename__ = "test"
    blt_int_id = Column(Integer, Identity(start=1), primary_key=True)
    bulk1 = Column(Integer)
    bulk2 = Column(Integer)


class Build(db_pc):
    """lists PCA builds"""

    __tablename__ = "build"
    build_int_id = Column(Integer, Identity(start=1), primary_key=True)
    builder = Column(String(48))
    build_time = Column(DateTime, index=True)
    model_load_start = Column(DateTime, nullable=True)
    model_load_complete = Column(DateTime, nullable=True)
    model_loaded = Column(Boolean)
    build_annotations = relationship("BuildAnnotation", backref="Build")
    contributing_basepos = relationship("ContributingBasePos", backref="Build")
    contributing_pos = relationship("ContributingPos", backref="Build")
    eigenvector = relationship("EigenVector", backref="Build")
    explained_variance_ratio = relationship("ExplainedVarianceRatio", backref="Build")
    sample = relationship("Sample", backref="Build")
    featassoc = relationship("FeatureAssociation", backref="build")
    statmodel = relationship("StatisticalModel", backref="build")

class BuildAnnotation(db_pc):
    """Metadata related to a PCA build
    example data:
    0|built|True|<class 'bool'>
    1|num_train_on|498669|<class 'int'>
    """

    __tablename__ = "build_metadata"
    bm_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    variable = Column(String(36))
    value = Column("value", String(1024))
    native_type = Column("native_type", String(16))


class ContributingBasePos(db_pc):
    """base positions contributing to the model"""

    __tablename__ = "contributing_basepos"
    cbp_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    basepos = Column(String(12))


class ContributingPos(db_pc):
    """positions contributing to the model"""

    __tablename__ = "contributing_pos"
    cp_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    pos = Column(Integer)


class EigenVector(db_pc):
    """weights at each positions in order to compute transformed coordinates"""

    __tablename__ = "eigenvector"
    ev_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    pc = Column(Integer)
    pos = Column(Integer)
    allele = Column(String(8))
    col = Column(String(8))
    weight = Column(Float)
    outside_3mad = Column(Boolean)


class ExplainedVarianceRatio(db_pc):
    """positions contributing to the model"""

    __tablename__ = "explained_variance_ratio"
    evr_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    pc = Column(Integer)
    explained_variance_ratio = Column(Float)
    pos_per_pc = Column(Integer)


class Sample(db_pc):
    """Samples considered for PCA & quality info"""

    __tablename__ = "analysed_sample"
    sample_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    sample_id = Column(String(38), index=True)
    m_in_model = Column(Integer)
    n_in_model = Column(Integer)
    model_positions = Column(Integer)
    reference_positions = Column(Integer)

    m_total = Column(Integer)
    m_expected_proportion = Column(Float)
    m_observed_proportion = Column(Float)
    m_p_value = Column(Float)

    n_total = Column(Integer)
    n_expected_proportion = Column(Float)
    n_observed_proportion = Column(Float)
    n_p_value = Column(Float)

    used_in_pca = Column(Boolean, index=True)
    suspect_quality = Column(Boolean)
    tcc = relationship("TransformedCoordinateCategory", backref="Build")


Index(
    "ix_SAMPLE",
    Sample.build_int_id,
    Sample.sample_id,
    Sample.used_in_pca
)


class TransformedCoordinateCategory(db_pc):
    """results of the pca"""

    __tablename__ = "transformed_coordinate_category"
    tcc_int_id = Column(Integer, Identity(start=1), primary_key=True)
    sample_int_id = Column(Integer, ForeignKey(Sample.sample_int_id), index=True)
    transformed_coordinate = Column(Float)
    pc = Column(Integer)
    cat = Column(Integer)
    pc_cat = Column(String(8))
    statmodel = relationship("StatisticalModel", backref="TCC")

Index(
    "ix_TCC",
    TransformedCoordinateCategory.sample_int_id,
    TransformedCoordinateCategory.pc_cat,
)

class StatisticalModel(db_pc):
    """types of statistical model applied to the data"""

    __tablename__ = "statistical_model"
    statmodel_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    tcc_int_id = Column(
        Integer, ForeignKey(TransformedCoordinateCategory.tcc_int_id), nullable=True, index=True
    )
    analysis_family_id = Column(String(38))
    analysis_id = Column(String(38))
    analysis_readable_title = Column(String(30))
    analysis_type = Column(String(30))
    analysis_software = Column(String(30))
    analysis_status = Column(String(8), index=True)
    readable_description = Column(Text)
    input_data = Column(Text, nullable= True)
    run_parameters = Column(Text, nullable=True)
    output_data = Column(Text, nullable=True)
    pc = Column(Integer, nullable = True)
    pc_cat = Column(String(8), nullable = True)
    date_start = Column(DateTime, nullable = True)
    date_end = Column(DateTime, nullable = True)
    interval_analysed_days = Column(Integer, nullable=True)
    fit = relationship("StatisticalModelFit", backref="StatModel")
    alert = relationship("Alert", backref="StatModel")

class Alert(db_pc):
    """alerts based on statistical models"""

    __tablename__ = "alert"
    alert_int_id = Column(Integer, Identity(start=1), primary_key=True)
    statmodel_int_id = Column(
        Integer, ForeignKey(StatisticalModel.statmodel_int_id), index=True
    )
    alert_rule = Column(String(18))
    alert_description = Column(String(100))
    pc_cat = Column(String(8))

class StatisticalModelFit(db_pc):
    """statistical models fitted"""

    __tablename__ = "statistical_model_fit"
    smf_int_id = Column(Integer, Identity(start=1), primary_key=True)
    statmodel_int_id = Column(
        Integer, ForeignKey(StatisticalModel.statmodel_int_id), index=True
    )
    analysis_family_id = Column(String(38))
    analysis_id = Column(String(38))
    param = Column(String(16))
    estimate = Column(Float)
    std_error = Column(Float)
    estimate_ci_low = Column(Float)
    estimate_ci_high = Column(Float)
    p_value = Column(Float, index=True)
    pc_cat = Column(String(8))
    param_desc = Column(String(255))
    is_reference = Column(Boolean)
    comments = Column(Text)
    has_ci = Column(Boolean)
    estimate2naturalspace = Column(String(8))

class ClinicalMetadata(db_pc):
    """holds clinical metadata, if it exists.
    Note that there is no guarantee or assumption made about the order in which the clinical metadata will
    be made available relative to sequence data.
    For this reason, FK constraints are not enforced between ClinicalMetaData and the Sample table."""

    __tablename__ = "clinical_metadata"
    cm_int_id = Column(Integer, Identity(start=1), primary_key=True)
    sample_id = Column(String(38), index=True, unique=True)
    sample_date = Column(Date, nullable=True)
    data_addition_datetime = Column(DateTime, index=True)
    country = Column(String(50), nullable=True)
    adm1 = Column(String(50), nullable=True)
    adm2 = Column(String(50), nullable=True)
    residence_longitude = Column(Float, nullable=True)
    residence_latitude = Column(Float, nullable=True)
    surname = Column(String(50), nullable=True)
    forename = Column(String(50), nullable=True)
    birthdate = Column(Date, nullable=True)
    addr1 = Column(String(128), nullable=True)
    postcode = Column(String(14), nullable=True)
    id = relationship("ClinicalIdentifier", backref="ClinicalMetadata")
    seqfeature = relationship("SequenceFeature", backref="ClinicalMetadata")


class ClinicalIdentifier(db_pc):
    """holds clinical metadata, if it exists.
    Note that there is no guarantee or assumption made about the order in which the clinical metadata will
    be made available relative to sequence data.
    For this reason, FK constraints are not enforced between ClinicalMetaData and the Sample table."""

    __tablename__ = "clinical_identifier"
    ci_int_id = Column(Integer, Identity(start=1), primary_key=True)
    cm_int_id = Column(Integer, ForeignKey(ClinicalMetadata.cm_int_id))
    administration = Column(String(8))
    id_type = Column(String(16))
    id_value = Column(String(50))

Index(
    "ix_clinical_ids",
    ClinicalIdentifier.administration,
    ClinicalIdentifier.id_type,
    ClinicalIdentifier.id_value,
)

class SequenceFeature(db_pc):
    """holds identifiers related to sequences, such as SNPs, lineages, etc.
    Note that there is no guarantee or assumption made about the order in which the clinical metadata will
    be made available relative to sequence data.
    For this reason, FK constraints are not enforced between ClinicalMetaData and the Sample table."""

    __tablename__ = "sequence_feature"
    sf_int_id = Column(Integer, Identity(start=1), primary_key=True)
    cm_int_id = Column(Integer, ForeignKey(ClinicalMetadata.cm_int_id), index=True)
    variable = Column(String(24))
    value = Column(String(80))
    sequencefeature = Column(String(105), index=True)
    version = Column(String(38), nullable=True)
    evidence_for = Column(Float, nullable=True)
    evidence_against = Column(Float, nullable=True)
    is_lineage = Column(Boolean, nullable=False)


Index("ix_sequence_ids", SequenceFeature.variable, SequenceFeature.value)


class FeatureAssociation(db_pc):
    """associations between pc_cats and sequencefeatures

    contains a 2x2 contingency table, see below

    as well as an estimated odds ratio.
    """

    __tablename__ = "feature_association"
    featassoc_int_id = Column(Integer, Identity(start=1), primary_key=True)
    build_int_id = Column(Integer, ForeignKey(Build.build_int_id), index=True)
    pc_cat = Column(String(8), index=True)
    sequencefeature = Column(String(50), index=True)
    a = Column(Integer)
    b = Column(Integer)
    c = Column(Integer)
    d = Column(Integer)
    df = Column(Integer)
    p_value = Column(Float, index=True)
    log_or = Column(Float)
    log_or_lower_ci = Column(Float)
    log_or_upper_ci = Column(Float)
    sens = Column(Float, index=True)
    spec = Column(Float)
    ppv = Column(Float)
    npv = Column(Float)


class ContingencyTable:
    """computes statistics on a 2x2 contingency table"""

    def __init__(self, build_int_id, pc_cat, sequencefeature, a, b, c, d):
        """expects as parameters count data for a 2x2 contingency table

        Feature    Present       Absent
        PC-CAT Y      a            b          a+b
               N      c            d          c+d
                    a+c           b+d       a+b+c+d

        build_int_id, pc_cat, sequence_feature, tokens used to identify the results

        The 'Feature' is the 'True result' and we are interested in whether the pc_cat predicts the 'True result'
        a: True positives
        b: False positive
        c: False negative
        d: True negative
        reports an estimated odds ratio.

        returns:
        a FeatureAssociation object"""

        self.results = {
            "a": int(a),
            "b": int(b),
            "c": int(c),
            "d": int(d),
            "build_int_id": build_int_id,
            "pc_cat": pc_cat,
            "sequencefeature": sequencefeature,
        }

        # shift zeros up by 0.5 to prevent infinite results.  Haldane-Anscombe correction.  See also:
        # Ruxton and Neuhäuser 2013 (Review of alternative approaches to calculation of a confidence interval for the odds ratio of a 2x2 contingency table)
        if b == 0 or d == 0:
            a = a + 0.5
            b = b + 0.5
            c = c + 0.5
            d = d + 0.5
        sens = a / (a + c)
        spec = b / (b + d)
        ppv = a / (a + b)
        npv = d / (c + d)

        # construct matrix
        self.mat = np.array([[a, c], [b, d]])
        self.tbl = sm.stats.Table2x2(self.mat)  # added 0.5 to any empty cells
        stbl = self.tbl.summary()
        self.results["log_or"] = float(stbl[2][1].data)
        self.results["log_or_lower_ci"] = float(stbl[2][3].data)
        self.results["log_or_upper_ci"] = float(stbl[2][4].data)
        self.results["p_value"] = float(stbl[2][5].data)
        self.results["df"] = 1
        self.results["sens"] = sens
        self.results["spec"] = spec
        self.results["ppv"] = ppv
        self.results["npv"] = npv

    def summary(self):
        return self.tbl.summary()

    def featureassociation(self):
        return FeatureAssociation(**self.results)

    # DEBUG: probably need to remove this |  not used at present
    t_trend_summary = Table(
        "trend_summary",
        metadata,
        Column("pc_cat", Text, index=True),
        Column("n_samples", Integer),
        Column("min_trans_coord", Float),
        Column("max_trans_coord", Float),
        Column("mean_trans_coord", Float),
        Column("n_gt_14_days_before", Float),
        Column("n_gt_30_days_before", Float),
        Column("n_gt_60_days_before", Float),
        Column("n_gt_90_days_before", Float),
        Column("n_gt_120_days_before", Float),
        Column("earliest_date", Text),
        Column("latest_date", Text),
        Column("n_days_observed", Integer),
        Column("analysis_id", Text),
        Column("analysis", Text),
        Column("analysis_type", Text),
        Column("date_end", Text),
        Column("interval_analysed_days", Integer),
        Column("param", Text),
        Column("param_desc", Text),
        Column("Estimate", Float),
        Column("Estimate_CI_low", Float),
        Column("Estimate_CI_high", Float),
        Column("p_value", Float),
        Column("is_reference", Float),
        Column("has_CI", Float),
        Column("Estimate2NaturalSpace", Text),
    )
    

class PCADatabaseManager:
    """manages an RDBMS containing PCA output"""

    def __init__(self, connection_config=None, debug=False, show_bar = True):
        """creates the RDBMS connection

        Parameters
        -----------
        connection_config:  
        One of 
        1. a key to a dictionary containing one or more database configuration details: (e.g. 'prod', 'test')
        2. a valid sqlalchemy database connection string (if this is sufficient for connections)  e.g. 'pyodbc+mssql://myserver'
        3. None.  This is considered to mean 'sqlite://' i.e. an in memory sqlite database, which is not persisted when the program stops.
        
        if it is not none, a variable called PCA_CONNECTION_CONFIG_FILE must be present.  This must point to a file containing credentials.
        the name of an environment variable containing (in json format) a dictionary, or None if it is not required.
        An example of such a dictionary is as below:
        {
            'prod':{'DBTYPE':'sqlite', 'ENGINE_NAME':'sqlite:///db/proddb.sqlite'}
            'dev': {'DBTYPE':'sqlite', 'ENGINE_NAME':'sqlite:///db/devdb.sqlite'},
            'test':{'DBTYPE':'sqlite', 'ENGINE_NAME':'sqlite://'}
        }
        The DBTYPE and ENGINE_NAME keys are essential.
        Other keys may also be present, and are required in some cases (for example by Oracle connections).
        {
            'prod':{'DBTYPE':'oracle', 
                    'ENGINE_NAME':''oracle+cx_oracle://PROD:97bxxxxxxxxX@(description: .........)))',
                    'TNS_ADMIN':'/secrets/oracle/pca_prod'
                    },

            'dev':{'DBTYPE':'oracle', 
                    'ENGINE_NAME':''oracle+cx_oracle://PROD:97bxxxxxxxxX@(description: .........)))',
                    'TNS_ADMIN':'/secrets/oracle/pca_prod'
                    }
        }
        Note, the bit after the @(description describes where your database is, and will be found in your cloud wallet, if you are using cloud databases.  See below.
        In this case, TNS_ADMIN is the value you wish the TNS_ADMIN environment variable to be set to. 
        The software will set TNS_ADMIN, and, if you are using a virtual environment, it will be scoped to the virtual environment.
        In summary, it will configure necessary settings to allow Oracle database connections.
                    
        Note, if you are using a python virtual environment, this environment variable should be included in the .env file in the root of your project.  The .env file should not be under source control.

        configuration engine_name: an SQLalchemy connect string, e.g. sqlite::// for temporary memory db, see https://docs.sqlalchemy.org/en/13/core/engines.html
        debug: if True, deletes any existing data on startup.
        show_bar: show a progress bar during long operations


        NOTE:
        This software has been tested with 
        (i) Sqlite 3.3.2+ 
        (ii) Oracle Autonomous Database (cloud)
        https://blogs.oracle.com/oraclemagazine/getting-started-with-autonomous
        
        To connect to Oracle there are several steps.
        0. Install dependencies, see
        https://cx-oracle.readthedocs.io/en/latest/user_guide/installation.html
        wget https://download.oracle.com/otn_software/linux/instantclient/211000/instantclient-basic-linux.x64-21.1.0.0.0.zip
        need to set the LD_LIBRARY_PATH variable, see
        https://www.oracle.com/database/technologies/instant-client/linux-x86-64-downloads.html
        e.g. export LD_LIBRARY_PATH=/data/software/instantclient_21_1:LD_LIBRARY_PATH
        
        these parameters have to go into the python .env file e.g.
        ----------------------------------------------------------------
        LD_LIBRARY_PATH="/software/instantclient_21_1"
        PCADB_CONNECTION_CONFIGFILE="/secret/config.json"
        
        Where config.json looks like
        {
            'prod':{'DBTYPE':'oracle', 
                    'ENGINE_NAME':''oracle+cx_oracle://PROD:97bxxxxxxxxX@(description: .........)))',
                    'TNS_ADMIN':'/secrets/oracle/pca_prod'
                    }
        }
        ** NOTE: as per normal json conventions, escape quotes (i.e. \" not " around the certificate name, otherwise SSL connections will fail)  **

        1. Download your OCI wallet, & unzip it somewhere
        2. Set the TNS_ADMIN env var to point to this directory
        3. Edit the WALLET_LOCATION in the sqlnet.ora file to point to the relevant directory, e.g. WALLET_LOCATION = (SOURCE = (METHOD = file) (METHOD_DATA = (DIRECTORY="/data/credentials/oci_test")))
        4. Check the 
        4. Create a user with relevant privileges see below)
        5. Set the OCI_ENGINE_NAME env var.
        An example of this is as below (redacted so not live)
        oracle+cx_oracle://scott:tigerX22@(description= (retry_count=20)(retry_delay=3)(address=(protocol=tcps)(port=1522)(host=host.oraclecloud.com))(connect_data=(service_name=redacted))(security=(ssl_server_cert_dn="redacted")))


        This data can be found in the tnsnames.ora file: for details, see 
         https://docs.sqlalchemy.org/en/14/dialects/oracle.html#dialect-oracle-cx_oracle-connect
         https://stackoverflow.com/questions/37471892/using-sqlalchemy-dburi-with-oracle-using-external-password-store
         https://stackoverflow.com/questions/14140902/using-oracle-service-names-with-sqlalchemy/35215324
         https://blogs.oracle.com/sql/how-to-create-users-grant-them-privileges-and-remove-them-in-oracle-database
         https://docs.oracle.com/en/database/oracle/oracle-database/19/sqlrf/GRANT.html#GUID-20B4E2C0-A7F8-4BC8-A5E8-BE61BDC41AC3

        Configuring interactions with external OCI databases
        ====================================================
        Your application will need to run as a user (we'll call it PCADB) will need some priviledges granted.
        The exact privileges required involving creating, dropping tables & indexes, as well as inserting and deleting data.
        CREATE USER PCADB IDENTIFIED BY 'MyPassword1234!';         
        GRANT CONNECT TO PCADB;
        GRANT CREATE SESSION TO PCADB;
        GRANT CREATE SEQUENCE TO PCADB;
        GRANT CREATE TABLE TO PCADB;
        GRANT CREATE SYNONYM TO PCADB;
        ALTER USER PCADB DEFAULT TABLESPACE DATA quota unlimited on DATA; 
        
        """

        # connect and create session.  Validate inputs carefully.
        if connection_config is None:       
            logging.info("Connection config is None: using in-memory sqlite.")
            self.engine_name = 'sqlite://'
        elif "://" in connection_config:
            logging.info("Connection config provided; using {0}".format(connection_config))
            self.engine_name = connection_config
        else:
            # PCA_CONNECTION_CONFIG_FILE should contain credentials
            
            try:
                conn_detail_file = os.environ['PCA_CONNECTION_CONFIG_FILE']
            except KeyError:
                raise PCADBManagerError("Environment variable PCA_CONNECTION_CONFIG_FILE does not exist; however, it is required.  If you are using a python virtual environment, you need to set it in .env, not globally") 
            
            if not os.path.exists(conn_detail_file):
                raise FileNotFoundError("Connection file specified but not found: {0}".format(conn_detail_file))

            # read the config file
            with open(conn_detail_file,'rt') as f:
                conn_detail = json.load(f)

            if not connection_config in conn_detail.keys():
                raise PCADBManagerError("Connection {0} does not correspond to one of the keys {1} of the configuration json file at {2}".format(connection_config, conn_detail.keys(), conn_detail_file))

            # configure engine
            this_configuration = conn_detail[connection_config]         # extract the relevant part of the config dictionary
            
            # two keys are always present
            essential_keys = set(['DBTYPE','ENGINE_NAME'])
            if len(essential_keys - set(this_configuration.keys()))>0:
                raise PCADBManagerError("Provided keys for {0} are not correct.  Required are {1}".format(connection_config,essential_keys))

            # if it's Oracle, then three keys are required.
            if this_configuration['DBTYPE'] == 'oracle':
                essential_keys = set(['DBTYPE','ENGINE_NAME','TNS_ADMIN'])
                if len(essential_keys - set(this_configuration.keys()))>0:
                    raise PCADBManagerError("Provided keys for oracle db in {0} are not correct.  Required are {1}".format(connection_config,essential_keys))

                # set the TNS_ADMIN variable.
                logging.info("Set TNS_ADMIN to value specified in config file {0}".format(this_configuration['TNS_ADMIN']))
                os.environ['TNS_ADMIN'] = this_configuration['TNS_ADMIN']
                
            logging.info("Set ENGINE_NAME configuration string from config file.")
            self.engine_name = this_configuration['ENGINE_NAME']

        ## DEBUG
        #print(os.environ)
        ## 

        # now we can start
        self.Base = db_pc
        logging.info("PCADatabaseManager: Connecting to database")
        self.engine = create_engine(self.engine_name)
        
        self.is_oracle = "oracle+cx" in self.engine_name
        self.is_sqlite = "sqlite://" in self.engine_name
        self.show_bar = show_bar

        logging.info(
            "PCADatabaseManager: Database connection made; there are {0} tables.  Oracle database = {1}".format(
                len(self._table_names()), self.is_oracle
            )
        )
        
        # drop existing tables if in debug mode
        if debug:
            self._drop_existing_tables()

        self.Base.metadata.create_all(self.engine)  # create the table(s)
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    def _table_names(self):
        """returns table names in the schema"""
        return inspect(self.engine).get_table_names()

    def _drop_existing_tables(self):
        """drops any existing tables"""

        # test - wait 10 seconds - for any existing transactions to cleartherwise, ORA-00054 can occur during oracle testing
        # happens only with multiple cpus
        if self.is_oracle:
            time.sleep(0)

        self.Base.metadata.create_all(self.engine)  # create the table(s)
        BulkLoadTest.__table__.drop(self.engine)
        FeatureAssociation.__table__.drop(self.engine)
        SequenceFeature.__table__.drop(self.engine)
        ClinicalIdentifier.__table__.drop(self.engine)
        ClinicalMetadata.__table__.drop(self.engine)
        Alert.__table__.drop(self.engine)
        StatisticalModelFit.__table__.drop(self.engine)
        StatisticalModel.__table__.drop(self.engine)
        TransformedCoordinateCategory.__table__.drop(self.engine)
        EigenVector.__table__.drop(self.engine)
        ContributingPos.__table__.drop(self.engine)
        ContributingBasePos.__table__.drop(self.engine)
        ExplainedVarianceRatio.__table__.drop(self.engine)
        Sample.__table__.drop(self.engine)
        BuildAnnotation.__table__.drop(self.engine)
        Build.__table__.drop(self.engine)
        Neighbour.__table__.drop(self.engine)
        FN4Sample.__table__.drop(self.engine)

        remaining = len(self._table_names())
        if remaining > 0:
            raise PCADBManagerError(
                "Failed to remove all tables in the database.  The following remain: {0}".format(
                    self._table_names()
                )
            )

        # wait 10 seconds - otherwise, ORA-00054 can occur during oracle testing
        if self.is_oracle:
            time.sleep(0)
        return


    def _bulk_load(self, upload_df, target_table, max_batch=100000):
        """bulk uploads pandas dataframes.
        If using an oracle database, uses Oracle's cx-oracle package to bulk upload data

        upload_df: a pandas dataframe.  Names **much match** the table target_table
        target_table: the name of the table to upload into
        (note: this is not the name of the class in the table definition, it's the name of the table in the SQL)

        Returns:
        number of items uploaded (integer)

        Background:
        - ORM is slow for large inserts
        - Oracle is not currently supported by pandas .to_sql method='multi', so a bespoke method is needed
        - Need custom code for bulk loading to Oracle.  Methods are provided in the cx_oracle package, see
        https://cx-oracle.readthedocs.io/en/latest/user_guide/batch_statement.html


        The maximum size that can be transmitted to the Oracle server is 2G.  If htis happens, a cx_Oracle.DatabaseError
        is raised with result DPI-1015.  Reduce the max_batch

        
        """

        # verify data : ensure that target_table is a table [Essential: otherwise, can be a vector for SQL injection]
        if target_table not in self._table_names():
            raise ValueError(
                "target_table {0} does not exist: valid tables are {1}".format(
                    target_table, self._table_names()
                )
            )

        # check that upload_df is pandas dataframe
        if not isinstance(upload_df, pd.DataFrame):
            raise TypeError("upload_df needs to be pandas DataFrame, not a {0}".format(type(upload_df)))

        # check that we have been passed data
        ncol = len(upload_df.columns)
        if ncol == 0:
            # we treat this as an error
            raise PCADBManagerError("Passed data frame to _bulk_upload to {0} contains no columns {1}".format(target_table, upload_df)) 

        if self.is_sqlite:
            # there is a max variable limit of 32,766 for Sqlite 3.32.0 on https://www.sqlite.org/limits.html
            # set max_batch to keep below this.
            max_batch = int(32000 / ncol)
            logging.info("Autoset max_batch to {0}, as running SQLite".format(max_batch))

        if self.is_oracle:

            # parse the engine_name into dsn, database & password
            e1 = self.engine_name.replace("oracle+cx_oracle://", "")
            up, dsn = e1.split("@")
            u, p = up.split(":")

            # get into the right format for loading: note: holds all data in ram
            loadvar = list(upload_df.itertuples(index=False, name=None))
            ncol = len(upload_df.columns.to_list())

            # estimate the maximum buffer size required and auto-reduce the max_batch if required.
            estimated_row_size = len(str(loadvar[0]))
            estimated_max_batch = int(1e7/(estimated_row_size))       # large margin of safety 10M batch size
            if estimated_max_batch < max_batch:
                max_batch = estimated_max_batch
                logging.info("Reduced max_batch to keep estimated buffer size within acceptable limited (aim to control to 10M).  max_batch is {0}".format(max_batch))

            ## disabled: does not work.  The heuristic for estimating max batch doesn't work properly - DPI-1015 errors occur
            #else: 
            #    print("Larger max_batch may be possible while keeping buffer size required < 2G.  max_batch set to {0}".format(estimated_max_batch))
            #    max_batch = estimated_max_batch
            
            # construct sql statment.
            # Should be injection-safe; we have checked the target_table is a table, and are incorporating integers only.
            collabels = [":" + str(x + 1) for x in list(range(ncol))]
            insert_cols = ",".join(collabels)
            target_cols = ",".join(upload_df.columns.to_list())
            sql_statement = "INSERT INTO {0} ({1}) VALUES ({2})".format(
                target_table, target_cols, insert_cols
            )
            start_n = len(loadvar)
            if self.show_bar:
                bar = progressbar.ProgressBar(max_value=start_n)

            with cx_Oracle.connect(user=u, password=p, dsn=dsn) as conn:
                cursor = conn.cursor()

                while len(loadvar) > 0:
                    if self.show_bar:
                        bar.update(start_n - len(loadvar))
                    

                    cursor.executemany(sql_statement, loadvar[0:max_batch])
                    loadvar = loadvar[max_batch:]

                conn.commit()
            if self.show_bar:
                bar.finish()
        else:

            # note: there may be limits to the complexity of the statement permitted: a maximum number of parameters.
            # therefore, we do this in batches as well.
            start_n = len(upload_df.index)
            if self.show_bar:
                bar = progressbar.ProgressBar(max_value=start_n)

            while len(upload_df.index) > 0:
                logging.info(
                    "Bulk upload of {0} : {1} remain".format(
                        target_table, len(upload_df)
                    )
                )
                to_upload = upload_df.head(n=max_batch)
                to_upload.to_sql(
                    target_table,
                    self.engine,
                    if_exists="append",
                    index=False,
                    method="multi",
                )
                upload_df = upload_df.iloc[max_batch:]
                if self.show_bar:
                    bar.update(start_n - len(upload_df.index))

        logging.info("Bulk upload complete")
        if self.show_bar:
            bar.finish()
        return len(upload_df.index)

    def store_variation_model(self, vm):
        """stores a VariationModel object's contents to database"""

        logging.info("Preparing variation model for storage")

        if not isinstance(vm, VariationModel):
            raise TypeError(
                "Need to be passed a VariationModel, not a {0}".format(
                    type(VariationModel)
                )
            )

        logging.info("Creating new build.")
        this_build = Build(
            builder="pcadb.PCADatabaseManager.store_variation_model",
            build_time= datetime.datetime.fromisoformat(vm.model['build_time']),
            model_load_start = datetime.datetime.now(),
            model_loaded = 0
        )
        self.session.add(this_build)  # add the build

        metadata = []
        for key in vm.model:
            native_type = type(vm.model[key])
            if native_type in [bool, int, float, str]:
                metadata.append(
                    {
                        "variable": key,
                        "value": str(vm.model[key]),
                        "native_type": str(native_type),
                    }
                )
            elif type(vm.model[key]) == np.int64:
                metadata.append({"variable": key, "value": str(int(vm.model[key]))})

        # add annotations
        logging.info("Loading annotations")
        for item in metadata:
            this_build.build_annotations.append(BuildAnnotation(**item))
        # add explained variance ratio and positions per pc
        logging.info("Loading explained variance ratio")
        for pc, evr in enumerate(vm.model["explained_variance_ratio"]):
            this_build.explained_variance_ratio.append(
                ExplainedVarianceRatio(
                    pc=pc,
                    explained_variance_ratio=evr,
                    pos_per_pc=vm.model["pos_per_pc"][pc],
                )
            )

        # add base positions
        logging.info("Loading base positions")
        for item in vm.model["contributing_basepos"]:
            this_build.contributing_basepos.append(ContributingBasePos(basepos=item))
        # add positions
        logging.info("Loading contributing positions")
        for item in vm.model["contributing_pos"]:
            this_build.contributing_pos.append(ContributingPos(pos=item))

        logging.info("Committing initial model data")
        self.session.commit()

        # bulk load of eigenvectors
        vm.model["eigenvectors"]["build_int_id"] = this_build.build_int_id
        self._bulk_load(vm.model["eigenvectors"], "eigenvector")

        # load transformed coordinates
        logging.info("Loading transformed_coordinate_categories")
        tcc_df = vm.model["transformed_coordinate_categories"]

        # load sample mixture data
        sample_df = vm.model["mix_quality_info"]
        sample_df["build_int_id"] = this_build.build_int_id
        sample_df["sample_id"] = vm.model["mix_quality_info"].index

        # -- some databases (Oracle) won't store very small numbers, code them as zero
        # -- oracle yields DPI-1044: value cannot be represented as an Oracle number
        small_N_p = sample_df.index[sample_df["n_p_value"] < 1e-30]
        sample_df.loc[small_N_p, "n_p_value"] = 0
        small_M_p = sample_df.index[sample_df["m_p_value"] < 1e-30]
        sample_df.loc[small_M_p, "m_p_value"] = 0

        sample_df["used_in_pca"] = sample_df.index.isin(vm.model["sample_id"])
        sample_df["suspect_quality"] = sample_df.index.isin(
            vm.model["suspect_quality_seqs"].index
        )

        self._bulk_load(sample_df, "analysed_sample")

        # query sample, recovering the sample_int_ids for the inserted data.
        # this has to be merged into the TransformedCoordinateCategory data.
        get_sample_int_ids_sql = (
            self.session.query(Sample.sample_id, Sample.sample_int_id)
            .filter(Sample.build_int_id == this_build.build_int_id)
            .statement
        )
        sample_id_df = pd.read_sql(get_sample_int_ids_sql, con=self.engine)
        tcc_df = tcc_df.merge(sample_id_df, on="sample_id", how="inner")
        tcc_df.drop(["sample_id", "initial_cat"], axis=1, inplace=True)
        self._bulk_load(tcc_df, "transformed_coordinate_category")

        logging.info("Marking build as complete")
        this_build.model_loaded = 1
        this_build.model_load_complete = datetime.datetime.now()
        self.session.commit()
        logging.info("Build is complete")

    def _wipe_clinical_and_sequence_metadata(self):
        """deleted all data from the following:
        ClinicalMetadata, ClinicalIdentifier, and SequenceFeature
        """
        self.session.query(ClinicalIdentifier).delete()
        self.session.query(SequenceFeature).delete()
        self.session.query(ClinicalMetadata).delete()
        self.session.commit()

    def existing_sample_ids_in_clinical_metadata(self):
        """returns a set containing sample_ids already present in the clinical metadata tables"""
        existing_sample_ids = [
            x for x, in self.session.query(ClinicalMetadata.sample_id)
        ]
        return set(existing_sample_ids)

    def store_cog_metadata(
        self, cogfile, date_end=datetime.datetime.now(), replace_all=False
    ):
        """extracts data from cog-uk format metadata files and imports them into RDBMS

        cogfile: the cog-uk metadata file.
        date_end: a date or datetime value.  don't add information with specimen dates after this date
        replace_all: if True, then deletes all stored metadata and replaces.  If False, will only add new files.

        uses bulk uploads

        returns:
        number of samples added
        """

        # check the file exists
        if not os.path.exists(cogfile):
            raise FileNotFoundError("COG metadata file not found: {0}".format(cogfile))

        # read the cog metadata file
        logging.info("Loading data from file")
        cogdf = pd.read_csv(cogfile, header=0)

        # extract the sequence_id from the sequence_name
        cogdf["sample_id"] = [x[1] for x in cogdf["sequence_name"].str.split("/")]
        
        # wipe data if replace_all:
        if replace_all:
            self._wipe_clinical_and_sequence_metadata()

        # determine existing sample_ids stored in the database.
        logging.info("Checking existing data")
        existing_sample_ids = self.existing_sample_ids_in_clinical_metadata()

        # drop existing data, and duplicates.
        logging.info("Clinical data found in database. Rows = {0}".format(len(existing_sample_ids)))
        logging.info("COG_UK data loaded. Rows = {0}".format(len(cogdf.index)))
        cogdf.drop_duplicates(['sample_id'], inplace=True)
        logging.info("After deduplicating, Rows = {0}".format(len(cogdf.index)))
        
        drop_ix = cogdf[cogdf['sample_id'].isin(existing_sample_ids)].index
        logging.info("There are {0} rows which have already been loaded".format(len(drop_ix)))
        cogdf.drop(drop_ix, inplace=True)
        logging.info("There are {0} rows remaining to load".format(len(cogdf.index)))

        # load rest
        n_added = 0
        cm_to_insert = list()
        sm_to_insert = list()
        for ix in cogdf.index:
            if (
                datetime.datetime.fromisoformat(cogdf.at[ix, "sample_date"])
                <= date_end
            ):
                n_added += 1
                if n_added % 50000 == 0:
                    logging.info("Parsing cog-uk data file n={0}. ".format(n_added))

                cm = dict(
                    sample_id=cogdf.at[ix, "sample_id"],
                    sample_date=datetime.datetime.fromisoformat(
                        cogdf.at[ix, "sample_date"]
                    ),
                    data_addition_datetime=datetime.datetime.now(),
                    country=cogdf.at[ix, "country"],
                    adm1=cogdf.at[ix, "adm1"],
                )
                cm_to_insert.append(cm)
                # add the original sequence_name as a sequence identifier.
                sm_to_insert.append(
                    dict(
                        sample_id=cogdf.at[ix, "sample_id"],
                        variable="COG_UK_sequence_name",
                        value=cogdf.at[ix, "sequence_name"],
                        is_lineage=False,
                        evidence_for=1,
                        evidence_against=0,
                        version="-",
                        sequencefeature="",
                    )
                )

                # add the lineage
                lm = dict(
                    sample_id=cogdf.at[ix, "sample_id"],
                    variable="lineage",
                    value=cogdf.at[ix, "lineage"],
                    version=cogdf.at[ix, "lineages_version"],
                    sequencefeature="{0}:{1}".format(
                        "pangolearn", cogdf.at[ix, "lineage"]
                    ),
                    is_lineage=True,
                    evidence_for=1,
                    evidence_against=0,
                )

                # add scores if they are provided
                if not np.isnan(cogdf.at[ix, "lineage_ambiguity_score"]):
                    lm["evidence_for"] = cogdf.at[ix, "lineage_ambiguity_score"]
                if not np.isnan(cogdf.at[ix, "lineage_conflict"]):
                    lm["evidence_against"] = cogdf.at[ix, "lineage_conflict"]
                # provided the lineage isn't blank we'll add the record
                if not cogdf.at[ix, "lineage"] == "nan":
                    sm_to_insert.append(lm)

                # add the mutations
                mutations = [
                    "e484k",
                    "t1001i",
                    "d614g",
                    "p323l",
                    "del_21765_6",
                    "a222v",
                    "p681h",
                    "q27stop",
                    "n501y",
                    "n439k",
                    "y453f",
                    "del_1605_3",
                ]
                for mutation in mutations:
                    sm_to_insert.append(
                        dict(
                            sample_id=cogdf.at[ix, "sample_id"],
                            variable=mutation,
                            value=cogdf.at[ix, mutation],
                            sequencefeature="{0}:{1}".format(
                                mutation, cogdf.at[ix, mutation]
                            ),
                            is_lineage=False,
                            evidence_for=1,
                            evidence_against=0,
                            version="-",
                        )
                    )

                # scorpio classification
                sm = dict(
                    sample_id=cogdf.at[ix, "sample_id"],
                    variable="scorpio",
                    value=cogdf.at[ix, "scorpio_call"],
                    sequencefeature="{0}:{1}".format(
                        "scorpio", cogdf.at[ix, "scorpio_call"]
                    ),
                    is_lineage=True,
                    evidence_for=1,
                    evidence_against=0,
                )
                # add scores if they are provided
                if not np.isnan(cogdf.at[ix, "scorpio_support"]):
                    sm["evidence_for"] = cogdf.at[ix, "scorpio_support"]
                if not np.isnan(cogdf.at[ix, "scorpio_conflict"]):
                    sm["evidence_against"] = cogdf.at[ix, "scorpio_conflict"]
                # provided the lineage isn't blank we'll add the record
                if not cogdf.at[ix, "scorpio_call"] == "nan":
                    sm_to_insert.append(lm)

        # convert to pandas
        logging.info("Preparing data for loading to database ..")
        cm_df = pd.DataFrame.from_records(cm_to_insert)
        if len(cm_df.index) == 0:
            logging.info("Nil to update, finished.")
            # nil to do
            return 0

        logging.info("Starting bulk loading ..")
        sm_df = pd.DataFrame.from_records(sm_to_insert)
        self._bulk_load(cm_df, "clinical_metadata")
        # read the clinical metadata, linking sample_id to cm_int_id
        get_cm_int_ids_sql = self.session.query(
            ClinicalMetadata.sample_id, ClinicalMetadata.cm_int_id
        ).statement
        cm_id_df = pd.read_sql(get_cm_int_ids_sql, con=self.engine)

        sm_df = sm_df.merge(cm_id_df, on="sample_id", how="inner")
        sm_df.drop(["sample_id"], axis=1, inplace=True)
        self._bulk_load(sm_df, "sequence_feature")

        return len(cm_df.index)

    def latest_build_int_id(self, only_if_model_loaded = True):
        """returns the most recent build id
        
        Parameters
        ----------
        only_if_model_loaded: boolean, default True, ignore any entries where the model did not load
        
        Returns
        -------
        None, if there are no builds, or the build_int_id of the latest build
        
        """

        model_loaded_acceptable_values = [1]
        if only_if_model_loaded is False:
            model_loaded_acceptable_values = [0,1]
        (retVal,) = self.session.query(func.max(Build.build_int_id)).filter(Build.model_loaded.in_(model_loaded_acceptable_values)).one_or_none()
        return retVal

    def sequence_features(self, is_lineage):
        """returns sequence_features supplied as part of clinical & sequence metadata.

        Parameters
        ----------
        is_lineage: either 1,0 or None.  If 1, only returns lineages.  If 0, only returns mutations.  If None, returns both.
        
        Returns
        -------
        A list of sequence features.  These are in the form featuretype:feature, eg. pangolineage:B1.1.7"""
        search_lineages = [1, 0]
        if is_lineage is not None:
            search_lineages = [is_lineage]

        retVal = [
            x
            for x, in self.session.query(SequenceFeature.sequencefeature)
            .filter(SequenceFeature.sequencefeature is not None)
            .filter(SequenceFeature.is_lineage.in_(search_lineages))
            .group_by(SequenceFeature.sequencefeature)
        ]
        return retVal

    def make_contingency_tables(self, only_pc_cats_less_than_days_old=120, build_int_id = None):
        """makes contingency tables, computing the relationships between various features and pc_cats.
        
        Parameters:
        -----------
        only_pc_cats_less_than_days_old: only reports on pc_cats which have recently emerged
        build_int_id: make report for build_int_id.  If none, reports on the latest build.
        
        Returns:
        --------
        None
        
        
        """
        logging.info("Making contingency tables")


        # for the latest build
        if build_int_id is not None:
            lbii = build_int_id
        else:
            lbii = self.latest_build_int_id()

        # if there are no builds, we return.
        if lbii is None:
            logging.info("No builds.  Cannot produce contingency tables. ")
            return 

        # test where already entered
        n_feature_associations_present, = self.session.query(func.count(FeatureAssociation.featassoc_int_id)).filter(FeatureAssociation.build_int_id == lbii).one()
        if n_feature_associations_present > 0:
            logging.info("Features already stored for build {0}".format(lbii))
            return 

        # we compute the total number of samples used in the model.
        logging.info("Computing total samples used in the model")
        (a_b_c_d,) = (
            self.session.query(func.count(Sample.sample_int_id))
            .join(ClinicalMetadata, Sample.sample_id == ClinicalMetadata.sample_id)
            .filter(Sample.build_int_id == lbii)
            .filter(Sample.used_in_pca == 1)
            .one_or_none()
        )
        logging.info("Computing sequence feature counts (marginal totals) used in the model")
        a_c_sql = (
            self.session.query(
                SequenceFeature.sequencefeature,
                func.count(Sample.sample_int_id.distinct()).label("a_c"),
                func.min(ClinicalMetadata.sample_date).label("earliest_sample_date"),
                func.max(ClinicalMetadata.sample_date).label("latest_sample_date"),
            )
            .join(
                ClinicalMetadata,
                SequenceFeature.cm_int_id == ClinicalMetadata.cm_int_id,
            )
            .join(Sample, Sample.sample_id == ClinicalMetadata.sample_id)
            .filter(Sample.build_int_id == lbii)
            .filter(Sample.used_in_pca == 1)
            .filter(SequenceFeature.is_lineage == 1)
            .group_by(SequenceFeature.sequencefeature)
            .statement
        )
        a_c = pd.read_sql(a_c_sql, con=self.engine)
        a_c["date_now"] = datetime.date.today()
        a_c["diff_days"] = a_c["date_now"] - a_c["latest_sample_date"]
        a_c["diff_days"] = a_c["diff_days"] / np.timedelta64(1, "D")

        # if appropriate, restrict to lineages seen in last epoch.
        if only_pc_cats_less_than_days_old is not None:
            a_c = a_c[a_c["diff_days"] <= only_pc_cats_less_than_days_old]

        logging.info("Computing pc_cat counts (marginal totals) used in the model")
        a_b_sql = (
            self.session.query(
                TransformedCoordinateCategory.pc_cat,
                func.count(Sample.sample_int_id.distinct()).label("a_b"),
                func.min(ClinicalMetadata.sample_date).label("earliest_sample_date"),
            )
            .join(
                Sample,
                TransformedCoordinateCategory.sample_int_id == Sample.sample_int_id,
            )
            .join(ClinicalMetadata, ClinicalMetadata.sample_id == Sample.sample_id)
            .filter(Sample.build_int_id == lbii)
             .filter(Sample.used_in_pca == 1)
            .group_by(TransformedCoordinateCategory.pc_cat)
            .statement
        )
       
        a_b = pd.read_sql(a_b_sql, con=self.engine)
        logging.info("Recovered details of {0} pc_cats".format(len(a_b.index)))
        
        a_b["date_now"] = datetime.date.today()
        a_b["diff_days"] = a_b["date_now"] - a_b["earliest_sample_date"]
        a_b["diff_days"] = a_b["diff_days"] / np.timedelta64(1, "D")
        
        if only_pc_cats_less_than_days_old is not None:
            a_b = a_b[a_b["diff_days"] <= only_pc_cats_less_than_days_old]

        logging.info(
            "Computing associations all for {1} pc_cats, from which results cat be extracted for {0} pairs, or {2} comparisons".format(
                len(a_c.index), len(a_b.index), len(a_b.index) * len(a_c.index)
            )
        )
        n_added = 0

        if self.show_bar:
            bar = progressbar.ProgressBar(max_value=len(a_b.index))

        for i, a_b_ix in enumerate(a_b.index):
            this_pc_cat = a_b.at[a_b_ix, "pc_cat"]
            this_a_b = a_b.at[a_b_ix, "a_b"]
            if self.show_bar:
                bar.update(i)
            
            # compute associations for the whole selected sequencefeature x pc_cat matrix.
            a_sql = (
                self.session.query(
                    SequenceFeature.sequencefeature,
                    func.count(Sample.sample_int_id.distinct()).label("a")
                )
                .join(
                    ClinicalMetadata, ClinicalMetadata.sample_id == Sample.sample_id
                )
                .join(
                    SequenceFeature,
                    SequenceFeature.cm_int_id == ClinicalMetadata.cm_int_id,
                )
                .join(
                    TransformedCoordinateCategory,
                    TransformedCoordinateCategory.sample_int_id
                    == Sample.sample_int_id,
                )
                .filter(Sample.build_int_id == lbii)
                .filter(Sample.used_in_pca == 1)
                .filter(SequenceFeature.is_lineage == 1)
                .filter(TransformedCoordinateCategory.pc_cat == this_pc_cat)
                .group_by(
                    SequenceFeature.sequencefeature,
                )
                .statement
            )
            a_df = pd.read_sql(a_sql, con=self.engine)

            aca = a_c.merge(a_df, how='left', on='sequencefeature')
           
            aca['a'].fillna(0, inplace=True)       
            aca['pc_cat'] = this_pc_cat         
            for a_c_ix in aca.index:
                this_a_c = aca.at[a_c_ix, "a_c"]
                this_a =aca.at[a_c_ix, "a"]
                this_sequencefeature = aca.at[a_c_ix, "sequencefeature"]
                this_c = this_a_c - this_a
                this_b = this_a_b - this_a
                this_d = a_b_c_d - (this_a + this_b + this_c)

                n_added += 1
                ct = ContingencyTable(
                    lbii,
                    this_pc_cat,
                    this_sequencefeature,
                    this_a,
                    this_b,
                    this_c,
                    this_d
                )

                self.session.add(ct.featureassociation())

                if n_added % 1 == 0:
                    self.session.commit()
        if self.show_bar:
            bar.finish()
        logging.info("Computing associations completed")
        self.session.commit()

    def counts_per_day_per_pc_cat(self, build_int_id=None, before_date=None):
        """returns counts per day, by
        - pc_cat
        - nation

        Parameters
        ----------
        build_int_id: the build_int_id of a particular build.  If None, uses the latest build
        before_date:  don't report dates after this date.  If none, uses today's date

        Returns
        -------
        A pandas dataframe containing the count data
        Suitable for subsequent linear modelling.  
        """

        # determine the latest date
        if date_end is None:
            date_end = datetime.datetime.now()

        # determine the 
        if build_int_id is not None:
            lbii = build_int_id
        else:
            lbii = self.latest_build_int_id()

        ## INCOMPLETE

    def fn4_bulk_upload(self, dumpdir):
        """loads the files, as produced by findNeighbour4_dumpneighbours.py,
        into tables.

        Note that this function is intended for mass upload of data from a running findNeighbour4 server.
        It is not suitable for 'dripfeeding' new links.  To do this latter operation, please see function
        XXXXXXXXXXXXXX [TODO]

        Suitable SQL to recover loaded links is as below.

        -------------------------------------------------------------------------------------------------
        select s1.sample_id sample_id_1, s2.sample_id sample_id_2, n.dist from pcadaemon.fn4_sample s1
        inner JOIN
        pcadaemon.neighbour n
        on s1.fn4s_int_id = n.fn4s_int_id_1

        inner JOIN

        pcadaemon.fn4_sample s2
        on s2.FN4S_INT_ID = n.fn4s_int_id_2

        where s1.sample_id = 'SGUH-12E89'; 
        --------------------------------------------------------------------------------------------------
        """

        # check whether there is existing data.
        # if so, it's not appropriate to use this function
        n_already, = self.session.query(func.count(FN4Sample.fn4s_int_id)).one()
        if n_already > 0:
            raise PCADBManagerError("Cannot do bulk upload of fn4 data if there is data already present")

        # load the samples file
        samples_file = os.path.join(dumpdir, 'samples.json')
        with open(samples_file, 'rt') as f:
            samples_dict = json.load(f)
        samples_df = pd.DataFrame.from_dict(samples_dict, orient='index')
        samples_df.columns = ['is_invalid']
        samples_df['sample_id'] = samples_df.index


        # load the samples into the database
        self._bulk_load(samples_df, 'fn4_sample')

        # recover the integer ids for each sample_id
        sample_sql = self.session.query(FN4Sample).statement
        samples_id_df = pd.read_sql(sample_sql, self.engine)
        
        samples_id_df_1 = samples_id_df.drop(['is_invalid'], axis =1)
        samples_id_df_1.columns = ['fn4s_int_id_1','sample_id_1']
        samples_id_df_2 = samples_id_df.drop(['is_invalid'], axis =1)
        samples_id_df_2.columns = ['fn4s_int_id_2','sample_id_2']
        

        n_links = 0
        neighbours_files = glob.glob(os.path.join(dumpdir, 'neighbours_*.json'))
      
        for i, neighbour_file in enumerate(neighbours_files):
            
            logging.info("Loading file {0}/{1}".format(i,len(neighbours_files)))
            with open(neighbour_file, 'rt') as f:
                neighbours = json.load(f)

            if len(neighbours)>0:       # something to load
                neighbours_df = pd.DataFrame.from_records(neighbours)
                neighbours_df.columns = ['sample_id_1', 'sample_id_2', 'dist']
                initial_n = len(neighbours_df.index)
                
                # merge in the integer ids
                neighbours_df = neighbours_df.merge(samples_id_df_1, how='inner', on='sample_id_1')
                neighbours_df = neighbours_df.merge(samples_id_df_2, how='inner', on='sample_id_2')
                post_merge_n = len(neighbours_df.index)
                assert(initial_n == post_merge_n)       # if this is not the case, then there not all the neighbours are in the samples_df, which is not allowed and indicates a software/database issue

                to_load_df = neighbours_df.drop(['sample_id_1','sample_id_2'], axis = 1)
                self._bulk_load(to_load_df, 'neighbour')
                n_links = n_links + len(to_load_df.index)
        
        logging.info("Load completed.  Loaded a total of {0} links.".format(n_links))


        