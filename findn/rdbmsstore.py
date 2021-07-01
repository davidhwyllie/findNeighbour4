#!/usr/bin/env python
""" provides a storage layer for meta-data and snv distances from the
findneighbour4 system in a RDBMS

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
import time
import os
import json
import pandas as pd
import logging

import numpy as np
import warnings
import progressbar
from sqlalchemy import (
    Integer,
    Column,
    Float,
    MetaData,
    Text,
    String,
    DateTime,
    Identity,
    ForeignKey,
    TIMESTAMP
)

#from sqlalchemy.dialects import oracle
from sqlalchemy import create_engine, inspect
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
import cx_Oracle

# global: definition of database structure
# classes mapping to persistence database inherit from this
db_pc = declarative_base()  # global
metadata = MetaData()


class RDBMSError(Exception):
    """a general purpose error used by the rdbms module.

    ## NOTE: not clear whether this needed at present"""

    pass


## define schema
class BulkLoadTest(db_pc):
    """used only for testing bulk uploads as part of unit testing"""

    __tablename__ = "fn4_bulk_load_test"
    blt_int_id = Column(Integer, Identity(start=1), primary_key=True)
    bulk1 = Column(Integer)
    bulk2 = Column(Integer)

class Config(db_pc):
    """stores config data"""

    __tablename__ = "config"
    cfg_int_id = Column(Integer, Identity(start=1), primary_key=True)
    config_key = Column(String(56), index= True, unique = True)
    config_value = Column(Text)     # blob

class RefCompressedSeq(db_pc):
    """stores reference compressed sequences, which are large character objects

    Note: the mongodb equivalent is the GridFS meta-collection refcompressedseq.
    Note: the mongodb equivalent stores the reference compressed sequence objects pickled.
    Note:  There is no good reason to do this, apart from speed (not tested).  JSON object storage is completely acceptable

    For details see:
    https://docs.mongodb.com/manual/core/gridfs/
    """

    __tablename__ = "refcompressedseq"
    seq_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    sequence_id = Column(
        String(38),
        index=True,
        unique=True,
        comment="the sample_id represented by the entry; sample_ids are typically guuids",
    )
    examination_date = Column(
        TIMESTAMP,
        index=True,
        comment="the date and time the record was examined and compressed",
    )
    annotations = Column(
        Text, comment="a json string, representing metadata about the sequence"
    )
    invalid = Column(
        Integer,
        index=True,
        comment="whether the sequence is of sufficient quality to be analysed (invalid = 0) or is not (invalid = 1).  Part of the annotations, but extracted into separate field for indexing.",
    )
    prop_actg = Column(
        Float,
        index=True,
        comment="the proportion of A,C,G,T (as opposed to N,-, or IUPAC codes).  Part of the annotations, but extracted into separate field for indexing.",
    )
    content = Column(
        Text, comment="a json string, representing the reference compressed sequence"
    )


class Edge(db_pc):
    """represents a pair of RefCompressedSeq which are similar.  Table is called 'edge' because it represents an edge in a network in which the vertices are RefCompressedSequences.
    Note
    - that an edge A -> B is represented twice in this data representation: as A -> B, and as B -> A.
    - This means that all the edges of A can be obtained by
    statements such as SELECT * from Edge where seq_int_id_1 = 217, where 217 is the seq_int_id of sequence A."""

    __tablename__ = "edge"
    edge_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    seq_int_id_1 = Column(
        Integer,
        ForeignKey(RefCompressedSeq.seq_int_id),
        comment="One of a pair of sequences. Note that for rapid loading, disabling of the FK could be considered.",
    )
    seq_int_id_2 = Column(
        Integer,
        comment="One of a pair of sequences.  Note: foreign key constraint not enforced at a database level ",
    )
    dist = Column(Integer, comment="the SNV distance between sequences")


class Cluster(db_pc):
    """stores clusters, which are large character objects (json)"""

    __tablename__ = "cluster"
    cl_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    upload_date = Column(
        DateTime, index=True, comment="the time the record was inserted"
    )
    cluster_build_id = Column(
        String(40),
        index=True,
        unique=True,
        comment="an identifier for the contents; this is typically and sha-1 hash on the contents",
    )
    upload_date = Column(
        DateTime, index=True, comment="the time the record was inserted"
    )
    content = Column(Text, comment="a json string describing the cluster")


class Monitor(db_pc):
    """stores server monitoring entries, which are large character objects (json)"""

    __tablename__ = "monitor"
    mo_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    upload_date = Column(
        DateTime, index=True, comment="the time the record was inserted"
    )
    content = Column(Text, comment="a json dictionary including ")


class MSA(db_pc):
    """stores multisequence alignments, which are large character objects (json)"""

    __tablename__ = "msa"
    msa_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    msa_id = Column(
        String(40),
        index=True,
        unique=True,
        comment="an identifier for the contents; this is typically and sha-1 hash on the contents",
    )
    upload_date = Column(
        DateTime, index=True, comment="the time the record was inserted"
    )
    content = Column(Text, comment="character large object containing a json string")


class TreeStorage(db_pc):
    """stores trees and related data, which are large character objects (json)"""

    __tablename__ = "tree"
    ts_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    ts_id = Column(
        String(40),
        index=True,
        unique=True,
        comment="an identifier for the contents; this is typically and sha-1 hash on the contents",
    )

    upload_date = Column(
        DateTime, index=True, comment="the time the record was inserted"
    )
    content = Column(Text, comment="character large object containing a json string")


class NPEncoder(json.JSONEncoder):
    """encodes Numpy types as jsonisable equivalents"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NPEncoder, self).default(obj)


class fn3persistence:
    """System for persisting results from  large numbers of sequences stored in FindNeighbour 3+.
    Uses a generic rdbms, with optimisations for Oracle databases when using the cx_oracle package.

    the associated schema is defined using SQLalchemy tables, as above.

    Note that parts of this code are taken from the class pcadb.PCADatabaseManager.
    In future, it may be possible to integrate this class with that, if successful
    implementation of findNeighbour4 functionality using an RDBMS is possible


    Note:
    #################################################################################################################################
    to what extent the api for this class can be made to exactly match the mongodb fn3persistence object deserves thought,
    because if it can't, then a lot of breaking changes are going to have to be made and mongo either removed or a fork maintained.
    Doing this is highly undesirable

    Mongodb persistence object creation method looks like

        CONFIG OF MONGOSTORE
            # code handling startup and shutdown.
            def __init__(
                self,
                connString,
                dbname="fn3_unittesting",
                debug=0,
                config_settings={},
                max_neighbours_per_document=100000,
                server_monitoring_min_interval_msec=0,
            ):
                Creates a connection to a MongoDb database.

                connString : the mongoDb connection string
                dbname: the name of the mongoDb database to use.
                if debug = 0 or 1, the database is opened or created.
                if debug = 2, any existing collections are deleted.
                config_settings: only used on db creation; optional dictionary to note items in the database's config collection.
    ##################################################################################################################################
    this work has not yet been carried out.

    """

    def __init__(self, connection_config=None, debug=0):

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
            self.engine_name = "sqlite://"
        elif "://" in connection_config:
            logging.info(
                "Connection config provided; using {0}".format(connection_config)
            )
            self.engine_name = connection_config
        else:
            # PCA_CONNECTION_CONFIG_FILE should contain credentials
            conn_detail_file = None
            try:
                conn_detail_file = os.environ["PCA_CONNECTION_CONFIG_FILE"]
            except KeyError:
                raise RDBMSError(
                    "Environment variable PCA_CONNECTION_CONFIG_FILE does not exist; however, it is required.  If you are using a python virtual environment, you need to set it in .env, not globally"
                )

            if conn_detail_file is None:
                # we failed to set it
                raise RDBMSError(
                    "Tried to set conn_detail_file from environment variable PCA_CONNECTION_CONFIG_FILE, but it is still None."
                )

            if not os.path.exists(conn_detail_file):
                raise FileNotFoundError(
                    "Connection file specified but not found: {0}".format(
                        conn_detail_file
                    )
                )

            # read the config file
            with open(conn_detail_file, "rt") as f:
                conn_detail = json.load(f)

            if connection_config not in conn_detail.keys():
                raise RDBMSError(
                    "Connection {0} does not correspond to one of the keys {1} of the configuration json file at {2}".format(
                        connection_config, conn_detail.keys(), conn_detail_file
                    )
                )

            # configure engine
            this_configuration = conn_detail[
                connection_config
            ]  # extract the relevant part of the config dictionary

            # two keys are always present
            essential_keys = set(["DBTYPE", "ENGINE_NAME"])
            if len(essential_keys - set(this_configuration.keys())) > 0:
                raise RDBMSError(
                    "Provided keys for {0} are not correct.  Required are {1}".format(
                        connection_config, essential_keys
                    )
                )

            # if it's Oracle, then three keys are required.
            if this_configuration["DBTYPE"] == "oracle":
                essential_keys = set(["DBTYPE", "ENGINE_NAME", "TNS_ADMIN"])
                if len(essential_keys - set(this_configuration.keys())) > 0:
                    raise RDBMSError(
                        "Provided keys for oracle db in {0} are not correct.  Required are {1}".format(
                            connection_config, essential_keys
                        )
                    )

                # set the TNS_ADMIN variable.
                logging.info(
                    "Set TNS_ADMIN to value specified in config file {0}".format(
                        this_configuration["TNS_ADMIN"]
                    )
                )
                os.environ["TNS_ADMIN"] = this_configuration["TNS_ADMIN"]

            logging.info("Set ENGINE_NAME configuration string from config file.")
            self.engine_name = this_configuration["ENGINE_NAME"]

        # now we can start
        self.Base = db_pc
        logging.info("DatabaseManager: Connecting to database")
        self.engine = create_engine(self.engine_name)
        self.is_oracle = "oracle+cx" in self.engine_name
        self.is_sqlite = "sqlite://" in self.engine_name
        self.show_bar = True  # maybe define a method to switch this off

        logging.info(
            "DatabaseManager: Database connection made; there are {0} tables.  Oracle database = {1}".format(
                len(self._table_names()), self.is_oracle
            )
        )

        # drop existing tables if in debug mode
        if debug == 2:
            print("Dropping existing tables")
            self._drop_existing_tables()

        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.Base.metadata.create_all(bind=self.engine)  # create the table(s)
        self.session.commit()

    def no_progressbar(self):
        """don't use progress bars"""
        self.show_bar = False

    def _table_names(self):
        """returns table names in the schema"""
        return inspect(self.engine).get_table_names()

    def _drop_existing_tables(self):
        """drops any existing tables"""

        # test - wait 10 seconds - for any existing transactions to cleartherwise, ORA-00054 can occur during oracle testing
        # happens only with multiple cpus
        if self.is_oracle:
            time.sleep(0)

        self.Base.metadata.create_all(
            self.engine
        )  # create the table(s) if they don't already exist
        BulkLoadTest.__table__.drop(self.engine)
        Config.__table__.drop(self.engine)
        Edge.__table__.drop(self.engine)
        Cluster.__table__.drop(self.engine)
        Monitor.__table__.drop(self.engine)
        MSA.__table__.drop(self.engine)
        TreeStorage.__table__.drop(self.engine)
        RefCompressedSeq.__table__.drop(self.engine)

        remaining = len(self._table_names())
        if remaining > 0:
            warnings.warn(
                "Failed to remove all tables in the database.  Is this database being used by another program?  The following remain: {0}".format(
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

        The maximum size that can be transmitted to the Oracle server is 2G.
        If this happens, a cx_Oracle.DatabaseError
        is raised with result DPI-1015.  R Code to prevent this by auto-setting max_batch is provided.

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
            raise TypeError(
                "upload_df needs to be pandas DataFrame, not a {0}".format(
                    type(upload_df)
                )
            )

        # check that we have been passed data
        ncol = len(upload_df.columns)
        if ncol == 0:
            # we treat this as an error
            raise RDBMSError(
                "Passed data frame to _bulk_upload to {0} contains no columns {1}".format(
                    target_table, upload_df
                )
            )

        logging.info("Bulk upload to {0} started".format(target_table))

        if self.is_sqlite:
            # there is a max variable limit of 32,766 for Sqlite 3.32.0 on https://www.sqlite.org/limits.html
            # set max_batch to keep below this.
            max_batch = int(32000 / ncol)
            logging.info(
                "Autoset max_batch to {0}, as running SQLite".format(max_batch)
            )

        if self.is_oracle:
            # ************* commit via cx_Oracle bulk upload syntax ****************
            # parse the engine_name into dsn, database & password
            e1 = self.engine_name.replace("oracle+cx_oracle://", "")
            up, dsn = e1.split("@")
            u, p = up.split(":")

            # get into the right format for loading: note: holds all data in ram
            loadvar = list(upload_df.itertuples(index=False, name=None))
            ncol = len(upload_df.columns.to_list())

            # estimate the maximum buffer size required and auto-reduce the max_batch if required.
            estimated_row_size = len(str(loadvar[0]))
            estimated_max_batch = int(
                1e7 / (estimated_row_size)
            )  # large margin of safety 10M batch size
            if estimated_max_batch < max_batch:
                max_batch = estimated_max_batch
                logging.info(
                    "Reduced max_batch to keep estimated buffer size within acceptable limits (<= 10M target).  max_batch is {0}".format(
                        max_batch
                    )
                )

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

            # *****************************  ALL DATABASES OTHER THAN ORACLE ************************
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
                )  # pandas method
                upload_df = upload_df.iloc[max_batch:]
                if self.show_bar:
                    bar.update(start_n - len(upload_df.index))

        logging.info("Bulk upload to {0} complete".format(target_table))
        if self.show_bar:
            bar.finish()
        return len(upload_df.index)
