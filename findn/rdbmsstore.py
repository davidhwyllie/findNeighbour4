#!/usr/bin/env python
""" provides a storage layer for meta-data and snv distances from the
findneighbour4 system in a RDBMS

Tested with:
- Oracle ADW (cloud service)
- Sqlite (but, sqlite can't be used as a server backend)

Not tested:
MS SQL server, PostgreSQL

Tested but doesn't work at present
MySQL - the issue is with storage of character large objects in TEXT fields.  The SQL alchemy version used make TEXT, not LARGETEXT fields.
      - Workarounds described https://stackoverflow.com/questions/47644739/what-column-type-does-sqlalchemy-use-for-text-on-mysql did not fix this
      - probably a soluble problem 
      - tested with MySQL 8 on Ubuntu 20.  Connection string was "mysql://root:root@localhost:3306/test_db" with user/password root

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""
import bson  # type: ignore
from datetime import datetime, timedelta, date
import hashlib
import os
import json
import pandas as pd
import psutil
import logging
import numpy as np
import warnings
import uuid
import cx_Oracle
from sentry_sdk import capture_message, capture_exception
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
    Index,
    TIMESTAMP,
    func,
    create_engine,
    inspect,
)
from findn.seq2json import SeqDictConverter
from sqlalchemy.sql.expression import desc
from sqlalchemy.orm import (
    sessionmaker,
    scoped_session,
)
from sqlalchemy.ext.declarative import declarative_base
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    NoReturn,
    Optional,
    Set,
    TypedDict,
    Union,
)

Guid2NeighboursFormat1 = List[Union[str, int]]
Guid2NeighboursFormat3 = Union[str]
Guid2NeighboursFormat4 = Dict[str, Union[str, int]]
Guid2NeighboursFormats = Union[
    Guid2NeighboursFormat1, Guid2NeighboursFormat3, Guid2NeighboursFormat4
]


class RecentDatabaseMonitoringRet(TypedDict, total=False):
    recompression_data: bool
    latest_stats: Dict[str, Union[int, np.float64]]
    trend_stats: List[Dict[str, Any]]


# global: definition of database structure
# classes mapping to persistence database inherit from this
db_pc = declarative_base()  # global
metadata = MetaData()


class RDBMSError(Exception):
    """a general purpose error used by the rdbms module."""

    pass


## define schema
class BulkLoadTest(db_pc):
    """used only for testing bulk uploads as part of unit testing"""

    __tablename__ = "fn4_bulk_load_test"
    blt_int_id = Column(Integer, Identity(start=1), primary_key=True)
    bulk1 = Column(Integer)
    bulk2 = Column(Integer)


class FNLock(db_pc):
    """used for storing details of one or more classes of lock"""

    __tablename__ = "fn4lock"
    lock_int_id = Column(
        Integer,
        primary_key=True,
        comment="an integer reflecting the kind of lock studied",
    )
    lock_set_date = Column(
        TIMESTAMP, index=True, comment="the date and time the lock was modified"
    )
    uuid = Column(String(32))
    lock_status = Column(
        Integer, comment="whether the lock is in place (1) or not in place (0)"
    )


class Config(db_pc):
    """stores config data"""

    __tablename__ = "config"
    cfg_int_id = Column(Integer, Identity(start=1), primary_key=True)
    config_key = Column(String(56), index=True, unique=True)
    config_value = Column(Text(50000000))  # 50M limit.  Required for Mysql


class RefCompressedSeq(db_pc):
    """stores reference compressed sequences, which are large character objects, and their annotations.

    Note: the mongodb equivalent is the GridFS meta-collection refcompressedseq and the standard collection guid2meta.
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
        comment="the sample_id represented by the entry; sample_ids are typically guids",
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
        default=-1,
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
    statements such as SELECT * from Edge where seq_int_id_1 = 217, where 217 is the seq_int_id of sequence A.

    - if insertion is fast enough, we could enable foreign key constraints here.
    - it is likely that insert speed will be the main determinant of server speed
    - and the faster the inserts work the better.
    - in the mongo implementation, FK constraints are not implemented, and the relationships are guaranteed by application logic, not at a database level .
    At present, this approach is being used here too.

    """

    __tablename__ = "edge"
    edge_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    sequence_id_1 = Column(
        String(38),
        comment="One of a pair of sequences. Note: foreign key constraint not enforced at a database level",
    )
    sequence_id_2 = Column(
        String(38),
        comment="One of a pair of sequences.  Note: foreign key constraint not enforced at a database level ",
    )
    dist = Column(Integer, comment="the SNV distance between sequences")


Index("ix_Edge_1", Edge.sequence_id_1, unique=False, oracle_compress=1)


class Cluster(db_pc):
    """stores clusters, which are large character objects (json)"""

    __tablename__ = "sequence_cluster"
    cl_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    cluster_build_id = Column(
        String(40),
        index=True,
        comment="an identifier for the contents; this is typically a sha-1 hash on the contents",
    )
    upload_date = Column(
        DateTime, index=True, comment="the time the record was inserted"
    )
    content = Column(Text, comment="a json string describing the cluster")


class Monitor(db_pc):
    """stores monitor entries, which are large character objects (html)"""

    __tablename__ = "monitor"
    mo_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    mo_id = Column(
        String(40),
        index=True,
        unique=True,
        comment="an identifier for the contents; this is typically and sha-1 hash on the contents",
    )
    content = Column(Text, comment="html data")


class ServerMonitoring(db_pc):
    """stores server monitoring entries, which are large character objects (json)"""

    __tablename__ = "server_monitoring"
    sm_int_id = Column(
        Integer,
        Identity(start=1),
        primary_key=True,
        comment="the primary key to the table",
    )
    sm_id = Column(
        String(40),
        index=True,
        unique=True,
        comment="an identifier for the contents; this is typically and sha-1 hash on the contents",
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
        String(60),
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
    """encodes Numpy  and datetime types as jsonisable equivalents"""

    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, date):
            return obj.isoformat()
        elif isinstance(obj, datetime):
            return obj.isoformat()
        else:
            return super(NPEncoder, self).default(obj)


class fn3persistence_r:
    """System for persisting results from  large numbers of sequences stored in FindNeighbour 3+.
    Uses a generic rdbms, with optimisations for Oracle databases when using the cx_oracle package.

    the associated schema is defined using SQLalchemy tables, as above.

    Note that parts of this code are taken from the class pcadb.PCADatabaseManager.
    In future, it may be possible to integrate this class with that, if successful
    implementation of findNeighbour4 functionality using an RDBMS is possible

    See also the Persistence class in the persistence module.
    This provides an API which will either use this class, or the mongo based fn3persistence class, depending on
    software settings.

    """

    def __init__(
        self, connection_config, debug=0, server_monitoring_min_interval_msec=0
    ):

        """creates the RDBMS connection

        Parameters
        -----------
        connection_config:
        One of
        1. a key to a dictionary containing one or more database configuration details: (e.g. 'prod', 'test')
        2. a valid sqlalchemy database connection string (if this is sufficient for connections)  e.g. 'pyodbc+mssql://myserver'
        3. None.  This is considered to mean 'sqlite://' i.e. an in memory sqlite database, which is not persisted when the program stops.

        if it is not none, a variable called DB_CONNECTION_CONFIG_FILE must be present.  This must point to a file containing credentials.
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

        self.sjc = SeqDictConverter()
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.INFO)

        self.debug = debug
        self.storage_technology = "rdbms"

        self.logger.info("Storage technology is {0}".format(self.storage_technology))

        self.server_monitoring_min_interval_msec = server_monitoring_min_interval_msec
        self.previous_server_monitoring_data = {}
        self.previous_server_monitoring_time = None

        # connect and create session.  Validate inputs carefully.
        if connection_config is None:
            self.logger.info("Connection config is None: using in-memory sqlite.")
            self.engine_name = "sqlite://"

        elif "://" in connection_config:
            # it's not None, and we assume what we are provided is an sqlalchemy database connection string
            self.logger.info(
                "Connection config provided; using {0}".format(connection_config)
            )
            self.engine_name = connection_config
        else:
            # we have been passed a token.  this should be a key to a dictionary, stored in
            # DB_CONNECTION_CONFIG_FILE which contains credentials
            conn_detail_file = None
            try:
                conn_detail_file = os.environ["DB_CONNECTION_CONFIG_FILE"]
            except KeyError:
                raise RDBMSError(
                    "Environment variable DB_CONNECTION_CONFIG_FILE does not exist; however, it is required.  If you are using a python virtual environment, you need to set it in .env, not globally"
                )

            if conn_detail_file is None:
                # we failed to set it
                raise RDBMSError(
                    "Tried to set conn_detail_file from environment variable DB_CONNECTION_CONFIG_FILE, but it is still None."
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
                self.logger.info(
                    "Set TNS_ADMIN to value specified in config file {0}".format(
                        this_configuration["TNS_ADMIN"]
                    )
                )
                os.environ["TNS_ADMIN"] = this_configuration["TNS_ADMIN"]

            self.logger.info("Set ENGINE_NAME configuration string from config file.")
            self.engine_name = this_configuration["ENGINE_NAME"]

        self.using_sqlite = self.engine_name.startswith("sqlite://")

        # now we can start
        self.Base = db_pc
        self.logger.info("DatabaseManager: Connecting to database")
        self.is_oracle = "oracle+cx" in self.engine_name
        self.is_sqlite = "sqlite://" in self.engine_name
        self.show_bar = True  # maybe define a method to switch this off

        # create engine
        self.engine = create_engine(self.engine_name)  # sqlalchemy generic pool manager

        # oracle pool manager code
        # use cx_Oracle pool manager
        # u, p, dsn = self.oracle_connstring_parts(self.engine_name)
        # self.oracle_pool = cx_Oracle.SessionPool(
        #    user = u,
        #    password = p,
        #    dsn = dsn,
        #   min = 4,
        #    max = 4,
        #    encoding="UTF-8",
        #    nencoding="UTF-8"
        # )
        # self.engine = create_engine("oracle://", creator = self.oracle_pool.acquire, poolclass = NullPool)

        self.logger.info(
            "DatabaseManager: Database connection made; there are {0} tables.  Oracle database = {1}".format(
                len(self._table_names()), self.is_oracle
            )
        )

        self.Base.metadata.create_all(
            bind=self.engine
        )  # create the table(s) if they don't exist

        # create thread-local sessions, see https://docs.sqlalchemy.org/en/14/orm/contextual.html#using-custom-created-scopes
        session_factory = sessionmaker(bind=self.engine)

        # this object will call session factory when we create/request a thread local session
        # e.g. thread_local_session = self.Session()
        self.Session = scoped_session(session_factory)

        # drop existing tables if in debug mode
        # delete any pre-existing data if we are in debug mode.
        if debug == 2:
            logging.warning(
                "Debug mode operational [DEBUG={0}]; deleting all data from tables.".format(
                    debug
                )
            )
            self._delete_existing_clustering_data()
            self._delete_existing_data()

        else:
            self.logger.info("Using stored data in rdbms")

    def closedown(self):
        """closes the session(s) & disposes of any engine.
        Is required for unit testing"""
        try:
            self.engine.dispose()
        except Exception:
            pass  # KeyError and other errors can occur if dispose() has already been called automatically

    def no_progressbar(self):
        """don't use progress bars"""
        self.show_bar = False

    def _table_names(self):
        """returns table names in the schema.
        If the schema's contents have not been created, returns an empty list"""
        return inspect(self.engine).get_table_names()

    def thread_local_session(self, n_retries=3, simulate_failure="no", log_errors=True):
        """generates, or selects a thread local session from a session factory or pool.

        For context, see https://docs.sqlalchemy.org/en/13/orm/contextual.html

        Checks that the session recovered from the session pool is still valid, since they can time out.
        If it is timed out, tries another.  Will retry up to n_retries times.

        Parameters:
        n_retries:  the number of attempts which will be made to generate a functional connection.
        simulate_failure: simulates the failure of a connection (closes the connection before returning it) - used only for unittesting.  Valid values:
                'no' : normal operation
                'once': fail once
                'always': fail every time, even on repeated attempts to connect
        log_errors: logs any errors to sentry & error log (recommended)"""

        tries = n_retries
        while tries > 0:
            tries = tries - 1
            tls = self.Session()

            # simulate failure if required to do so
            if (simulate_failure == "once" and tries - 1 == n_retries) or (
                simulate_failure == "always"
            ):
                tls = (
                    None  # not a session object.  Attempts to use it as such will fail.
                )

            # test whether it is working
            try:
                tls.query(Config).filter_by(
                    config_key="config"
                ).first()  # try to connect

                # if execution continues here, the session works
                return tls

            except Exception as e:
                logging.info(
                    "Failed to connect on trial {0}/{1}".format(
                        n_retries - tries, n_retries
                    )
                )
                if log_errors:
                    capture_exception(e)
                    logging.error(e)
                try:
                    tls.remove()
                except Exception as e:
                    capture_message("Failed to remove session")
                    logging.info("Failed to remove session")
                    logging.error(e)

        raise RDBMSError(
            "Could not connect to database.  Tried {0} times with different sessions".format(
                n_retries
            )
        )

    def _drop_existing_tables(self):
        """empties, but does not drop, any existing tables"""

        self.Base.metadata.create_all(
            self.engine
        )  # create the table(s) if they don't already exist

        BulkLoadTest.__table__.drop(self.engine)
        Config.__table__.drop(self.engine)
        Edge.__table__.drop(self.engine)
        Cluster.__table__.drop(self.engine)
        Monitor.__table__.drop(self.engine)
        ServerMonitoring.__table__.drop(self.engine)
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

        return

    def oracle_connstring_parts(self, connstring):
        """splits an oracle connection string into username, password and DSN"""
        if connstring.startswith("oracle+cx_oracle://"):
            e1 = self.engine_name.replace("oracle+cx_oracle://", "")
            up, dns = e1.split("@")
            u, p = up.split(":")
            return u, p, dns
        else:
            return None

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

        self.logger.info("Bulk upload to {0} started".format(target_table))

        if self.is_sqlite:
            # there is a max variable limit of 32,766 for Sqlite 3.32.0 on https://www.sqlite.org/limits.html
            # set max_batch to keep below this.
            max_batch = int(32000 / ncol)
            self.logger.info(
                "Autoset max_batch to {0}, as running SQLite".format(max_batch)
            )

        if self.is_oracle:
            # ************* commit via cx_Oracle bulk upload syntax ****************
            # parse the engine_name into dsn, database & password

            u, p, dsn = self.oracle_connstring_parts(self.engine_name)

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
                self.logger.info(
                    "Reduced max_batch to keep estimated buffer size within acceptable limits (<= 10M target).  max_batch is {0}".format(
                        max_batch
                    )
                )

            # construct sql statment.
            # Should be injection-safe; we have checked the target_table is a table, and are incorporating integers and verified strings only.
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
                self.logger.info(
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

        self.logger.info("Bulk upload to {0} complete".format(target_table))
        if self.show_bar:
            bar.finish()
        return len(upload_df.index)

    def delete_server_monitoring_entries(self, before_seconds: int) -> None:
        """deletes server monitoring entries more than before_seconds ago"""
        tls = self.thread_local_session()
        earliest_allowed = datetime.now() - timedelta(seconds=before_seconds)
        tls.query(ServerMonitoring).filter(
            ServerMonitoring.upload_date < earliest_allowed
        ).delete()
        tls.commit()
        # finished

    def summarise_stored_items(self) -> Dict[str, Any]:
        """counts how many sequences exist of various types"""
        return {}

    def connect(self) -> None:
        """test whether the database is connected, and if not, tries to connect.
        Does nothing here, just a stub."""
        pass

    def rotate_log(self) -> None:
        """forces rotation of the mongo log file; a stub here, does nothing"""
        pass

    def raise_error(self, token: str) -> NoReturn:
        """raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging"""
        raise ZeroDivisionError(token)

    def _delete_existing_data(self) -> None:
        """deletes existing data from the databases"""

        tls = self.thread_local_session()

        tls.query(Config).delete()
        tls.query(Edge).delete()
        tls.query(RefCompressedSeq).delete()
        tls.query(Monitor).delete()
        tls.query(ServerMonitoring).delete()
        tls.query(BulkLoadTest).delete()
        tls.query(Cluster).delete()
        tls.query(MSA).delete()
        tls.query(TreeStorage).delete()
        tls.query(FNLock).delete()
        tls.commit()
        # finished

    def _delete_existing_clustering_data(self) -> None:
        """deletes any clustering data from the databases"""
        tls = self.thread_local_session()
        tls.query(Cluster).delete()
        tls.query(MSA).delete()
        tls.query(TreeStorage).delete()

    def first_run(self) -> bool:
        """if there is no config entry, it is a first-run situation"""
        tls = self.thread_local_session()
        row = tls.query(Config).filter_by(config_key="config").first()
        return row is None

    def __del__(self) -> None:
        """closes any session"""
        self.closedown()

    def memory_usage(self) -> Dict[str, Union[int, float]]:
        """returns memory usage by current python3 process
        Uses the psutil module, as the resource module is not available in windows.
        """
        memdict = psutil.virtual_memory()._asdict()
        sm = {"server|mstat|" + k: v for k, v in memdict.items()}
        return sm

    # methods for the config collection
    def config_store(self, key: str, object: Dict[str, Any]) -> Any:
        """stores object into config collection
        It is assumed object is a dictionary
        """

        if not isinstance(object, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(object))
            )
        tls = self.thread_local_session()
        if row := tls.query(Config).filter_by(config_key=key).first():
            row.config_value = json.dumps(object, cls=NPEncoder).encode("utf-8")
        else:
            row = Config(
                config_key=key,
                config_value=json.dumps(object, cls=NPEncoder).encode("utf-8"),
            )
            tls.add(row)
        tls.commit()
        # finished

    def config_read(self, key: str) -> Any:
        """loads object from config.
        It is assumed object is a dictionary"""
        tls = self.thread_local_session()
        row = tls.query(Config).filter_by(config_key=key).first()
        if row is None:
            return None
        else:
            return dict(_id=key, **json.loads(row.config_value))

    # methods for the server and database monitoring
    def recent_database_monitoring(
        self, max_reported: int = 100
    ) -> RecentDatabaseMonitoringRet:
        """computes trends in the number of records holding pairs (guid2neighbours) vs. records.
        This ratio is a measure of database health.  Ratios > 100 indicate the database may become very large, and query slowly"""
        return {"recompression_data": False, "latest_stats": {"storage_ratio": 1}}

    def _to_string(self, x=[str, bytes]) -> str:
        """returns a string version of x; if x is bytes, returns a string, assuming utf-8 encoding

        This function is required because some databases return bytes objects
        from 'text' fields (like sqlite) while others return str objects (like oracle)
        """

        if isinstance(x, bytes):
            return x.decode("utf-8")
        else:
            return x

    def recent_server_monitoring(
        self,
        max_reported: int = 100,
        selection_field: Optional[str] = None,
        selection_string: Optional[str] = None,
    ) -> List[dict]:
        """returns a list containing recent server monitoring, in reverse order (i.e. tail first).
        The _id field is an integer reflecting the order added.  Lowest numbers are most recent.

        Inputs
        max_reported - return this number of lines, at most.
        selection_field - if not None, will only return lines containing selection_string
                          in the 'selection_field' key of the returned dictionary.
        selection_string -if selection_field is not None, only returns rows if
                          selection_string is present in the 'selection_field' key of the
                          monitoring element. If None, this constraint is ignored.
        """

        if not isinstance(max_reported, int):
            raise TypeError(
                f"limit must be an integer, but it is a {type(max_reported)}"
            )
        if not max_reported >= 0:
            raise ValueError("limit must be more than or equal to zero")

        if max_reported == 0:
            return []

        def row_to_dict(res: ServerMonitoring) -> dict:
            d = json.loads(res.content)
            d["_id"] = res.sm_int_id
            if selection_field is None:
                return d
            else:
                if d[selection_field] == selection_string:
                    return d
                else:
                    return None

        tls = self.thread_local_session()
        return [
            d
            for res in tls.query(ServerMonitoring)
            .order_by(desc(ServerMonitoring.sm_int_id))
            .limit(max_reported)
            for d in (row_to_dict(res),)
            if d is not None
        ]

    def server_monitoring_store(
        self,
        message: str = "No message provided",
        what: Optional[str] = None,
        guid: Optional[str] = None,
        content: Dict[str, Any] = {},
    ) -> bool:
        """stores content, a dictionary, into the server monitoring log"""

        if not isinstance(content, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(content))
            )
        tls = self.thread_local_session()
        now = dict(content)
        if what is not None:
            now["content|activity|whatprocess"] = what
        if guid is not None:
            now["content|activity|guid"] = guid
        now["context|info|message"] = message
        current_time = datetime.now()
        now["context|time|time_now"] = str(current_time.isoformat())
        now["context|time|time_boot"] = datetime.fromtimestamp(
            psutil.boot_time()
        ).strftime("%Y-%m-%d %H:%M:%S")

        # should we write this data?  We have the option not to log all messages, to prevent the store getting very full.
        write_content = False
        if self.previous_server_monitoring_time is None:
            write_content = True  # yes if this is the first record written.
        else:
            time_since_last_write = (
                current_time - self.previous_server_monitoring_time
            )  # yes if it's after the server_moni
            t = (
                1000 * float(time_since_last_write.seconds)
                + float(time_since_last_write.microseconds) / 1000
            )
            if t >= self.server_monitoring_min_interval_msec:
                write_content = True

        if write_content:
            json_now = json.dumps(now).encode("utf-8")
            row = ServerMonitoring(
                sm_id=hashlib.sha1(json_now).hexdigest(),
                upload_date=current_time,
                content=json_now,
            )
            tls.add(row)
            self.previous_server_monitoring_time = current_time
            self.previous_server_monitoring_data = now
            tls.commit()
            # finished
            return True
        else:
            return False

    # methods for monitor, which store the contents of an html file
    # in a gridFS store.
    def monitor_store(self, monitoring_id: str, html: str) -> str:
        """stores the monitor output string html.  Overwrites any prior object."""

        if not isinstance(html, str):
            raise TypeError("Can only store string objects, not {0}".format(type(html)))
        tls = self.thread_local_session()
        tls.query(Monitor).filter_by(mo_id=monitoring_id).delete()

        row = Monitor(mo_id=monitoring_id, content=html.encode("utf-8"))
        tls.add(row)
        tls.commit()
        # finished

    def monitor_read(self, monitoring_id: str) -> Optional[str]:
        """loads stored string (e.g. html object) from the monitor collection."""
        tls = self.thread_local_session()
        if res := tls.query(Monitor).filter_by(mo_id=monitoring_id).first():
            return self._to_string(res.content)
        else:
            return None

    # methods for multisequence alignments
    def msa_store(self, msa_token: str, msa: dict) -> Optional[str]:
        """stores the msa object msa under token msa_token."""
        tls = self.thread_local_session()
        if not isinstance(msa, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(msa))
            )

        # we don't replace.  These entries are write once.
        if tls.query(MSA).filter_by(msa_id=msa_token).one_or_none() is None:

            res = MSA(
                msa_id=msa_token,
                upload_date=datetime.now(),
                content=json.dumps(msa).encode("utf-8"),
            )
            tls.add(res)
            tls.commit()

    def msa_read(self, msa_token: str) -> Optional[dict]:
        """loads object from msa collection.
        It is assumed object is a dictionary"""
        tls = self.thread_local_session()
        if res := tls.query(MSA).filter_by(msa_id=msa_token).first():
            return json.loads(res.content)
        else:
            return None

    def msa_delete(self, msa_token: str) -> None:
        """deletes the msa with token msa_token"""
        tls = self.thread_local_session()
        tls.query(MSA).filter_by(msa_id=msa_token).delete()
        tls.commit()
        # finished

    def msa_stored_ids(self) -> List[str]:
        """returns a list of msa tokens of all objects stored"""
        tls = self.thread_local_session()
        return [res.msa_id for res in tls.query(MSA)]

    def msa_delete_unless_whitelisted(self, whitelist: Iterable[str]) -> None:
        """deletes the msa unless the id is in whitelist"""
        tls = self.thread_local_session()
        tls.query(MSA).filter(MSA.msa_id.not_in(whitelist)).delete()
        tls.commit()
        # finished

    # methods for trees

    def tree_store(self, tree_token: str, tree: dict) -> Optional[str]:
        """stores the tree object tree under token tree_token.

        Will not overwrite; requests to do so fail, silently."""

        if not isinstance(tree, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(tree))
            )

        tls = self.thread_local_session()

        if tls.query(TreeStorage).filter_by(ts_id=tree_token).one_or_none() is None:
            row = TreeStorage(
                ts_id=tree_token,
                upload_date=datetime.now(),
                content=json.dumps(tree).encode("utf-8"),
            )
            tls.add(row)
            tls.commit()

    def tree_read(self, tree_token: str) -> Optional[dict]:
        """loads object from tree collection.
        It is assumed object is a dictionary"""
        tls = self.thread_local_session()
        if res := tls.query(TreeStorage).filter_by(ts_id=tree_token).first():
            return json.loads(res.content)
        else:
            return None

    def tree_delete(self, tree_token: str) -> None:
        """deletes the tree with token tree_token"""
        tls = self.thread_local_session()
        tls.query(TreeStorage).filter_by(ts_id=tree_token).delete()
        tls.commit()

    def tree_stored_ids(self) -> List[str]:
        """returns a list of tree tokens of all objects stored"""
        tls = self.thread_local_session()
        return [res.ts_id for res in tls.query(TreeStorage)]

    def tree_delete_unless_whitelisted(self, whitelist: Iterable[str]) -> None:
        """deletes the tree unless the id is in whitelist"""
        tls = self.thread_local_session()
        tls.query(TreeStorage).filter(TreeStorage.ts_id.not_in(whitelist)).delete()
        tls.commit()
        # finished

    # methods for clusters
    def cluster_store(self, clustering_key: str, obj: dict) -> int:
        """stores the clustering object obj.  retains previous version.

        obj: a dictionary to store
        clustering_key: the name of the clustering, e.g. TBSNP12-graph

        Returns:
        current cluster version

        Note: does not replace previous version, but stores a new one.
        Note: to delete legacy versions, call cluster_delete_legacy().
        """
        tls = self.thread_local_session()
        if not isinstance(obj, dict):
            raise TypeError(f"Can only store dictionary objects, not {type(obj)}")

        json_repr = json.dumps(obj, cls=NPEncoder).encode("utf-8")
        cluster = Cluster(
            cluster_build_id=clustering_key,
            upload_date=datetime.now(),
            content=json_repr,
        )
        tls.add(cluster)
        tls.commit()
        # finished
        return cluster.cl_int_id

    def cluster_read(self, clustering_key: str) -> Optional[dict]:
        """loads object from clusters collection corresponding to the most recent version of
        the clustering identified by 'clustering_key'.

        Parameters:
        clustering_key: a string identifying a clustering result

        Returns:
        the clustering information, in the form of a dictionary if it exists, or None if it does not
        """
        tls = self.thread_local_session()
        if (
            res := tls.query(Cluster)
            .filter_by(cluster_build_id=clustering_key)
            .order_by(desc(Cluster.cl_int_id))
            .first()
        ):
            return json.loads(res.content)
        else:
            return None

    def cluster_read_update(
        self, clustering_key: str, current_cluster_version: int
    ) -> Optional[dict]:
        """loads object from clusters collection corresponding to the most recent version
        of the clustering, saved with cluster_build_id = 'clustering_key'.
        it will read only if the current version is different from current_cluster_version; other wise, it returns None


        Parameters:
        clustering_key: a string identifying the cluster
        current_cluster_version: an integer identifying a legacy cluster version

        Returns:
        the clustering information, in the form of a dictionary if it exists, or None if it does not
        """
        tls = self.thread_local_session()
        if (
            res := tls.query(Cluster)
            .filter_by(cluster_build_id=clustering_key)
            .filter(Cluster.cl_int_id != current_cluster_version)
            .order_by(desc(Cluster.cl_int_id))
            .first()
        ):
            return json.loads(res.content)
        return None

    def cluster_latest_version(self, clustering_key: str) -> int:
        """returns id of latest version, which is the maximum number

        Parameters:
        clustering_key: a string identifying the cluster

        Returns:
        cl_int_id, the primary key to the cluster table"""

        tls = self.thread_local_session()
        if (
            res := tls.query(func.max(Cluster.cl_int_id))
            .filter(Cluster.cluster_build_id == clustering_key)
            .first()
        ):
            retVal = res[0]  # it's a tuple
            return retVal
        else:
            return None

    def cluster_keys(self, clustering_name: Optional[str] = None) -> List[str]:
        """lists  clustering keys beginning with clustering_name.  If clustering_name is none, all clustering keys are returned."""

        tls = self.thread_local_session()
        if clustering_name:
            return list(
                sorted(
                    set(
                        res.cluster_build_id
                        for res in tls.query(Cluster).filter(
                            Cluster.cluster_build_id.startswith(clustering_name)
                        )
                    )
                )
            )
        else:
            return list(sorted(set(res.cluster_build_id for res in tls.query(Cluster))))

    def cluster_versions(self, clustering_key: str) -> List[bson.objectid.ObjectId]:
        """lists ids and storage dates corresponding to versions of clustering identifed by clustering_key.
        the newest version is first.
        """
        tls = self.thread_local_session()
        return list(
            tls.query(Cluster)
            .filter_by(cluster_build_id=clustering_key)
            .order_by(desc(Cluster.upload_date))
        )

    def cluster_delete_all(self, clustering_key: str) -> None:
        """delete all clustering objects, including the latest version, stored under clustering_key"""
        tls = self.thread_local_session()
        tls.query(Cluster).filter(Cluster.cluster_build_id == clustering_key).delete()
        tls.commit()
        # finished

    def cluster_delete_legacy_by_key(self, clustering_key: str) -> None:
        """delete all clustering objects, except latest version, stored with key clustering_key"""
        tls = self.thread_local_session()
        cl_int_ids = set()
        for (cl_int_id,) in tls.query(Cluster.cl_int_id).filter_by(
            cluster_build_id=clustering_key
        ):
            cl_int_ids.add(cl_int_id)
        if len(cl_int_ids) == 0:
            return
        else:
            latest_cl_int_id = max(cl_int_ids)
            for this_cl_int_id in cl_int_ids:
                if not this_cl_int_id == latest_cl_int_id:
                    tls.query(Cluster).filter_by(cl_int_id=this_cl_int_id).delete()
        tls.commit()
        # finished

    def cluster_delete_legacy(self, clustering_name: str) -> None:
        """delete all clustering objects, except latest version, stored with  clustering_name"""
        for clustering_key in self.cluster_keys(clustering_name=clustering_name):
            self.cluster_delete_legacy_by_key(clustering_key)

    def refcompressedseq_store(self, guid: str, obj: dict) -> str:
        """stores the json object obj with guid guid.

        Parameters:
        guid:  the sequence identifer
        obj:   a reference compressed sequence representation, as produced by seqComparer.compress().
        Here is an example:

        {
            'A':set([1,2,3]), 'C':set([6]), 'T':set([4]), 'G':set([5]), 'M':{11:'Y', 12:'k'}, 'invalid':0
        }

        If the guid already exists in the database, ignores the request silently."""
        tls = self.thread_local_session()
        if not isinstance(obj, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(obj))
            )

        if "invalid" not in obj.keys():
            raise KeyError(
                "An invalid key must be present.  Keys are: {0}".format(obj.keys())
            )

        # if the record already exists, we don't re-add it
        res = (
            tls.query(RefCompressedSeq.seq_int_id)
            .filter_by(sequence_id=guid)
            .one_or_none()
        )
        if res is None:
            tls.add(
                RefCompressedSeq(
                    sequence_id=guid,
                    invalid=obj["invalid"],
                    examination_date=datetime.now(),
                    content=self.sjc.to_json(obj),
                    prop_actg=None,
                    annotations=json.dumps({}),
                )
            )
        tls.commit()
        # finished

    def refcompressedsequence_read(self, guid: str) -> Any:
        """loads object from refcompressedseq collection.
        It is assumed object stored is a dictionary"""
        tls = self.thread_local_session()
        if (
            rcs := tls.query(RefCompressedSeq.content)
            .filter_by(sequence_id=guid)
            .first()
        ):
            return self.sjc.from_json(rcs.content)
        else:
            return None

    def refcompressedsequence_guids(self) -> Set[str]:
        """loads guids from refcompressedseq collection."""
        tls = self.thread_local_session()
        return set(res.sequence_id for res in tls.query(RefCompressedSeq.sequence_id))

    def guid_annotate(self, guid: str, nameSpace: str, annotDict: dict) -> None:
        """adds multiple annotations of guid from a dictionary;
        all annotations go into a namespace.
        updates the annotation if it exists"""
        tls = self.thread_local_session()
        if not isinstance(annotDict, dict):
            raise TypeError(
                "Can only store dictionary objects, not {0}".format(type(annotDict))
            )

        # The reference compressed sequence must exist.
        rcs = (
            tls.query(RefCompressedSeq)
            .filter(RefCompressedSeq.sequence_id == guid)
            .one_or_none()
        )
        if rcs is None:
            raise RDBMSError(
                "Asked to annotate a record {0} but it does not exist".format(guid)
            )

        if nameSpace == "DNAQuality":
            # coerce examination date to string
            if "examinationDate" in annotDict:
                rcs.examination_date = annotDict["examinationDate"]
                if isinstance(annotDict["examinationDate"], datetime):
                    # convert to isoformat pre-jsonisation
                    annotDict["examinationDate"] = annotDict[
                        "examinationDate"
                    ].isoformat()

            if "propACTG" in annotDict:
                rcs.prop_actg = annotDict["propACTG"]

        annotations = json.loads(rcs.annotations)  # what's there now
        annotations[nameSpace] = annotDict  # replace or add the new namespace
        rcs.annotations = json.dumps(annotations).encode("utf-8")

        tls.commit()
        # finished

    def guids(self) -> Set[str]:
        """returns all registered guids"""
        return self.refcompressedsequence_guids()

    def guids_considered_after(self, addition_datetime: datetime) -> Set[str]:
        """returns all registered guid added after addition_datetime
        addition_datetime: a date of datetime class."""
        tls = self.thread_local_session()
        if not isinstance(addition_datetime, datetime):
            raise TypeError(
                "addition_datetime must be a datetime value.  It is {0}.  Value = {1}".format(
                    type(addition_datetime), addition_datetime
                )
            )
        retVal = []
        for (guid,) in tls.query(RefCompressedSeq.sequence_id).filter(
            RefCompressedSeq.examination_date > addition_datetime
        ):
            retVal.append(guid)
        return set(retVal)

    def _guids_selected_by_validity(self, validity: int) -> Set[str]:
        """returns  registered guids, selected on their validity

        0 = guid is valid
        1 = guid is invalid

        """
        tls = self.thread_local_session()
        return set(
            res.sequence_id
            for res in tls.query(RefCompressedSeq.sequence_id).filter_by(
                invalid=validity
            )
        )

    def singletons(
        self, method: str = "approximate", return_top: int = 1000
    ) -> pd.DataFrame:
        """

        This method is not important in the RDBMS implementation of the fn3persistence store.

        Returns:
        An empty data frame.
        """
        return pd.DataFrame()

    def guids_valid(self) -> set:
        """return all registered valid guids.

        Validity is determined by the contents of the DNAQuality.invalid field, on which there is an index"""
        return self._guids_selected_by_validity(0)

    def guids_invalid(self) -> set:
        """return all invalid guids

        Validity is determined by the contents of the DNAQuality.invalid field, on which there is an index"""
        return self._guids_selected_by_validity(1)

    def guid_exists(self, guid: str) -> bool:
        """checks the presence of a single guid"""
        tls = self.thread_local_session()
        return (
            tls.query(RefCompressedSeq.sequence_id).filter_by(sequence_id=guid).first()
            is not None
        )

    def guid_valid(self, guid: str) -> int:
        """checks the validity of a single guid

        Parameters:
        guid: the sequence identifier

        Returns
        -1   The guid does not exist
        0    The guid exists and the sequence is valid
        1    The guid exists and the sequence is invalid
        """
        tls = self.thread_local_session()
        if (
            res := tls.query(RefCompressedSeq.invalid)
            .filter_by(sequence_id=guid)
            .first()
        ):
            if res.invalid == 0:
                return 0
            elif res.invalid == 1:
                return 1
            else:
                raise ValueError(
                    "invalid is neither 1 nor 0 but {0}".format(res.invalid)
                )
        else:
            return -1

    def guid_examination_time(self, guid: str) -> Optional[datetime]:
        """returns the examination time for a single guid

        Parameters:
        guid: the sequence identifier

        Returns either
        The examination datetime value for this guid OR
        None if the guid does not exist
        """
        tls = self.thread_local_session()
        if (
            res := tls.query(RefCompressedSeq.examination_date)
            .filter_by(sequence_id=guid)
            .first()
        ):
            return res.examination_date
        else:
            return None

    def guids_considered_after_guid(self, guid: str) -> Set[str]:
        """returns all registered guids added after guid
        guid: a sequence identifier"""
        if addition_datetime := self.guid_examination_time(guid):
            return self.guids_considered_after(addition_datetime)
        else:
            raise ValueError("guid is not valid: {0}".format(guid))

    def guid_quality_check(
        self, guid: str, cutoff: Union[float, int]
    ) -> Optional[bool]:
        """Checks whether the quality of one guid exceeds the cutoff.

        If the guid does not exist, returns None.
        If the guid does exist and has quality< cutoff, returns False.
        Otherwise, returns True.
        """
        tls = self.thread_local_session()
        # test input
        if not type(cutoff) in [float, int]:
            raise TypeError(
                "Cutoff should be either floating point or integer, but it is %s"
                % type(cutoff)
            )
        if not type(guid) == str:
            raise TypeError("The guid passed should be as string, not %s" % str(guid))

        # recover record, compare with quality

        res = tls.query(RefCompressedSeq).filter_by(sequence_id=guid).first()
        if res is None:  # no entry for this guid
            return None
        else:
            # report whether it is larger or smaller than cutoff
            return res.prop_actg >= cutoff

    def _guid2seq(self, guidList: Optional[List[str]]) -> Iterable[RefCompressedSeq]:
        """returns the annotations, sequence_id and prop_actg from each RefCompressedSeq for each guid in guidList
        If guidList is None, all items are returned.
        """
        tls = self.thread_local_session()
        if guidList is None:  # rreturn everything
            return tls.query(
                RefCompressedSeq.sequence_id,
                RefCompressedSeq.annotations,
                RefCompressedSeq.prop_actg,
                RefCompressedSeq.examination_date,
            )
        else:
            return tls.query(RefCompressedSeq).filter(
                RefCompressedSeq.sequence_id.in_(guidList)
            )

    def guid2item(
        self, guidList: Optional[List[str]], namespace: str, tag: str
    ) -> dict:
        """returns the annotation (such as sequence quality, which is stored as an annotation)
        in namespace:tag for all guids in guidlist.
        If guidList is None, all items are returned.
        An error is raised if namespace and tag is not present in each record.
        """
        return {
            res.sequence_id: json.loads(res.annotations)[namespace][tag]
            for res in self._guid2seq(guidList)
        }

    def guid2ExaminationDateTime(self, guidList: Optional[List[str]] = None) -> dict:
        """returns quality scores and examinationDate for all guids in guidlist.  If guidList is None, all results are returned."""

        return {
            res.sequence_id: res.examination_date for res in self._guid2seq(guidList)
        }

    def guid2quality(self, guidList: Optional[List[str]] = None) -> Optional[dict]:
        """returns quality scores for all guids in guidlist (or all samples if guidList is None)
        potentially expensive query if guidList is None."""

        return {res.sequence_id: res.prop_actg for res in self._guid2seq(guidList)}

    def guid2propACTG_filtered(self, cutoff: Union[int, float] = 0) -> Dict[str, float]:
        """recover guids which have good quality, > cutoff.
        These are in the majority, so we run a table scan to find these.

        This query is potentially very inefficient- best avoided
        """
        tls = self.thread_local_session()
        query = tls.query(
            RefCompressedSeq.sequence_id, RefCompressedSeq.prop_actg
        ).filter(RefCompressedSeq.prop_actg >= cutoff)

        return {res.sequence_id: res.prop_actg for res in query}

    def guid2items(
        self, guidList: Optional[List[str]], namespaces: Optional[Set[str]]
    ) -> Dict[Any, Dict[str, Any]]:
        """returns all annotations in namespaces, which is a list
        If namespaces is None, all namespaces are returned.
        If guidList is None, all items are returned.
        To do this, a table scan is performed - indices are not used.
        """

        def select_namespaces(annotations: dict) -> dict:
            if namespaces:
                return {ns: annotations[ns] for ns in annotations.keys() & namespaces}
            else:
                return annotations

        return {
            res.sequence_id: select_namespaces(json.loads(res.annotations))
            for res in self._guid2seq(guidList)
        }

    def guid_annotations(self) -> Optional[Dict[Any, Dict[str, Any]]]:
        """return all annotations of all guids"""

        return self.guid2items(None, None)  # no restriction by namespace or by guid.

    def guid_annotation(self, guid: str) -> Optional[Dict[Any, Dict[str, Any]]]:
        """return all annotations of one guid"""

        return self.guid2items([guid], None)  # restriction by guid.

    def guid2neighbour_add_links(
        self,
        guid: str,
        targetguids: Dict[str, Dict[str, int]],
        use_update: bool = False,
    ) -> Dict[str, int]:
        """adds links between guid and their neighbours ('targetguids')

        Parameters:
        guid: the 'source' guid for the matches eg 'guid1'
        targetguids: what is guid linked to, eg
        {
                'guid2':{'dist':12},
                'guid3':{'dist':2}
        }
        use_update -  currently ignored, always False.  Setting True yields NotImplementedError

        This stores links in the guid2neighbour collection.

        Returns:
        The number of records written

        Note:
        uses bulk upload methodology to write fast, as some samples may have thousands or tens of thousands of neighbours

        """

        load_list = []
        for guid2, dist in targetguids.items():
            load_list.append(
                {"sequence_id_1": guid, "sequence_id_2": guid2, "dist": dist["dist"]}
            )
            load_list.append(
                {"sequence_id_1": guid2, "sequence_id_2": guid, "dist": dist["dist"]}
            )
        load_df = pd.DataFrame.from_records(load_list)

        if len(load_df.index) > 0:
            self._bulk_load(load_df, "edge")

    class Guid2NeighbourRepackRet(TypedDict):
        guid: str
        finished: str
        pre_optimisation: dict
        actions_taken: Dict[str, int]

    def guid2neighbour_repack(
        self, guid: str, always_optimise: bool = False, min_optimisable_size: int = 1
    ) -> Guid2NeighbourRepackRet:
        """In the mongo implementation, alters the mongodb representation of the links of guid.

        In the rdbms implementation, this does nothing, and just returns a status report.
        Parameters:
        guid : the identifier of a sample to repack
        always_optimise: consider for optimisation even if there are no 's' (single item) records
        """
        return {
            "guid": guid,
            "finished": datetime.now().isoformat(),
            "pre_optimisation": {
                "s_records": 0,
                "msg": "Repacking not necessary on RDBMS",
            },
            "actions_taken": {},
        }

    class Guid2NeighboursRet(TypedDict):
        guid: str
        neighbours: List[Guid2NeighboursFormats]

    def guid2neighbours(
        self, guid: str, cutoff: int = 20, returned_format: int = 2
    ) -> Guid2NeighboursRet:
        """returns neighbours of guid with cutoff <=cutoff.

        Parameters:
        guid: the sequence identifier
        cutoff: a SNV distance cutoff
        returned_format: see below.

        Returns links either as

        format 1 [[otherGuid, distance],[otherGuid2, distance2],...]
        or as
        format 2 [[otherGuid, distance, N_just1, N_just2, N_either],[],...]
        or as
        format 3 [otherGuid1, otherGuid2, otherGuid3]
        or as
        format 4 [{'guid':otherguid, 'snv':distance}]

            Internally, the documents in guid2neighbour are of the form
            {'guid':'guid1', 'rstat':'s', 'neighbours':{'guid2':{'dist':12, ...}}} OR
            {'guid':'guid1', 'rstat':'m', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}} OR
            {'guid':'guid1', 'rstat':'f', 'neighbours':{'guid2':{'dist':12, ...}, 'guid3':{'dist':5, ...}}

            However, irrespective of their internal representation, this function always returns
            exactly one item for each link of 'guid'; duplicates are not possible.
            The last example occurs when the maximum number of neighbours permitted per record has been reached.
        """

        def f(res):
            otherGuid = res.sequence_id_2
            dist = res.dist
            if returned_format == 1:
                return [otherGuid, dist]
            elif returned_format == 2:
                raise NotImplementedError("format 2 is no longer supported")
            elif returned_format == 3:
                return otherGuid
            elif returned_format == 4:
                return {"guid": otherGuid, "snv": dist}
            else:
                raise ValueError(
                    f"Unable to understand returned_format = {returned_format}"
                )

        tls = self.Session()
        return {
            "guid": guid,
            "neighbours": [
                f(res)
                for res in tls.query(Edge)
                .filter_by(sequence_id_1=guid)
                .filter(Edge.dist <= cutoff)
            ],
        }

    def _set_lock_status(self, lock_int_id, lock_status):
        """locks or unlocks resources identified by lock_int_id, allowing cross- process sequential processing (e.g. insertion)

        To lock, set lock_status =1 ; to unlock, set lock_status =0
        To return the relevant row, set lock_status to None

        See the acquire_lock() method for more details

        returns:
        True if update succeeded, false if it did not

        Technical notes:
        https://docs.sqlalchemy.org/en/14/orm/session_transaction.html
        https://www.amazon.com/Expert-Oracle-Database-Architecture-Thomas-dp-1430262982/dp/1430262982/ref=dp_ob_title_bk

        """
        # make sure there is an entry for this lock
        tls = self.Session()

        lock_row = (
            tls.query(FNLock).filter(FNLock.lock_int_id == lock_int_id).one_or_none()
        )

        # if the row doesn't exist, we add it, with the lock not set.
        if lock_row is None:
            lock_row = FNLock(
                lock_int_id=lock_int_id,
                lock_status=0,
                lock_set_date=datetime.now(),
                uuid=uuid.uuid4().hex,
            )
            tls.add(lock_row)
            tls.commit()

        # analyse the record for this row
        lock_row = (
            tls.query(FNLock)
            .filter(FNLock.lock_int_id == lock_int_id)
            .with_for_update()
            .one()
        )
        if lock_status is None:
            retval = lock_row

        if lock_row.lock_status == 0 and lock_status == 0:
            # it's already unlocked
            retval = True

        elif lock_row.lock_status == 1 and lock_status == 1:
            # it's already locked and we're asked to acquire a lock.  We can't.
            retval = False

        elif lock_row.lock_status == 0 and lock_status == 1:
            # it's already unlocked, we can lock
            lock_row.lock_status = 1
            lock_row.lock_set_date = datetime.now()
            lock_row.uuid = uuid.uuid4().hex
            retval = True

        elif lock_row.lock_status == 1 and lock_status == 0:
            # it's already locked, we can unlock
            lock_row.lock_status = 0
            lock_row.lock_set_date = datetime.now()
            lock_row.uuid = uuid.uuid4().hex
            retval = True

        tls.commit()
        return retval

    def lock_status(self, lock_int_id):
        """determine whether a database-based lock is open (0) or closed (1).

        Parameters:
        lock_int_id: an integer identifier to the lock of interest

        Returns:
        0 if the lock is open
        1 if it is locked"""

        return self._set_lock_status(lock_int_id, None)

    def lock(self, lock_int_id):
        """obtains a database-based lock.

        Parameters:
        lock_int_id: an integer identifier to the lock of interest

        Returns:
        True if the lock is acquired
        False if it is not"""

        return self._set_lock_status(lock_int_id, 1)

    def unlock(self, lock_int_id, force=False):
        """obtains a database-based lock.

        Parameters:
        lock_int_id: an integer identifier to the lock of interest
        force: if True, will unlock irrespective of current status, returning True

        Returns:
        True if the lock is acquired
        False if it is not"""

        res = self._set_lock_status(lock_int_id, 0)
        if force:
            return True
        else:
            return res
