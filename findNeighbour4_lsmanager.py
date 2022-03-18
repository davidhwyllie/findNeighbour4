#!/usr/bin/env python3
""" backs up sequence data to a local tar file, which allows rapid restarting.

Background
==========

findneighbour4 by default keeps a copy of the reference compressed sequence data on the server on which findneighbour4 runs.

Why does it do this?
--------------------
When large numbers (> 100k samples) are present in the server, if the server has to load the reference compressed sequence 
data from the database then database access can slow operations, delay operationsa and consume large amounts of network bandwidth.

It only does this under two circumstances:
* When it is restarting: if catwalk (a relatedness system used by findneighbour4) is not running, then it will populate catwalk from sequences in the database
* When performing PCA on sequence data

Sequence data, once added to findneighbour4, is immutable.  Therefore, findneighbour4 can keep a local copy of the data in the database.  By default it does, keeps it up to date, and uses it when performing restarts and PCA.

The storage location is determined from the specified LOGFILE in the config file supplied on startup.
```
./fn4_startup.sh config.json
```

findneighbour4 uses the directory in which LOGFILE is specified to create a directory called 
```localcache/{SERVERNAME}```.  So if the config.json file looks like this
```
{
"DESCRIPTION":"A test server operating in on localhost for unit testing using mapped TB data",
"IP":"127.0.0.1",
"INPUTREF":"reference/TB-ref.fasta",
"EXCLUDEFILE":"reference/TB-exclude-adaptive.txt",
"DEBUGMODE":2,
"SERVERNAME":"fn4server1",      
..
"LOGFILE":"/data/fn4storage/fn4/fn4_server1.log",
...
}
```
then it will create a directory ```/data/fn4storage/fn4/localcache/fn4server1```

Inside this will be two subdirectories
```rcs```
```mfa```

This directory can contain large amounts of data, so you should ensure there is sufficient space.  It does not need to be backed up, because it will be automatically regenerated.


A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

# import libraries
import os
import argparse
import logging
import logging.handlers
import warnings
import pathlib
import sentry_sdk
from version import version
import time
from findn.common_utils import ConfigManager
from findn.cw_seqComparer import cw_seqComparer
from findn.persistence import Persistence
from localstore.localstoreutils import LocalStore

# startup
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_lsmanager, a service which backs up sequence data into a local tar file ('localstore').""",
    )
    parser.add_argument(
        "--path_to_config_file",
        type=str,
        action="store",
        nargs="?",
        help="the path to the configuration file",
        default=None,
    )
    parser.add_argument(
        "--run_once_only",
        action="store_true",
        help="run once; do not keep running.  Default False.  Mainly useful for unit testing.",
    )
    args = parser.parse_args()

    ############################ LOAD CONFIG ######################################

    config_file = args.path_to_config_file
    if config_file is None:
        config_file = os.path.join("config", "default_test_config.json")
        warnings.warn(
            "No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. "
        )
    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config(not_debug_mode=True)  # override debug settings

    ########################### SET UP LOGGING #####################################
    # create a log file if it does not exist.
    print("Starting logging")
    logdir = os.path.dirname(CONFIG["LOGFILE"])
    pathlib.Path(os.path.dirname(CONFIG["LOGFILE"])).mkdir(parents=True, exist_ok=True)

    # set up logger
    loglevel = logging.INFO
    if "LOGLEVEL" in CONFIG.keys():
        if CONFIG["LOGLEVEL"] == "WARN":
            loglevel = logging.WARN
        elif CONFIG["LOGLEVEL"] == "DEBUG":
            loglevel = logging.DEBUG

    # configure logging object
    logger = logging.getLogger()
    logger.setLevel(loglevel)
    logfile = os.path.join(
        logdir, "monitor-{0}".format(os.path.basename(CONFIG["LOGFILE"]))
    )
    print("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(
        logfile, mode="a", maxBytes=1e7, backupCount=7
    )
    formatter = logging.Formatter(
        "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s "
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(logging.StreamHandler())

    # launch sentry if API key provided
    if "SENTRY_URL" in CONFIG.keys():
        logger.info("Launching sentry client")
        sentry_sdk.init(CONFIG["SENTRY_URL"], release=version)

    # #######################  START Operations ###################################

    while True:

        tarfilename = os.path.join(cfm.rcscache, "rcs.tar")
        logging.info("Creating localstore in {0}".format(tarfilename))
        # create database connection
        pm = Persistence()
        PERSIST = pm.get_storage_object(
            dbname=CONFIG["SERVERNAME"],
            connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=0,
            server_monitoring_min_interval_msec=CONFIG[
                "SERVER_MONITORING_MIN_INTERVAL_MSEC"
            ],
            verbose=True,
        )

        # create tar file connection
        localstore = LocalStore(
            tarfilename
        )

        # create catwalk connection
        hc = cw_seqComparer(
            reference=CONFIG["reference"],
            maxNs=CONFIG["MAXN_STORAGE"],
            snpCeiling=CONFIG["SNPCEILING"],
            excludePositions=set(CONFIG["excludePositions"]),
            preComparer_parameters=CONFIG["PRECOMPARER_PARAMETERS"],
            PERSIST=PERSIST,
            unittesting=False,
            disable_insertion=True,
            localstore= localstore,
            update_localstore_only = True  
        )

        if args.run_once_only:
            logger.info("Exiting as run_once_only specified")
            exit(0)

        time.sleep(3600)  # rerun in 1 hour
