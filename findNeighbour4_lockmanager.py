#!/usr/bin/env python3
""" continuously monitors the findneighbour4 insert lock and resets it if it 
is held for > 90 seconds, since this is indicative of a web worker having crashed.

Background
==========
When running findNeighbour4 with gunicorn, multiple web workers are instantiated.  
Each has a timeout period (90sec at present when fn4 is started using fn4_startup.sh).

If a web worker process, including a process inserting samples, crashes or does not respond, gunicorn will kill it and restart another.
This happens very rarely (~ 1 / 250,000 samples during testing with SARS-CoV-2 genomes).  
However, if the process was inserting before it was killed, this gunicorn initiated process kill can result in

- the insert semaphore remaining set, blocking future inserts
- depending on at what state the insert failed, catwalk and the underlying database may be in an inconsistent state

This process continuously monitors to check whether such an eventuality has occurred, and if so
- checks the database and catwalk to ensure it is in a consistent state
- released the lock

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
import datetime
import sentry_sdk
from version import version
import time
from findn.common_utils import ConfigManager
from findn.cw_seqComparer import cw_seqComparer
from findn.persistence import Persistence

# startup
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_monitor, a service which every 1 minute checks the server insert lock, and releases it if appropriate.""",
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
    parser.add_argument(
        "--max_run_time",
        action="store",
        default=90,
        type=int,
        help="The maximum time a gunicorn web worker can run before it times out.",
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
    logger.info("Monitoring lock")

    logger.info("Connecting to backend data store")
    pm = Persistence()

    last_check_time = datetime.datetime.now()
    logging.info(
        "findNeighbour4_lockmanager startup; will verify insertion & unlock if insertion takes > {0} seconds".format(
            args.max_run_time
        )
    )

    while True:
        PERSIST = pm.get_storage_object(
            dbname=CONFIG["SERVERNAME"],
            connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=0,
            server_monitoring_min_interval_msec=CONFIG[
                "SERVER_MONITORING_MIN_INTERVAL_MSEC"
            ],
            verbose=True,
        )

        # get details of the lock
        lock_details = PERSIST.lock_details(1)
        current_time = datetime.datetime.now()

        if lock_details is not None:
            # a lock is in place; determine whether it has been held for > 90 seconds
            time_difference = current_time - lock_details["lock_set_date"]
            time_difference_seconds = time_difference.total_seconds()
            logging.info(
                "Insert Lock status checked; lock is held for {0} seconds".format(
                    time_difference_seconds
                )
            )
            if time_difference_seconds > args.max_run_time:
                logging.warning(
                    "Releasing lock on {0}".format(lock_details["sequence_id"])
                )
                # release lock; create hybridcomparer object
                hc = cw_seqComparer(
                    reference=CONFIG["reference"],
                    maxNs=CONFIG["MAXN_STORAGE"],
                    snpCeiling=CONFIG["SNPCEILING"],
                    excludePositions=set(CONFIG["excludePositions"]),
                    preComparer_parameters=CONFIG["PRECOMPARER_PARAMETERS"],
                    PERSIST=PERSIST,
                    unittesting=False,
                )

                hc.verify_insertion(lock_details["sequence_id"])
                PERSIST.unlock(1, force=True)
        else:
            logging.info("Insert Lock status checked; lock is not held")
        PERSIST.closedown()

        if args.run_once_only:
            logger.info("Exiting as run_once_only specified")
            exit(0)

        last_check_time = current_time

        time.sleep(60)  # rerun in 1 minute
