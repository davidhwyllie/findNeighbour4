#!/usr/bin/env python3
""" produces static depictions of server activity, based on findNeighbour monitoring data """

# import libraries
import os
import sys
import logging
import warnings
import pandas as pd
import pathlib
import sentry_sdk
import json

# fn3 storage module
from mongoStore import fn3persistence


if __name__ == "__main__":
    # command line usage.  Pass the location of a config file as a single argument.
    # an example config file is default_test_config.json

    ############################ LOAD CONFIG ######################################
    print("findNeighbour4_static_performance_report .. reading configuration file.")

    max_batch_size = 100

    if len(sys.argv) == 2:
        configFile = sys.argv[1]
    else:
        configFile = os.path.join("..", "config", "default_test_config.json")
        warnings.warn(
            "No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. "
        )

    # open the config file
    try:
        with open(configFile, "r") as f:
            CONFIG = f.read()

    except FileNotFoundError:
        raise FileNotFoundError(
            "Passed one parameter, which should be a CONFIG file name; tried to open a config file at {0} but it does not exist ".format(
                sys.argv[1]
            )
        )

    if isinstance(CONFIG, str):
        CONFIG = json.loads(CONFIG)  # assume JSON string; convert.

    # check CONFIG is a dictionary
    if not isinstance(CONFIG, dict):
        raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))

    # check that the keys of config are as expected.

    required_keys = set(["IP", "REST_PORT", "DEBUGMODE", "LOGFILE", "MAXN_PROP_DEFAULT"])
    missing = required_keys - set(CONFIG.keys())
    if not missing == set([]):
        raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))

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
    file_handler = logging.FileHandler("perfreport-{0}".format(os.path.basename(CONFIG["LOGFILE"])))
    formatter = logging.Formatter("%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(logging.StreamHandler())

    # launch sentry if API key provided
    if "SENTRY_URL" in CONFIG.keys():
        logger.info("Launching logger")
        sentry_sdk.init(CONFIG["SENTRY_URL"])

    # set min logging interval if not supplied
    if "SERVER_MONITORING_MIN_INTERVAL_MSEC" not in CONFIG.keys():
        CONFIG["SERVER_MONITORING_MIN_INTERVAL_MSEC"] = 0

    # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuration file.
    if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
        CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
        print("Set connection string to mongodb from environment variable")
    else:
        print("Using connection string to mongodb from configuration file.")

    # determine whether a FN_SENTRY_URL environment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        print("Set Sentry connection string from environment variable")
    else:
        print("Using Sentry connection string from configuration file.")

    ########################  START Operations ###################################
    logger.info("Preparing to produce visualisations")

    print("Connecting to backend data store at {0}".format(CONFIG["SERVERNAME"]))
    try:
        PERSIST = fn3persistence(
            dbname=CONFIG["SERVERNAME"],
            connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=CONFIG["DEBUGMODE"],
            server_monitoring_min_interval_msec=CONFIG["SERVER_MONITORING_MIN_INTERVAL_MSEC"],
        )
    except Exception as e:
        logger.exception("Error raised on creating persistence object")
        if e.__module__ == "pymongo.errors":
            logger.info("Error raised pertains to pyMongo connectivity")
            raise

    print("Recovering data from server ..")
    insert_data = PERSIST.recent_server_monitoring(
        selection_field="context|info|message", selection_string="About to insert", max_reported=40000
    )
    res = pd.DataFrame.from_records(insert_data, index="_id")
    res.to_excel("inserts.xlsx")
    print("Finished")
    exit(0)
