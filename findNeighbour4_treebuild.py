""" tree building for findNeighbour4
assumes a findNeighbour4 server is running, with the connection string stated in demos/AC587/config/config_cl.json.

./fn4_startup.sh demos/AC587/config/config_cl.json

This code performs clustering.

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

# import libraries
import os
import warnings
import pandas as pd
import pathlib
import sentry_sdk
import argparse
import progressbar
import time
from version import version

# logging
import logging
import logging.handlers

# startup
from findn.persistence import Persistence
from findn.msa import MSAStore
from findn.cw_seqComparer import cw_seqComparer
from findn.common_utils import ConfigManager


if __name__ == "__main__":

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_treebuild, a findNeighbour4 component.
                                     
Example usage: 
============== 

## does not require findNeighbour4_server to be running
Minimal example:
python findNeighbour4_treebuild.py config/myconfig_file.json  

if a config file is not provided, it will run once is debug mode: it will run once, and then terminate.  This is useful for unit testing.  If a config file is specified, the clustering will  run until terminated.  

Checks for new sequences are conducted once per minute.


""",
    )
    parser.add_argument(
        "path_to_config_file",
        type=str,
        action="store",
        nargs="?",
        help="the path to the configuration file",
        default="",
    )
    parser.add_argument(
        "--run_once",
        help="run once only.  The program will stop after the initial run if this option is used, and will need to be restarted without this option to function normally",
        action="store_true",
    )
    args = parser.parse_args()


    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 tree build .. reading configuration file.")

    if len(args.path_to_config_file) > 0:
        config_file = args.path_to_config_file
        debugmode = False
    else:
        config_file = os.path.join("config", "default_test_config.json")
        debugmode = True
        warnings.warn(
            "No config file name supplied ; using a configuration ('default_test_config.json') suitable only for testing, not for production. "
        )
    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config(
        not_debug_mode=True
    )  # never start in debug mode, which wipes all data; as we need data in the server to do any tests

    ########################### SET UP LOGGING #####################################
    # create a log file if it does not exist.
    print("Starting logging")
    logdir = os.path.dirname(CONFIG["LOGFILE"])
    pathlib.Path(os.path.dirname(CONFIG["LOGFILE"])).mkdir(parents=True, exist_ok=True)

    # set up logger
    logger = logging.getLogger()
    loglevel = logging.INFO
    if "LOGLEVEL" in CONFIG.keys():
        if CONFIG["LOGLEVEL"] == "WARN":
            loglevel = logging.WARN
        elif CONFIG["LOGLEVEL"] == "DEBUG":
            loglevel = logging.DEBUG

    # configure logging object
    logger.setLevel(loglevel)
    logfile = os.path.join(
        logdir, "clustering-{0}".format(os.path.basename(CONFIG["LOGFILE"]))
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

    # determine whether a FN_SENTRY_URLenvironment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        print("Set Sentry connection string from environment variable")
    else:
        print("Using Sentry connection string from configuration file.")

    # launch sentry if API key provided
    if "SENTRY_URL" in CONFIG.keys():
        logger.info("Launching communication with Sentry bug-tracking service")
        sentry_sdk.init(CONFIG["SENTRY_URL"], release=version)

  
    ########################### prepare to start tree building ####################################
    # construct the required global variables

    logger.info("Connecting to backend data store")
    pm = Persistence()
    PERSIST = pm.get_storage_object(
        dbname=CONFIG["SERVERNAME"],
        connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
        debug=0,
        verbose=True,
    )

    ##################### open a catwalk connection.  this is used to compute expected Ns ####
    hc = cw_seqComparer(
        reference=CONFIG["reference"],
        maxNs=CONFIG["MAXN_STORAGE"],
        snpCeiling=CONFIG["SNPCEILING"],
        excludePositions=CONFIG["excludePositions"],
        preComparer_parameters=CONFIG["PRECOMPARER_PARAMETERS"],
        PERSIST=PERSIST,
        disable_insertion=True,
    )

    ################################# tree building #############################################
    # open PERSIST and cw_seqComparer object used by all samples
    # this is only used for data access and msa.
    
    # get a clustering object's settings
    logger.info("Iterating over stored multiple sequence alignments ...")

    ms = MSAStore(
        PERSIST=PERSIST, in_ram_persistence_time=60
    )  # persist result for up to 60 seconds

    # recover existing msas
    stored_msa = ms.existing_tokens()

    bar = progressbar.ProgressBar(max_value=len(stored_msa))

    for i,token in enumerate(stored_msa):
        print(token)

        msa_object = ms.load(token)
        print(type(msa_object))     # MSAResult

        ## various access methods exist
        print(msa_object.msa_df())      # pandas
        print(msa_object.fconst)        # constant sites (required for iqtree)
        #print(msa_object.msa_fasta())   # fasta output
        #print(msa_object.msa_dict())    # dictionary

        #print(msa_object.msa_html())    # readable
        # load multisequence alignment

        #bar.update(i + 1)
        # token identifies cluster contents, nature of analysis, and outgroup
        
    bar.finish()

   