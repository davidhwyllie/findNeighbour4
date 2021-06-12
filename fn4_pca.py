#!/usr/bin/env python
""" 
A component of a findNeighbour4 server which provides relatedness information for bacterial genomes.
It does so using PCA, and supports PCA based cluster generation.

The associated classes compute a variation model for samples in a findNeighbour4 server.
Computation uses data in MongoDb, and is not memory intensive, using configuration information in a 
config file. 

If no config file is provided, it will run in  'testing' mode with the  parameters
in default_test_config.json.  This expects a mongodb database to be running on
the default port on local host. 

Functionality is provided in following classes:
* VariationModel - stores results of producing variant matrix and running PCA
* VariantMatrix - computes sample x variant matrix (requires: PERSIST object for mongodb access; server configuration file) 
* PCARunner - runs PCA on VariantMatrix

These classes are not yet properly unit tested.

Example usage:

pipenv run python3 pca/fn4_pca.py demos/covid/covid_config_v3.json  --outputdir /data/data/pca --analysis_name 2021-04-08  # analyse all, 200 components
pipenv run python3 pca/fn4_pca.py demos/covid/covid_config_v3.json  --outputdir /data/data/pca --analysis_name 2021-04-08 --train_on 50 --n_components 5

see also main()

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.
cronta --
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

"""

# import libraries
import os
import logging
import logging.handlers
import warnings
import datetime
from pathlib import Path
import sentry_sdk
import argparse

# reference based compression storage and clustering modules
from findn.mongoStore import fn3persistence
from findn import DEFAULT_CONFIG_FILE
from findn.common_utils import ConfigManager
from pca.pca import VariantMatrix, PCARunner
from pca.pcadb import PCADatabaseManager, FeatureAssociation

def main():
    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findneighbour4_pca, which extract pcs from a sequence collection
                                     
Example usage: 
============== 
# show command line options 
python fn4_pca.py --help  

# run with debug settings; only do this for unit testing.
python fn4_pca.py     

# run using settings in myConfigFile.json.  
python fn4_pca.py ../config/myConfigFile.json  # all samples, 200 pcs, output to default location
python fn4_pca.py ../config/myConfigFile.json  --outputdir /data/data/pca --analysis_name 2021-04-08 --num_train_on 50 --n_components 5

# covid
pipenv run python3 pca/fn4_pca.py demos/covid/covid_config_v3.json --outputdir /data/data/pca --analysis_name 2021-04-10 --n_components 400
pipenv run python3 pca/fn4_pca.py demos/covid/covid_config_v3.json --outputdir /data/data/pca/realtime_400 --analysis_name 2021-04-10 --n_components 400

""",
    )
    parser.add_argument(
        "path_to_config_file", 
        type=str, 
        action="store", 
        nargs="?", 
        help="the path to the fn4 server configuration file", 
        default=None
    )
    parser.add_argument(
        "connection_config", 
        type=str, 
        action="store", 
        nargs="?", 
        help="the key in the database credentials found in the json file at PCA_CONNECTION_CONFIG_FILE. OR a database configuration string e.g. sqlite:///mydb"
    )
    parser.add_argument(
        "--analysis_name",
        type=str,
        action="store",
        nargs="?",
        help="the name of the analysis.  used as the stem of the .sqlite file into which the output is written.",
        default=datetime.datetime.now().isoformat()[0:10],
    )
    parser.add_argument(
        "--num_train_on",
        type=int,
        action="store",
        nargs="?",
        help="the number of samples to train on.  Omitting this parameter will run an analysis on the whole dataset.",
        default=None,
    )
    parser.add_argument(
        "--n_components",
        type=int,
        action="store",
        nargs="?",
        help="the number of pcs to extract.  If this is smaller than the numbers of samples trained_on, instability is likely.  Default 200",
        default=200,
    )
    parser.add_argument(
        "--cogfile",
        type=str,
        action="store",
        nargs="?",
        help="the name of the analysis.  used as the stem of the .sqlite file into which the output is written.",
        default="/data/data/inputfasta/cog_metadata.csv"
    )
    parser.add_argument(
        "--remove_existing_data",  
        action="store_true", 
        default=False, 
        help="if this option is specified, it will delete any existing data from the database. "
    )
    args = parser.parse_args()
    
    ############################ LOAD CONFIG ######################################
    logging.info("findNeighbour4 PCA modelling .. reading configuration file.")

    config_file = args.path_to_config_file
    if config_file is None:
        config_file = DEFAULT_CONFIG_FILE
        warnings.warn(
            f"No config file name supplied; using configuration in {DEFAULT_CONFIG_FILE}, suitable only for testing, not for production."
        )

    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()

    # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
        CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
        logging.info("Set mongodb connection string  from environment variable")
    else:
        logging.info("Using mongodb connection string from configuration file.")

    # determine whether a FN_SENTRY_URLenvironment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        logging.info("Set Sentry connection string from environment variable")
    else:
        logging.info("Using Sentry connection string from configuration file.")

    ########################### SET UP LOGGING #####################################
    # create a log file if it does not exist.
    logging.info("Starting logging")
    logdir = os.path.dirname(CONFIG["LOGFILE"])
    Path(os.path.dirname(CONFIG["LOGFILE"])).mkdir(parents=True, exist_ok=True)

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
    logfile = os.path.join(logdir, "pca-{0}".format(os.path.basename(CONFIG["LOGFILE"])))
    logging.info("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode="a", maxBytes=1e7, backupCount=7)

    formatter = logging.Formatter("%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(logging.StreamHandler())

    ######################### launch sentry if API key provided ##################################
    # determine whether a FN_SENTRY_URLenvironment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        logging.info("Set Sentry connection string from environment variable")
    else:
        logging.info("Using Sentry connection string from configuration file.")

    if "SENTRY_URL" in CONFIG.keys():
        logger.info("Launching logger")
        sentry_sdk.init(CONFIG["SENTRY_URL"])

    # prepare to connection
    logging.info("Connecting to backend data store")
    try:
        PERSIST = fn3persistence(
            dbname=CONFIG["SERVERNAME"],
            connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=CONFIG["DEBUGMODE"],
            server_monitoring_min_interval_msec=0,
        )
    except Exception:
        logging.info("Error raised on creating persistence object")
        raise

    if args.remove_existing_data:
        logging.warning("fn4_pca is set to remove all existing data from database")
    else:
        logging.warning("fn4_pca is set to add build to existing data")

    # instantiate builder for PCA object
    rebuild = True
    logging.info("Rebuild is {0}".format(rebuild))
    if rebuild:
        try:
            var_matrix = VariantMatrix(CONFIG, PERSIST)
        except Exception:
            logging.info("Error raised on instantiating Variant Matrix object")
            raise

        logging.info("Building snp matrix")
        var_matrix.build(num_train_on=args.num_train_on)
        logging.info("Running PCA on snp matrix")
        pca_runner = PCARunner(var_matrix)
        pca_runner.run(n_components=args.n_components, pca_parameters={})
        vm = pca_runner.cluster()

        logging.info("Exporting to database")
        pdm = PCADatabaseManager(
            connection_config=args.connection_config, debug=args.remove_existing_data)
        logger.info("Loading cog-uk metadata")
        pdm.store_cog_metadata(cogfile=args.cogfile)
        logger.info("Storing variation model and PCA")
        pdm.store_variation_model(vm)
        pdm.make_contingency_tables(only_pc_cats_less_than_days_old = 120)
        

# startup
if __name__ == "__main__":
    main()
