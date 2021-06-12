""" dumps the contents of the findneighbour4 server's database into flatfiles compatible with loading into an RDBMS

It does this in an atomic manner- it is guaranteed that, even if the server is in use, only the samples (guids) 
present when the routine starts will be used.

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

# import libraries
import os
import pathlib
import logging
import sqlalchemy
import logging.handlers
import datetime
import pandas as pd
import sentry_sdk
import argparse
import progressbar
import subprocess
import json
from io import StringIO

# config
from findn.common_utils import ConfigManager

# data
from findn.mongoStore import fn3persistence
from pca.pcadb import PCADatabaseManager

if __name__ == "__main__":

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_dumpneighbours, a findNeighbour4 component.
                                     
Example usage: 
============== 

## does not require findNeighbour4_server to be running
Example dumping all data to output
python findNeighbour4_dumpneighbours.py demos/covid/covid_config_v3.json  /data/output/ 

Test usage
# no upload of data
python findNeighbour4_dumpneighbours.py demos/covid/covid_config_v3.json  /tmp/test/ --debug
# upload to the database identified by 'unittest_oracle'
pipenv run python findNeighbour4_dumpneighbours.py demos/covid/covid_config_v3.json  testdata/fn4snp/dumpneighbours --debug --connection_config unittest_oracle

if a config file is not provided, it will terminate

Checks for new sequences are conducted once per minute.

""",
    )
    parser.add_argument(
        "path_to_config_file", type=str, action="store", nargs="?", help="the path to the configuration file", default=""
    )
    parser.add_argument(
        "outputdir",
        type=str,
        action="store",
        nargs="?",
        help="the output directory.  Will try to make the directory if it does not exist",
    )
    
    parser.add_argument(
        "upload_",
        type=str,
        action="store",
        nargs="?",
        help="the output directory.  Will try to make the directory if it does not exist",
    )

    parser.add_argument(
        "--connection_config", 
        type=str, 
        action="store", 
        nargs="?", 
        help="If supplied, the key in the database credentials found in the json file at PCA_CONNECTION_CONFIG_FILE. OR a database configuration string e.g. sqlite:///mydb.  Will be used to upload the data obtained."
    )
    parser.add_argument("--debug", default=False, help="only analyse the first 2000 samples", action="store_true")
    args = parser.parse_args()

    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 guidetree .. reading configuration file.")

    if len(args.path_to_config_file) > 0:
        config_file = args.path_to_config_file
        debugmode = False
        logging.info(config_file)
    else:
        raise FileNotFoundError("No config file name supplied ; this is required ")

    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()

    ########################### SET UP LOGGING & BUG TRACKING #################################
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
    logfile = os.path.join(logdir, "guidetree_{0}".format(os.path.basename(CONFIG["LOGFILE"])))
    print("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode="a", maxBytes=1e7, backupCount=7)
    formatter = logging.Formatter("%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # launch sentry if API key provided
    if "SENTRY_URL" in CONFIG.keys():
        logger.info("Launching communication with Sentry bug-tracking service")
        sentry_sdk.init(CONFIG["SENTRY_URL"])

    ##################################################################################################
    # open PERSIST 
    logger.info("Connecting to backend data store")
    try:
        PERSIST = fn3persistence(
            dbname=CONFIG["SERVERNAME"], connString=CONFIG["FNPERSISTENCE_CONNSTRING"], debug=0
        )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception:
        logger.exception("Error raised on creating persistence object")
        raise

    logging.info("Loading all samples & annotations.")
    sample_annotations = PERSIST.guid2items(None, None)
    all_samples = list(sample_annotations.keys())
    print("There are {0} samples.".format(len(all_samples)))
    
    #
    ########################### create output directory if it does not exist ##########
    if args.outputdir is not None:
        outputdir = args.outputdir
        logging.info("Writing output to {0}".format(outputdir))
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    ################################### DEBUG : ONLY CONSIDER 200 for testing #####################################
    export_batch_size = 100000
    logging.info("Debug mode is {0}".format(args.debug))
    
    if args.debug:
        max_samples = 25
        export_batch_size = 4
        if len(all_samples)> max_samples:
            select_keys = all_samples[0:max_samples]
            new_subset =  { x: sample_annotations[x] for x in  select_keys}
            sample_annotations = new_subset
            all_samples = list(sample_annotations.keys())
            logging.warning("Running in debug mode.  Only {0} samples will be considered".format(max_samples))
    else:
        logging.info("Running in normal mode.  Will add data to empty fn4_sample & neighbour tables, but will not delete data.")
    #################################  prepare to iterate #########################################################
    # write the annotations
    outputfile = os.path.join(outputdir, 'samples.json')
    sample_annotation_df = pd.DataFrame.from_dict(sample_annotations, orient='index')
    sample_annotation_df['invalid'] = sample_annotation_df['DNAQuality:invalid']
    sample_annotation_df = sample_annotation_df['invalid']
    with open(outputfile,'wt') as f:
        json.dump(sample_annotation_df.to_dict(),f)

    all_samples = set(all_samples)
    bar = progressbar.ProgressBar(max_value=len(all_samples))

    # iterate over all samples
    export_batch = 0
    to_export = []
    for i, this_sample in enumerate(all_samples):
        test_sample_neighbour_info = PERSIST.guid2neighbours(this_sample, 1e6, returned_format=1)  

        # filter these neighbours.   Only include neighbours which are part of all_samples.
        # because samples get added to the server all the time, this is not guaranteed.
        test_sample_neighbours = []
        for item in test_sample_neighbour_info["neighbours"]:
            if item[0] in all_samples:
                to_export.append([this_sample, item[0], item[1]])
        if i % 500 == 0:
            bar.update(i)

        if i % export_batch_size== 0:           # export a batch of records
            export_batch +=1
            outputfile = os.path.join(outputdir, 'neighbours_{0}.json'.format(export_batch))
            with open(outputfile,'wt') as f:
                json.dump(to_export,f)
            to_export = []

    bar.finish() 

    outputfile = os.path.join(outputdir, 'neighbours_final.json')
    with open(outputfile,'wt') as f:
        json.dump(to_export,f) 
    logging.info("Dump finished.  Output is at {0}".format(outputdir))

    #--------------------------------------------------------------------------------------------------------------------------------

    if args.connection_config is not None:
        logging.info("Connection config {0} provided to upload data".format(args.connection_config))
        pdm = PCADatabaseManager(connection_config=args.connection_config, debug=args.debug)      # if not in debug mode, won't add data if exists
        pdm.fn4_bulk_upload(outputdir)
            

    else:
        logging.info("No connection_config specified; not uploading to database.")
logging.info("Complete")

