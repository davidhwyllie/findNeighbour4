#!/usr/bin/env python
""" loads samples into a findNeighbour4 server
from a multi-fasta file containing a simulated phylogeny.

Example usage: 
============== 
# show command line options 
python covid_querytest.py --help  

# example usage
# set up server if not already runing
./fn4_startup.sh demos/covidsim/config_performancetest.json

# load data
pipenv run python3 demo/covidsim.py demos/covidsim/config_performancetest.json /data/data/pca/sim/fasta  sim1.fasta 

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

# imports
import os
import glob
import datetime
import Bio
import logging
import logging.handlers
import argparse
import warnings
import pandas as pd
import sentry_sdk
from collections import Counter
import uuid
from findn.common_utils import ConfigManager
from fn4client import fn4Client
import aiohttp
import asyncio

http_results = {}

async def fetch(session, url):
    async with session.get(url) as response:
        return await response.text()

async def main(urls):
    async with aiohttp.ClientSession() as session:
        for url in urls:
            http_response = await fetch(session, url)
            token = uuid.uuid4().hex
            http_results[token] = dict(
                url=url, response=http_response, response_time=datetime.datetime.now()
            )


## define functions and classes
class DatabaseMonitorInoperativeError(Exception):
    """insert failed"""

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


if __name__ == "__main__":

    # command line usage.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Queries findneighbour server
                                     

                            Example usage: 
                            ============== 
                            # show command line options 
                            pipenv run python3 demo/covid_querytest.py --help  

                            # run tests
                            pipenv run python3 demo/covid_querytest.py demos/covidsim/config.json /data/data/test 
                            pipenv run python3 demo/covid_querytest.py  demos/covid/covid_config_v3.json /data/data/test
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
        "outputdir",
        type=str,
        action="store",
        nargs="?",
        help="the directory in which output files will appear",
        default="",
    )

    args = parser.parse_args()

    ####################################    STARTUP ###############################
    # validate input
    if not os.path.exists(args.outputdir):
        # that's an error
        raise FileNotFoundError(
            "The directory specified for output files does not exist: '{0}'".format(
                args.outputdir
            )
        )
    else:
        if not os.path.isdir(args.outputdir):
            raise FileNotFoundError(
                "The path specified for output files is not a directory: '{0}'".format(
                    args.outputdir
                )
            )

    config_file = args.path_to_config_file
    if config_file is None:
        exit("No config file name supplied.")

    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()

    # the fasta dir exists.  Make sure we have log directories.
    logdir = os.path.join(args.outputdir, "logs")
    completedir = os.path.join(args.outputdir, "completed")
    for testdir in [logdir, completedir]:
        os.makedirs(testdir, exist_ok=True)

    # launch logger
    logger = logging.Logger("sim_load_multifasta")
    logger.setLevel(logging.INFO)
    timenow = datetime.datetime.now().isoformat()
    logfile = os.path.join(logdir, "fn4_load_multifasta_v2.log")
    file_handler = logging.handlers.RotatingFileHandler(
        logfile, mode="a", maxBytes=1e7, backupCount=7
    )
    formatter = logging.Formatter(
        "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s "
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info("Startup with arguments: {0}".format(args))

    print("To see what is happening, do watch tail {0}".format(logfile))
    # launch comms with Sentry
    if os.environ.get("FN_SENTRY_URL") is not None:
        logger.info("Launching communication with Sentry bug-tracking service")
        sentry_sdk.init(os.environ.get("FN_SENTRY_URL"))
        logger.info("Sentry comms established")
    else:
        logger.info(
            "No error monitoring via sentry.  Set environment variable FN_SENTRY_URL to do so."
        )

    # instantiate client
    server_url = "http://localhost:{0}".format(CONFIG["REST_PORT"])
    fn4c = fn4Client(
        server_url
    )  # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    logger.info("There are {0} samples in the server".format(len(existing_guids)))

   
    # benchmarking snippet
    one_guid = min(existing_guids)
    urls = []
    for i, guid in enumerate(fn4c.guids()):
        # add urls; can choose what to use
        #urls.append("{0}/api/v2/{1}/exists".format(server_url, guid)) # one database access
        
        urls.append("{0}/api/v2/server_time".format(server_url))  # no database access
        #urls.append("{0}/api/v2/{1}/{2}/exact_distance".format(server_url, one_guid, guid))     # 2 database access needed + a 
        if i > 1000:
            break

    # fire them at the server using asyncio
    start_time = datetime.datetime.now()
    loop = asyncio.get_event_loop()
    loop.run_until_complete(main(urls))

    res = pd.DataFrame.from_dict(http_results, orient="index")
    res["start_time"] = start_time
    res["delta"] = res["response_time"] - res["start_time"]
    res["msec"] = [1000 * x.total_seconds() for x in res["delta"]]
    print(res)
    logging.info("Finished, terminating program.")
