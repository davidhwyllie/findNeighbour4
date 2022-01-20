#!/usr/bin/env python
""" extracts insert performance and runs queries against a findneighbour server as part of a benchmark; 

Example usage: 
============== 
# show command line options 
python fn4_benchmark.py --help

# example
pipenv run python3 fn4_benchmark.py ../fn4_atp_test/demos/covid/atp.json /data/data/test/atp

Part of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

"""

# imports
import os
import datetime
import logging
import logging.handlers
import argparse
import pandas as pd
import random
import uuid
import json
from findn.common_utils import ConfigManager
from fn4client import fn4Client
import aiohttp
import asyncio
import time


async def fetch(session, url):
    """[Coroutine to send HTTP request to URLs provided]

    Args:
        session ([aiohttp Client Session object]): [session object]
        url ([string]): [URL to query]

    Returns:
        [list]: [response url, status and time taken]

    https://stackoverflow.com/questions/64564373/python-asyncio-response-time

    """
    tic = time.perf_counter()  # Start timer
    try:
        response = await session.request(method="GET", url=url, timeout=30)
        toc = time.perf_counter()  # Stop timer
        time_taken = toc - tic  # Calculate time taken to get response
        response.raise_for_status()
        key = uuid.uuid4().hex

        results[key] = {"url": url, "status": response.status, "time_taken": time_taken}

        if "neighbours_within" in url:
            txt = await response.text()
            nneighbours = len(json.loads(txt))
        else:
            nneighbours = None
        results[key]["nneighbours"] = nneighbours
        return str(response.url), response.status, time_taken

    except aiohttp.web.HTTPError as http_err:
        logging.error(http_err)
    except Exception as err:
        raise
        logging.error(err)


async def main(urls):

    async with aiohttp.ClientSession() as session:
        htmls = await asyncio.gather(*[fetch(session, url) for url in urls])
        return htmls


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

    # instantiate client
    server_url = "http://localhost:{0}".format(CONFIG["REST_PORT"])
    fn4c = fn4Client(
        server_url
    )  # expects operation on local host; pass baseurl if somewhere else.

    print("Recovering insert rates")
    df = fn4c.server_memory_usage(nrows=int(3e5))
    df = df[df["detail"] == "catwalk comparison engine n_samples"]
    df = df[df["info_message"] == "About to insert"]

    df["delta_usec"] = None
    prior_ix = None
    for i, ix in enumerate(df.index):
        if i > 0:
            delta = datetime.datetime.fromisoformat(
                df.at[ix, "event_time"]
            ) - datetime.datetime.fromisoformat(df.at[prior_ix, "event_time"])
            df.at[prior_ix, "delta_usec"] = delta.microseconds
        prior_ix = ix
    print(df)

    outfile = os.path.join(completedir, "insert_timings.txt")
    df.to_csv(outfile)
    print("Data written to ", outfile)

    print("Recovering guids")
    existing_guids = set(fn4c.guids())
    logger.info("There are {0} samples in the server".format(len(existing_guids)))

    # benchmarking snippet
    one_guid = max(existing_guids)

    random_ordered = list(existing_guids)

    # fire them at the server using asyncio
    print("Starting benchmarking")
    start_time = datetime.datetime.now()
    loop = asyncio.get_event_loop()

    dfs = []

    for j in range(1000):  # replicates
        # constructs 80 requests of four different types, and send them to the server as fast as possible using asyncio
        urls = {
            "guids": [],
            "valid_guids": [],
            "invalid_guids": [],
            "exists": [],
            "neighbours": [],
            "server_time": [],
            "exact_distance": [],
        }

        for i in range(3):
            for endpoint in ["guids", "invalid_guids"]:
                urls[endpoint].append("{0}/api/v2/{1}".format(server_url, endpoint))

        random.shuffle(random_ordered)

        for i, guid in enumerate(random_ordered):
            # add urls; can choose what to use
            this_url = "{0}/api/v2/{1}/exists".format(server_url, guid)
            urls["exists"].append(this_url)
            this_url = "{0}/api/v2/{1}/neighbours_within/4".format(server_url, guid)
            urls["neighbours"].append(this_url)
            urls["server_time"].append(
                "{0}/api/v2/server_time".format(server_url)
            )  # no database access
            this_url = "{0}/api/v2/{1}/{2}/exact_distance".format(
                server_url, one_guid, guid
            )  # 2 database access needed + a

            urls["exact_distance"].append(this_url)
            if i > 1:
                break

        for url_type in ["exists", "neighbours", "server_time", "exact_distance"]:
            print(j, url_type)
            results = {}
            loop.run_until_complete(main(urls[url_type]))
            for key in results.keys():
                results[key]["url_type"] = url_type
                results[key]["replicate"] = j
            dfs.append(pd.DataFrame.from_dict(results, orient="index"))
        df = pd.concat(dfs)
    outfile = os.path.join(completedir, "read_benchmarks.txt")
    df.to_csv(outfile)
    print("Data written to ", outfile)

    loop.close()
