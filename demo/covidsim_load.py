#!/usr/bin/env python
""" loads samples into a findNeighbour4 server
from a multi-fasta file containing a simulated phylogeny.

pipenv run python3 findNeighbour4_server.py demos/covidsim/config.json

Example usage: 
============== 
# show command line options 
python covidsim.py --help  

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
        description="""Loads simulated data into findneighbour server
                                     

                            Example usage: 
                            ============== 
                            # show command line options 
                            python covidsim.py --help  

                            # load into server 
                            pipenv run python3 demo/covidsim.py demos/covidsim/config.json /data/data/pca/sim/fasta  sim1.fasta 
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
        "fastadir",
        type=str,
        action="store",
        nargs="?",
        help="the directory in which fasta files will appear",
        default="",
    )
    parser.add_argument(
        "fileglob",
        type=str,
        action="store",
        nargs="?",
        help="a pattern to glob for",
        default="",
    )

    args = parser.parse_args()

    ####################################    STARTUP ###############################
    # validate input
    if not os.path.exists(args.fastadir):
        # that's an error
        raise FileNotFoundError(
            "The directory specified for fasta files does not exist: '{0}'".format(
                args.fastadir
            )
        )
    else:
        if not os.path.isdir(args.fastadir):
            raise FileNotFoundError(
                "The path specified for fasta files is not a directory: '{0}'".format(
                    args.fastadir
                )
            )

    config_file = args.path_to_config_file
    if config_file is None:
        exit("No config file name supplied.")

    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()

    # the fasta dir exists.  Make sure we have log directories.
    logdir = os.path.join(args.fastadir, "logs")
    completedir = os.path.join(args.fastadir, "completed")
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
    clustering_created = False
    logger.info("There are {0} existing guids".format(len(existing_guids)))

    fastafiles = sorted(glob.glob(os.path.join(args.fastadir, args.fileglob)))
    logger.info(
        "There are {0} fasta files waiting to be loaded".format(len(fastafiles))
    )

    print("Found the following fasta files:")
    print(fastafiles)
    print("Will analyse sequentially.")

    ## optionally raise an error if this list is too long
    for fastafile in fastafiles:
        existing_guids = fn4c.guids()
        print(
            "At startup, there are {0} samples loaded.  Note, fn4_startup.sh will prepopulate the catwalk component from the database.".format(
                len(existing_guids)
            )
        )

        if len(existing_guids) > 0:
            warnings.warn(
                "Note: guids remain on startup: {0}.  This may be normal if you have run this script before.".format(
                    len(existing_guids)
                )
            )

        # add the reference sequence as the root if not already present
        ref_guid = "--Reference--"
        ref_guid_present = fn4c.guid_exists(ref_guid)
        for record in Bio.SeqIO.parse("reference/nc_045512.fasta", "fasta"):
            refseq = str(record.seq).upper()
            reflen = len(refseq)

        if not ref_guid_present:
            logger.info("Adding reference")
            res = fn4c.insert(guid=ref_guid, seq=refseq)
        else:
            logger.info("Reference already present")

        logger.info("Scanning {0}".format(fastafile))
        fastastem = os.path.basename(fastafile)
        fastastem = fastastem.replace(".fasta", "")
        print("Analysing:", fastastem)

        nSkipped = 0
        nBad = 0
        nGood = 0
        failed = []
        i = 0
        records = []
        seqdict = {}
        for record in Bio.SeqIO.parse(fastafile, "fasta"):

            i = i + 1
            t1 = datetime.datetime.now()

            sample_id, dmy_date, is_variant = record.id.split(":")
            d_date, m_date, y_date = dmy_date.split("-")
            sample_date = datetime.date(int(y_date), int(m_date), int(d_date))
            # print(sample_id, dmy_date, is_variant)

            seq = str(record.seq).upper()
            seq = seq.replace("?", "N")
            seq = seq.replace(" ", "N")
            counter = Counter(list(seq))
            seqdict[sample_id] = seq
            res = {
                "record_id": record.id,
                "sample_date": sample_date,
                "guid": sample_id,
                "seqlen": len(seq),
                **counter,
                "is_variant": is_variant,
            }
            if sample_id not in existing_guids:
                records.append(res)

        logger.info(
            "Complete read of simulated data not currently in server into RAM.  There are {0} sequences.  Sorting ...".format(
                len(records)
            )
        )

        seq_df = pd.DataFrame.from_records(records)
        if len(seq_df.index) > 0:
            seq_df.sort_values(by="sample_date", ascending=True, inplace=True)

            logging.info("Loading data into server")

            addition_dates = sorted(seq_df["sample_date"].unique())
            n_added = 0
            for addition_date in addition_dates:
                print(addition_date)
                for ix in seq_df.index[seq_df["sample_date"] == addition_date]:
                    sample_id = seq_df.at[ix, "guid"]
                    seq = seqdict[sample_id]
                    n_added += 1
                    print("Insertion", n_added, ix, len(seq))
                    res = fn4c.insert(guid=sample_id, seq=seq)

    # finished load
    logging.info("Load completed.")

    # benchmarking snippet
    urls = []
    for i, guid in enumerate(fn4c.guids()):
        # add urls; can choose what to use
        # urls.append("{0}/api/v2/{1}/exists".format(server_url, guid)) # one database access
        urls.append("{0}/api/v2/server_time".format(server_url))  # no database access
        # urls.append("{0}/api/v2/--Reference--/{1}/exact_distance".format(server_url, guid))     # 2 database access needed + a computation

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
