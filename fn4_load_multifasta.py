""" loads samples into a findNeighbour4 server
from a muliti-fasta file (e.g. the COVID-19 alignment produced by COG-UK)

assumes a findNeighbour4 server is running

Paths expected are currently hard-coded

Example usage: 
============== 
# show command line options 
python fn4_load.py --help  

python3 fn4_load_multifasta.py [server_url] [directory to look for fastas in] [filename to glob within directory]
# example usage
pipenv run python3 fn4_load_multifasta.py http://localhost:5023 /srv/data/covid COVID_MSA*.fasta

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
import dateutil
import shutil
import Bio
import logging
import logging.handlers
import argparse
import time
from collections import Counter
from fn4client import fn4Client
import sentry_sdk
from version import version

## define functions and classes
class DatabaseMonitorInoperativeError(Exception):
    """ insert failed """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


def measure_sr(fn4_client, ensure_database_monitor=False):
    """checks the most recently recorded server ratio (ratio of  rows in guid2neighbour to guid2meta),
    and notes whether this was in the last 10 minutes

    it it wasn't it either
    - raises a DatabaseMonitorInoperativeError (if ensure_database_monitor is True) or
    - issues are warning, and continues, reporting storage_ratio as 1

    Parameters
    ensure_database_monitor - either True or False.  If true, raises an error if
    """
    seconds_ago = None
    server_database_usage = fn4_client.server_database_usage(nrows=1)
    server_time_now = fn4_client.server_time()
    if "trend_stats" in server_database_usage.keys():
        report_time = server_database_usage["trend_stats"].loc[0, "context|time|time_now"]
        print(server_database_usage)

        td = dateutil.parser.parse(server_time_now["server_time"]) - dateutil.parser.parse(
            report_time
        )  # how long ago was the report on the database?
        seconds_ago = td.total_seconds()
        logging.info("Checked storage ratio.  Last measurement was {0} seconds ago".format(seconds_ago))

    if seconds_ago is None or seconds_ago > 600:  # 10 mins ago
        if ensure_database_monitor:
            raise DatabaseMonitorInoperativeError(
                "No recent measurements of server health indicate findNeighbour4_dbmanager is not operational",
                "Last measurement was {0} seconds ago (None = no record of ever measurement)".format(seconds_ago),
            )
        else:
            logging.warning(
                "No recent measurements of server health indicate findNeighbour4_dbmanager is not operational",
                "Last measurement was {0} seconds ago (None = no record of ever measurement)".format(seconds_ago),
            )

            server_database_usage["latest_stats"]["storage_ratio"] = 1  # not sure what it is now - continue
    return server_database_usage["latest_stats"]["storage_ratio"]


if __name__ == "__main__":

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_server, a service for bacterial relatedness monitoring.
                                     

Example usage: 
============== 
# show command line options 
python fn4_load_multifasta.py --help  
"""
    )
    parser.add_argument(
        "server_url",
        type=str,
        action="store",
        nargs="?",
        help="the url (and port) of the findNeighbour4 server into which we should insert",
        default="",
    )
    parser.add_argument(
        "fastadir", type=str, action="store", nargs="?", help="the directory in which fasta files will appear", default=""
    )
    parser.add_argument("fileglob", type=str, action="store", nargs="?", help="a pattern to glob for", default="")
    args = parser.parse_args()

    ####################################    STARTUP ###############################
    # validate input
    if not os.path.exists(args.fastadir):
        # that's an error
        raise FileNotFoundError("The directory specified for fasta files does not exist: '{0}'".format(args.fastadir))
    else:
        if not os.path.isdir(args.fastadir):
            raise FileNotFoundError("The path specified for fasta files is not a directory: '{0}'".format(args.fastadir))

    # the fasta dir exists.  Make sure we have log directories.
    logdir = os.path.join(args.fastadir, "logs")
    completedir = os.path.join(args.fastadir, "completed")
    for testdir in [logdir, completedir]:
        os.makedirs(testdir, exist_ok=True)

    # launch logger
    logger = logging.Logger("fn4_load_multifasta")
    logger.setLevel(logging.INFO)
    timenow = datetime.datetime.now().isoformat()
    logfile = os.path.join(logdir, "fn4_load_multifasta_v2.log")
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode="a", maxBytes=1e7, backupCount=7)
    formatter = logging.Formatter("%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info("Startup with arguments: {0}".format(args))

    print("To see what is happening, do watch tail {0}".format(logfile))
    # launch comms with Sentry
    if os.environ.get("FN_SENTRY_URL") is not None:
        logger.info("Launching communication with Sentry bug-tracking service")
        sentry_sdk.init(os.environ.get("FN_SENTRY_URL"), release=version)
        logger.info("Sentry comms established")
    else:
        logger.info("No error monitoring via sentry.  Set environment variable FN_SENTRY_URL to do so.")

    # instantiate client
    fn4c = fn4Client(args.server_url)  # expects operation on local host; pass baseurl if somewhere else.

    existing_guids = set(fn4c.guids())
    clustering_created = False
    logger.info("There are {0} existing guids".format(len(existing_guids)))

    # add the reference sequence as the root if not already present
    ref_guid = "--Wuhan-Reference--"
    ref_guid_present = fn4c.guid_exists(ref_guid)
    for record in Bio.SeqIO.parse("reference/nc_045512.fasta", "fasta"):
        refseq = str(record.seq).upper()
        reflen = len(refseq)

    if not ref_guid_present:
        logger.info("Adding reference")
        res = fn4c.insert(guid=ref_guid, seq=refseq)
    else:
        logger.info("Reference already present")

    fastafiles = sorted(glob.glob(os.path.join(args.fastadir, args.fileglob)))
    logger.info("There are {0} fasta files waiting to be loaded".format(len(fastafiles)))

    ## optionally raise an error if this list is too long
    for fastafile in fastafiles:
        logger.info("Scanning {0}".format(fastafile))

        nSkipped = 0
        nBad = 0
        nGood = 0
        failed = []
        i = 0

        for record in Bio.SeqIO.parse(fastafile, "fasta"):

            i = i + 1
            t1 = datetime.datetime.now()
            namebits = record.id.split("/")

            if len(namebits) == 1:
                guid = namebits[0]  # format CAMC-12CD1B14
            else:
                guid = namebits[1]  # format England/CAMC-12C1B14/2021

            guid = guid.replace("/", "-")
            guid = guid.replace(":", "-")

            seq = str(record.seq).upper()
            seq = seq.replace("?", "N")
            seq = seq.replace(" ", "N")
            counter = Counter(list(seq))

            res = {"record_id": record.id, "guid": guid, "seqlen": len(seq), **counter}
            if guid not in existing_guids:

                if len(seq) == reflen:
                    res = fn4c.insert(guid=guid, seq=seq)
                    if not res.status_code == 200:
                        # failed to add
                        failed.append((i, guid))
                        msg = "** FAILED **"

                    else:
                        msg = "succeeded"
                        nGood += 1

                    t2 = datetime.datetime.now()
                    i1 = t2 - t1
                    s = i1.total_seconds()
                    logger.info("Scanned {0} Adding #{1} ({2}) {3} in {4:.2f} secs.".format(i, nGood, guid, msg, s))

                    # build in pause if high storage ratio ('fragmentation')
                    if nGood % 50 == 0 and nGood > 0:

                        # check whether database is keeping repacked adequately.  If it isn't, insertion will pause.
                        # If we don't know, because monitoring is off, then an error will be raised.
                        sr = measure_sr(fn4c)
                        logger.info(
                            "Examined {0} / skipped {1}.  Database neighbour fragmentation is {2:.1f} (target: 1)".format(
                                i, nSkipped, sr
                            )
                        )
                        # ratio of records containing neighbours to those containing samples - target is 1:1
                        while sr > 100:  
                            logger.warning(
                                "Waiting 6 minutes to allow repacking operations.  Will resume when fragmentation, which is now {0:.1f}, is < 100.".format(
                                    sr
                                )
                            )
                            time.sleep(360)  # restart in 6 mins  if below target
                            sr = measure_sr(fn4c)
                else:
                    logger.info("{0} Wrong length {1} not {2}".format(guid, len(seq), reflen))

            else:
                nSkipped += 1

            if i % 2500 == 0:
                logger.info("Examined {0} / skipped {1}".format(i, nSkipped))

        logger.info(
            "Complete.  Skipped {0} guids which already exist in the server.  There are {1} bad sequences".format(
                nSkipped, nBad
            )
        )

        if len(failed) > 0:
            # do not move the file
            logger.warning("findneighbour_load | not all fasta files could be uploaded; {0} failed".format(len(failed)))
            logger.info("Failed samples are: {0}".format(failed))
        else:
            logger.info("Sample load succeeded.  Moving fasta file to /completed directory")
            shutil.move(fastafile, completedir)

        # finished
        logging.info("Finished, terminating program.")
