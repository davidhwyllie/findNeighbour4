#!/usr/bin/env python
""" loads samples into a findNeighbour4 server
from a multi-fasta file containing a simulated phylogeny.

performs PCA serially.

assumes a findNeighbour4 server is running
you can start one with

pipenv run python3 findNeighbour4_server.py demos/covidsim/config.json

Paths expected are currently hard-coded

Example usage: 
============== 
# show command line options 
python covidsim.py --help  

# example usage
pipenv run python3 demo/covidsim.py demos/covidsim/config.json /data/data/pca/sim/fasta  sim0.fasta --outputdir /data/data/pca/sim/sqlite --n_components 20

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
from collections import Counter
from fn4client import fn4Client
import sentry_sdk
from findn.mongoStore import fn3persistence
from findn.common_utils import ConfigManager
from pca.pca import VariantMatrix, PCARunner

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
                            pipenv run python3 demo/covidsim.py demos/covidsim/config.json /data/data/pca/sim/fasta  sim0.fasta --outputdir /data/data/pca/sim/sqlite --n_components 20


                            """,
    )

    parser.add_argument(
        "path_to_config_file",
        type=str,
        action="store",
        nargs="?",
        help="the path to the configuration file",
        default=None,
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

    parser.add_argument(
        "--outputdir",
        type=str,
        action="store",
        nargs="?",
        help="the directory in which the output will appear",
        default="",
    )
    parser.add_argument(
        "--n_components",
        type=int,
        action="store",
        nargs="?",
        help="the number of pcs to extract.  If this is smaller than the numbers of samples trained_on, instability is likely.  Default 20",
        default=20,
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

    print("Connecting to backend data store")
    try:
        PERSIST = fn3persistence(
            dbname=CONFIG["SERVERNAME"],
            connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=CONFIG["DEBUGMODE"],
            server_monitoring_min_interval_msec=0,
        )
    except Exception:
        print("Error raised on creating persistence object")
        raise
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
        logger.info("Resetting")
        fn4c.reset()
        existing_guids = fn4c.guids()
        print("After resetting, there remain {0} samples".format(len(existing_guids)))

        if len(existing_guids) > 0:
            print("Error.  More than 0 guids remain: {0}".format(len(existing_guids)))
            exit()

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
            records.append(res)

        logger.info(
            "Complete read of simulated data into RAM.  There are {0} sequences.  Sorting ...".format(
                len(records)
            )
        )

        seq_df = pd.DataFrame.from_records(records)
        seq_df.sort_values(by="sample_date", ascending=True, inplace=True)
        print(seq_df)

        destfile = os.path.join(args.outputdir, fastastem + ".meta.txt")
        with open(destfile, "wt") as f:
            seq_df.to_csv(f, date_format="YYYY-MM-DD", index=False, sep="\t")
        logging.info("Wrote metadata to {0}".format(destfile))

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

            if n_added > 200:
                # run pca
                try:
                    var_matrix = VariantMatrix(CONFIG, PERSIST)
                except Exception:
                    print("Error raised on instantiating Variant Matrix object")
                    raise

                print("Building snp matrix")
                try:
                    var_matrix.build()
                    # print("Running PCA on snp matrix")
                    pca_runner = PCARunner(var_matrix)
                    pca_runner.run(n_components=args.n_components, pca_parameters={})
                    vm = pca_runner.cluster()
                    vm.to_sqlite(
                        outputdir=args.outputdir,
                        analysis_name=fastastem + "_" + str(addition_date),
                    )
                except ValueError:
                    warnings.warn(
                        "Failed to run PCA - usually reflects inadequate amount of data"
                    )

    # finished
    logging.info("Finished, terminating program.")
