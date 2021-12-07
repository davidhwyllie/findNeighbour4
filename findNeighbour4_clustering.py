""" clustering for findNeighbour4
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
from findn.py_seqComparer import py_seqComparer
from findn.common_utils import ConfigManager

from snpclusters.ma_linkage import (
    MixtureAwareLinkage,
    MixPOREMixtureChecker,
    MixtureAwareLinkageResult,
)
from snpclusters.clusternomenclature import ClusterNameAssigner, ClusterNomenclature


if __name__ == "__main__":

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_clustering, a findNeighbour4 component.
                                     
Example usage: 
============== 

## does not require findNeighbour4_server to be running
Minimal example:
python findNeighbour4_clustering.py config/myconfig_file.json  

# Relabel clusters
nohup pipenv run python3 findNeighbour4_clustering.py ../phe_dev/config_phe_dev.json --rebuild_clusters_debug --label_clusters_using=reference/guid2cluster.xlsx &

if a config file is not provided, it will run (as does findNeighbour4_server) is debug mode: it will run once, and then terminate.  This is useful for unit testing.  If a config file is specified, the clustering will  run until terminated.  

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
        "--rebuild_clusters_debug",
        help="delete existing, and rebuild, clusters.  Only for use in a debug setting",
        action="store_true",
    )
    parser.add_argument(
        "--label_clusters_using",
        help="label clusters based on existing membership, as in the two column excel file referred to.  The program will stop after the initial run if this option is used, and will need to be restarted without this option to function normally",
        action="store",
    )
    parser.add_argument(
        "--run_once",
        help="run once only.  The program will stop after the initial run if this option is used, and will need to be restarted without this option to function normally",
        action="store_true",
    )
    args = parser.parse_args()

    # an example config file is default_test_config.json

    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 clustering .. reading configuration file.")

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

    ########################### read file containing labels, if it exists ##########
    relabel = False
    if args.label_clusters_using is not None:
        if os.path.exists(args.label_clusters_using):
            print(
                "File of cluster labels {0} exists; reading it".format(
                    args.label_clusters_using
                )
            )
            label_df = pd.read_excel(args.label_clusters_using)
            relabel = True
            existing_labels = set()
            # set column headers
            label_df.columns = ["original_label", "guid"]
            # drop anything which empty labels
            label_df.dropna(inplace=True)
            label_df["label"] = [x[:5] for x in label_df["original_label"]]
            previous_guid2cluster_label = {}
            for ix in label_df.index:
                existing_labels.add(label_df.at[ix, "label"])
                previous_guid2cluster_label[label_df.at[ix, "guid"]] = [
                    label_df.at[ix, "label"]
                ]
            existing_labels = list(existing_labels)
            logging.info(
                "Relabelling samples based on existing labels. In existing data, {0} guids have {1} labels".format(
                    len(previous_guid2cluster_label), len(existing_labels)
                )
            )

    ########################### prepare to start clustering ####################################
    # construct the required global variables

    logger.info("Connecting to backend data store")
    pm = Persistence()
    PERSIST = pm.get_storage_object(
        dbname=CONFIG["SERVERNAME"],
        connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
        debug=0,
        verbose=True,
    )

    if args.rebuild_clusters_debug:
        logger.warning(
            "Wiping existing clustering data as --rebuild_clusters_debug is set"
        )
        PERSIST._delete_existing_clustering_data()
        logger.warning("Wiped existing clustering data")
    else:
        logger.info("Working with existing data ... ")

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

    ################################# clustering #############################################
    # open PERSIST and cw_seqComparer object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed

    # get a clustering object's settings
    logger.info("Creating clustering objects ...")
    clusterers = {}
    clusternomenclature = {}
    clusternameassigner = {}
    for clustering_name in CONFIG["CLUSTERING"].keys():
        clustering_setting = CONFIG["CLUSTERING"][clustering_name]

        mpmc = MixPOREMixtureChecker(
            hc, **clustering_setting
        )  # uses cw_seqComparer to load samples and compute msas

        # check update adds remaining guids
        logger.info("Creating clustering object {0}".format(clustering_name))
        clusterers[clustering_name] = MixtureAwareLinkage(
            PERSIST=PERSIST,
            MIXCHECK=mpmc,
            mixed_sample_management=clustering_setting["mixed_sample_management"],
            snv_threshold=clustering_setting["snv_threshold"],
            serialisation=None,
            parameters=clustering_setting,
            name=clustering_name,
        )
        clusterers[clustering_name].remove_legacy()  # remove any old versions
        logger.info("Created clustering object {0}".format(clustering_name))
        # if applicable, make a cluster nomenclature object;
        if "cluster_nomenclature_method" in clusterers[clustering_name].parameters:
            if not relabel:  # otherwise, we use the labels provided to us
                existing_labels = clusterers[clustering_name].existing_labels()
            else:
                logging.info("Using pre-supplied existing labels")
            # we are instructed to do cluster naming
            clusternomenclature[clustering_name] = ClusterNomenclature(
                cluster_nomenclature_method=clusterers[clustering_name].parameters[
                    "cluster_nomenclature_method"
                ],
                existing_labels=existing_labels,
            )
            clusternameassigner[clustering_name] = ClusterNameAssigner(
                clusternomenclature[clustering_name]
            )
            logger.info(
                "Created name assigner {0} with method {1}".format(
                    clustering_name,
                    clusterers[clustering_name].parameters[
                        "cluster_nomenclature_method"
                    ],
                )
            )

    # now iterate - on an endless loop
    while True:
        whitelist = set()
        nbuilt = 0
        for clustering_name in CONFIG["CLUSTERING"].keys():  # DEBUG
            clustering_setting = CONFIG["CLUSTERING"][clustering_name]
            clusterers[clustering_name].update()
            clusterers[clustering_name].cluster()

            # if applicable, generate and store cluster labels
            if "cluster_nomenclature_method" in clusterers[clustering_name].parameters:
                # we are instructed to do cluster naming
                cluster2guid = clusterers[clustering_name].cluster2names
                if (
                    not relabel
                ):  # then we use the existing labels as our source, otherwise we use the labels already computed
                    previous_guid2cluster_label = clusterers[
                        clustering_name
                    ].guid2cluster_labels()
                else:
                    logging.info("Using pre-supplied guid to cluster lookup")
                clusterid2clusterlabel = clusternameassigner[
                    clustering_name
                ].assign_new_clusternames(
                    clusterid2guid=clusterers[clustering_name].cluster2names,
                    previous_guid2cluster_label=previous_guid2cluster_label,
                )
                clusterers[clustering_name].apply_cluster_labels(clusterid2clusterlabel)

            # store the output
            clusterers[clustering_name].persist(what="graph")
            clusterers[clustering_name].persist(what="output")
            clusterers[clustering_name].remove_legacy()  # remove old versions

            # read the result
            malr = MixtureAwareLinkageResult(PERSIST=PERSIST, name=clustering_name)

            ms = MSAStore(
                PERSIST=PERSIST, in_ram_persistence_time=60
            )  # persist result for up to 60 seconds

            # recover existing msas
            stored_msa = ms.existing_tokens()

            # build multisequence alighments for existing tokens.  These can be built on the fly, but precomputing thom speeds up GUI performance
            logger.info("Precomputing clusters for {0}".format(clustering_name))
            cluster_contents = malr.cluster2guid.values()
            bar = progressbar.ProgressBar(max_value=len(cluster_contents))

            for i, guids in enumerate(cluster_contents):

                bar.update(i + 1)
                # token identifies cluster contents, nature of analysis, and outgroup
                token = ms.get_token(
                    malr.parameters["uncertain_base_type"], False, guids
                )
                if len(guids) > 2:
                    whitelist.add(token)  # we need to retain this msa, if it exists

                    if token not in stored_msa:  # if we haven't already computed it

                        # estimate p1 if we have not done so already.
                        hc.update_p1_estimate(
                            sample_size=100,
                            uncertain_base_type=malr.parameters["uncertain_base_type"],
                        )

                        # construct a seqcomparer to do an msa
                        sc = py_seqComparer(
                            reference=CONFIG["reference"],
                            maxNs=CONFIG["MAXN_STORAGE"],
                            snpCeiling=CONFIG["SNPCEILING"],
                            excludePositions=CONFIG["excludePositions"],
                        )

                        for guid in guids:
                            seq = PERSIST.refcompressedsequence_read(guid)
                            sc.persist(seq, guid)

                        msa_result = sc.multi_sequence_alignment(
                            guids,
                            expected_p1=hc.estimated_p1,
                            uncertain_base_type=malr.parameters["uncertain_base_type"],
                        )

                        ms.persist(token, msa_result)
                        nbuilt += 1

            bar.finish()

            # cleanup anything we don't need, including old clustering versions and msas
            ms.unpersist(whitelist=whitelist)
            logger.info(
                "Cleanup complete.  Stored data on {0} MSAs; Built {1} new clusters".format(
                    len(whitelist), nbuilt
                )
            )

        if debugmode | relabel | args.run_once:
            PERSIST.closedown()
            exit(0)

        logger.info("Waiting 60 seconds")
        time.sleep(60)
