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
* VariatioremovenModel - stores results of producing variant matrix and running PCA
* VariantMatrix - computes sample x variant matrix (requires: PERSIST object for mongodb access; server configuration file) 
* PCARunner - runs PCA on VariantMatrix

Example usage:
# write to sqlite db
pipenv run python3 fn4_pca.py demos/covid/covid_config_v3.json sqlite:////data/data/pca/fn4_pca3/2020-06-01.sqlite  --n_components 100 --focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca3 --analysis_date 2020-06-01 --remove_existing_data

# write to an oracle db identified by 'prod'
pipenv run python3 fn4_pca.py demos/covid/covid_config_v3.json prod  --n_components 100 --focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca3 --analysis_date 2020-06-01 --remove_existing_data
nohup pipenv run python3 fn4_pca.py demos/covid/covid_config_v3.json prod  --n_components 400 --focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca3 --analysis_date 2021-07-10 --remove_existing_data > realtime.out --fdr 0.05 &

see also main()

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.
cronta --
 

"""

#
# due to current scikit-learn/ open blas limitations,
# need to ensure that too many jobs are not started
# https://github.com/scikit-learn/scikit-learn/issues/20539; crashes if > 64 cpus
# unless set OMP_NUM_THREADS=64.  One test machine we tried has 104 cores.
# https://stackoverflow.com/questions/55531880/does-partial-fit-runs-in-parallel-in-sklearn-decomposition-incrementalpca
# need to do this before loading numpy.  Set it in the .env file
import os
import logging
import logging.handlers
import warnings
import datetime
from pathlib import Path
import sentry_sdk
from version import version
import argparse
import shutil
import progressbar
import platform

# reference based compression storage and clustering modules
import numpy as np
from findn.persistence import Persistence
from findn import DEFAULT_CONFIG_FILE
from findn.common_utils import ConfigManager
from findn.cw_seqComparer import cw_seqComparer
from pca.pca_scalable import VariantMatrix, PCARunner
from pca.pcadb import PCADatabaseManager
from pca.fittrend import ModelCounts
from localstore.localstoreutils import LocalStore
from tree.tree_utils import IQTree, ManipulateTree, DepictTree

# os.environ['OMP_NUM_THREADS']='64'          # required for KMeans with large numbers of samples, otherwise OpenBLAS crashes; set this in .env


def export_trees(metadata, target_dir, mt, pdm, has_controls, population_id, newick_tree):
    """export trees to database and file
    metadata: a pandas dataframe containing output data
    target_dir: where to write to
    mt: a ManipulateTree object
    pdm: a PCA data manager object
    has_controls: whether controls are included.  Either + or -
    population_id: integer reflecting population in database
    newick_tree: newick format tree
    """

    # output and write to database
    targetfile = os.path.join(target_dir, "rectangular.svg")
    mt.render(targetfile, mode="r")

    add_infos = []

    if os.path.exists(targetfile):
        with open(targetfile, "rt") as f:
            svgfile_content_r = f.read()
        add_infos.append(
            dict(
                pop_int_id=population_id,
                info_tag=has_controls + "r_svg",
                info_description="A tree including trending samples in a population in radial svg format.  Not rooted",
                mime_type="image/svg+xml",
                info_class="svg",
                info=svgfile_content_r,
            )
        )

    targetfile = os.path.join(target_dir, "circular.svg")
    mt.render(targetfile, mode="c")
    if os.path.exists(targetfile):
        with open(targetfile, "rt") as f:
            svgfile_content_c = f.read()
            add_infos.append(
                dict(
                    pop_int_id=population_id,
                    info_tag=has_controls + "c_svg",
                    info_description="A tree including trending samples in a population in radial svg format. Not rooted",
                    mime_type="image/svg+xml",
                    info_class="svg",
                    info=svgfile_content_c,
                )
            )

    # write data to database
    targetfile = os.path.join(target_dir, "meta.csv")
    with open(targetfile, "wt") as f:
        metadata.to_csv(f, index=True, index_label="sample_id")
    with open(targetfile, "rt") as f:
        csvfile_content = f.read()

    add_infos.append(
        dict(
            pop_int_id=population_id,
            info_tag=has_controls + "c_newick",
            info_description="A tree including all the trending samples in a population in newick format.  Not rooted",
            mime_type="text/x-nh",
            info_class="newick",
            info=newick_tree,
        )
    )
    add_infos.append(
        dict(
            pop_int_id=population_id,
            info_tag=has_controls + "c_metadata",
            info_description="Metadata on all the trending samples in a population.  Not rooted",
            mime_type="text/csv",
            info_class="csv",
            info=csvfile_content,
        )
    )

    for add_info in add_infos:
        pdm.add_PopulationStudiedExtraInfo(**add_info)


def main():
    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""Runs findneighbour4_pca, which extract pcs from a sequence collection

If provided with a COG-UK format sequence file, will estimates rates of change of different sequence components as new variant detection
technology.

Example usage: 
============== 
# show command line options 
python fn4_pca.py --help  

# run with debug settings; only do this for unit testing.
python fn4_pca.py     

# run using settings in a configuration file

# example with data written sqlite
pipenv run python3 fn4_pca.py demos/covid/atp.json sqlite:////data/data/pca/fn4_pca3/2020-06-01TEST.sqlite /data/logs/findNeighbour4/localcache/fndev_atptest/pca --n_components 100 --focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca3 --analysis_date 2020-06-01 --remove_existing_data

# example with data written to oracle database identified by 'prod'
pipenv run python3 fn4_pca.py demos/covid/atp.json prod /data/logs/findNeighbour4/localcache/fndev_atptest/pca --n_components 100 --focus_on_most_recent_n_days 60 --compute_slope_over 30 --analysis_dir /data/data/pca/fn4_pca3 --analysis_date 2020-06-01 --remove_existing_data

########## End of old examples


""",
    )
    parser.add_argument(
        "path_to_config_file",
        type=str,
        action="store",
        nargs="?",
        help="the path to the fn4 server configuration file",
        default=None,
    )
    parser.add_argument(
        "connection_config",
        type=str,
        action="store",
        nargs="?",
        help="the key in the database credentials found in the json file at PCA_CONNECTION_CONFIG_FILE. OR a database configuration string e.g. sqlite:///mydb",
    )
    parser.add_argument(
        "storage_dir",
        action="store",
        help="a directory to store pre-processed sequences in.",
    )
    parser.add_argument(
        "--analysis_date",
        type=str,
        action="store",
        nargs="?",
        help="the date of the analysis.  Later samples are not considered ",
        default=datetime.date.today().isoformat(),
    )
    parser.add_argument(
        "--analysis_window_start_date",
        type=str,
        action="store",
        nargs="?",
        help="the samples earlier than this are not considered ",
        default=datetime.date(2020, 1, 1),
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
        "--min_variant_frequency",
        type=float,
        action="store",
        nargs="?",
        help="the minimum variant frequency observed to include in models.  If an integer, is assumed to be the minimum variant count acceptable.",
        default=10,
    )
    parser.add_argument(
        "--focus_on_most_recent_n_days",
        type=int,
        action="store",
        nargs="?",
        help="only analyse pc_cats which arose this number of days prior to the analysis date",
        default=120,
    )
    parser.add_argument(
        "--compute_slope_over",
        type=int,
        action="store",
        nargs="?",
        help="compute per-pc_cat incidence rate ratio over this number of dates prior to analysis_date",
        default=30,
    )
    parser.add_argument(
        "--fdr",
        type=float,
        action="store",
        nargs="?",
        help="Benjamini Hochberg False Discovery Rate for multiple testing of trends in pc_cats",
        default=0.01,
    )
    parser.add_argument(
        "--cogfile",
        type=str,
        action="store",
        nargs="?",
        help="the name of the cog-uk format data file to load.",
        default="/data/data/inputfasta/cog_metadata.csv",
    )
    parser.add_argument(
        "--remove_existing_data",
        action="store_true",
        default=False,
        help="if this option is specified, it will delete any existing data from the database. ",
    )
    parser.add_argument(
        "--analysis_dir",
        action="store",
        default="/tmp",
        help="a temporary directory to write trees etc into ",
    )
    parser.add_argument(
        "--remove_temporary_trees",
        action="store_true",
        default=False,
        help="if this option is specified, it will delete any temporary files generated in --analysis_dir ",
    )
    parser.add_argument(
        "--only_produce_tree_output",
        action="store_true",
        default=False,
        help="debug setting.  Doesn't do pca, just runs output (including tree depiction) on the last run",
    )

    args = parser.parse_args()

    ############################ LOAD CONFIG ######################################
    fdr = args.fdr
    if fdr < 0 or fdr > 1:
        raise ValueError("FDR must be between 0 and 1 not {0}".format(fdr))

    logging.info("findNeighbour4 PCA modelling .. system info.")
    logging.info(platform.platform())
    logging.info(platform.machine())
    logging.info(platform.processor())
    logging.info(platform.python_version())
    logging.info("Numpy config:")
    logging.info(np.show_config())

    logging.info("reading findneighbour configuration file.")

    config_file = args.path_to_config_file
    if config_file is None:
        config_file = DEFAULT_CONFIG_FILE
        warnings.warn(
            f"No config file name supplied; using configuration in {DEFAULT_CONFIG_FILE}, suitable only for testing, not for production."
        )

    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()

    # create an analysis dir if not present
    os.makedirs(args.analysis_dir, exist_ok=True)

    # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
        CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
        logging.info("Set mongodb connection string  from environment variable")
    else:
        logging.info("Using mongodb connection string from configuration file.")

    # determine whether a FN_SENTRY_URL environment variable is present,
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
    logfile = os.path.join(
        logdir, "pca-{0}".format(os.path.basename(CONFIG["LOGFILE"]))
    )
    logging.info("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(
        logfile, mode="a", maxBytes=1e7, backupCount=7
    )

    formatter = logging.Formatter(
        "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s "
    )
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
        sentry_sdk.init(CONFIG["SENTRY_URL"], release=version)

    # prepare connection to database
    logging.info("Connecting to backend data store")
    pm = Persistence()
    PERSIST = pm.get_storage_object(
        dbname=CONFIG["SERVERNAME"],
        connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
        debug=CONFIG["DEBUGMODE"],
        verbose=True,
    )

    # determine whether there is a localstore (tar file) of sequence data available.  If there is a localstore will use it because it's faster
    tarfile_name = os.path.join(cfm.rcscache, "rcs.tar")
    if os.path.exists(tarfile_name):
        SPERSIST = LocalStore(tarfile_name)
        logging.info(
            "Using local sequence store {0} as a source of sequence data ".format(
                tarfile_name
            )
        )
    else:
        SPERSIST = PERSIST  # use database connection
        logging.info(
            "Using database as a source of sequence data, as local tar file not found"
        )

    if args.remove_existing_data:
        logging.info("fn4_pca is set to remove all existing data from database")
    else:
        logging.info("fn4_pca is set to add build to existing data")

    # note today's date.  Used for sample selection & poisson modelling
    analysis_date = datetime.date.fromisoformat(args.analysis_date)
    analysis_window_start_date = datetime.date.fromisoformat(
        args.analysis_window_start_date
    )

    logging.info("Connecting to database")
    pdm = PCADatabaseManager(
        connection_config=args.connection_config, debug=args.remove_existing_data
    )

    logging.info("Setting up a cw_seqComparer object for sequence data access")
    # this is used for multisequence alignments late in the process
    hc = cw_seqComparer(
        reference=CONFIG["reference"],
        maxNs=CONFIG["MAXN_STORAGE"],
        snpCeiling=CONFIG["SNPCEILING"],
        excludePositions=set(CONFIG["excludePositions"]),
        preComparer_parameters={},
        PERSIST=PERSIST,
        disable_insertion=True,
    )

    if args.only_produce_tree_output is False:
        samples_added = None
        logger.info("Loading cog-uk metadata")
        samples_added = pdm.store_cog_metadata(
            cogfile=args.cogfile,
            date_start=analysis_window_start_date,
            date_end=analysis_date,
        )

        try:
            var_matrix = VariantMatrix(
                CONFIG, SPERSIST, args.storage_dir, show_bar=True
            )
        except Exception:
            logging.info("Error raised on instantiating Variant Matrix object")
            raise

        logging.info("Preparing matrix-format sequence data pre-analysis")
        var_matrix.prepare_to_analyse()

        logging.info("Running PCA on snp matrix")
        pca_runner = PCARunner(var_matrix)
        pca_runner.run(
            n_components=args.n_components,
            select_from=samples_added,
            min_variant_freq=args.min_variant_frequency,
            pca_parameters={},
        )

        pca_runner.cluster()

        logging.info("Storing variation model and PCA")
        pdm.store_variation_model(pca_runner.vm.vm)

        logging.info("Building contingency tables, relating Lineage to pc_cat")
        pdm.make_contingency_tables(
            only_pc_cats_less_than_days_old=args.focus_on_most_recent_n_days,
            today=analysis_date,
        )

        logging.info("Storing PCA summary")
        pdm.store_pca_summary()  # store a summary

        pcas_df = pdm.pca_summary(
            only_pc_cats_less_than_days_old=args.focus_on_most_recent_n_days,
            today=analysis_date,
        )
        pcas_df = pcas_df[pcas_df["n_days_observed"] >= 3]
        n_pc_cats = len(pcas_df["pc_cat"].unique())
        logging.info(
            "Fitting poisson models.  There are {0} models to fit over {1} pc_cats originating in the last {2} days".format(
                len(pcas_df.index), n_pc_cats, args.focus_on_most_recent_n_days
            )
        )

        bar = progressbar.ProgressBar(max_value=len(pcas_df.index))
        for i, pcas_int_id in enumerate(pcas_df.index):
            bar.update(i)

            # ---------------------------------

            pcas_obj = pdm.single_pcas_summary(pcas_int_id)
            cntdata = pdm.pcas_count_table(pcas_obj, output_format=1)

            # compute the date range over which the slope should be computed.
            # the last date is addition_date
            # the first date is either
            #  focus_on_most_recent_n_days before addition_date, OR
            #  the date the pc_cat was first seen, if that is later.
            analysis_start = analysis_date - datetime.timedelta(
                days=args.compute_slope_over
            )

            if cntdata["earliest_date"] > analysis_start:
                analysis_start = cntdata["earliest_date"]

            # specify latest date to be modelled
            mc = ModelCounts(
                counts=cntdata["counts"],
                denominators=cntdata["denominators"],
                pcas_int_id=pcas_int_id,
                pc_cat=pcas_obj.pc_cat,
                earliest_date=analysis_start,
                latest_date=analysis_date,
            )

            res = mc.fit_poisson()
            pdm.store_pcas_model_output(res)

            # --------------------------------------

        bar.finish()

    trending_details = pdm.trending_samples_metadata(
        max_size_of_trending_pc_cat=500,
        p_value_cutoff=fdr,
        filtering_method="bh",
        analysis_type="GLM:Poisson regression",
    )

    if trending_details is None:
        # nothing found to analyse
        logging.info("Nothing found to analyse")
        exit(0)

    logging.info("Identified trending pc_cats in the following:")
    logging.info(trending_details["population_meta"])

    # for each trending population, build a tree
    # optional depictions
    logging.info("Generation depictions for each population")

    iq = IQTree(genome_length=len(CONFIG["reference"]))
    # within each population, we get examples of the expanding samples, and controls
    for population_id in trending_details["population_annotations"].keys():
        population_id = int(population_id)  # not np.int64
        population_info = pdm.single_population_studied_from_pop_int_id(population_id)
        this_build_int_id = int(population_info.build_int_id)

        # make analysis dir
        print(
            "**************************************************************************************"
        )
        analysis_dir = os.path.join(
            args.analysis_dir, str(this_build_int_id), str(population_id)
        )
        os.makedirs(analysis_dir, exist_ok=True)

        ## just those selected, up to 1000 samples
        to_tree_build = trending_details["population_annotations"][population_id].copy()
        to_tree_build = to_tree_build.sort_values(
            "sample_date", ascending=False
        )  # get most recent samples
        print(to_tree_build)
        to_tree_build["expanding"] = "Y"

        exp_only = to_tree_build.head(1000)
        exp_only = exp_only.fillna("-")  # metadata
        logging.info(
            "Population # {1}: Selected expanding population (max 1000 samples) for tree build n = {0}".format(
                len(exp_only.index), population_id
            )
        )

        # don't root
        sample_ids = exp_only.index.to_list()
        sample_ids.append("--Reference--")
        msa_result = hc.multi_sequence_alignment(sample_ids)

        target_dir = os.path.join(analysis_dir, "tree_no_controls")
        iqr = iq.build(
            msa_result.msa_fasta(), msa_result.fconst, targetdir=target_dir
        )  # note: if directory is hard coded, will generate conflict if two processes run at the same time
        newick_tree = iqr["output"]["newick"]
        rt = ManipulateTree(newick_tree)
        rt.reroot("--Reference--")
        rt.remove_outgroup()
        newick_tree = rt.newick()
        iqr["output"]["newick"] = newick_tree  # store the re-rooted tree back

        title_info = "{0} {1}; {2} {3}".format(
            population_info.level_1_category_type,
            population_info.level_1_category,
            population_info.level_2_category_type,
            population_info.level_2_category,
        )
        mt = DepictTree(
            newick_tree,
            exp_only,
            title=[title_info],
            genome_length=len(CONFIG["reference"]),
        )

        # export
        export_trees(
            metadata=exp_only, 
            target_dir = target_dir, 
            mt= mt, 
            pdm = pdm, 
            has_controls="-",
            population_id = population_id, 
            newick_tree= newick_tree)
        
        ## include controls as well
        print(
            "**************************************************************************************"
        )
        pop_samples = pdm.population_members(
            population_id, max_rows=5000
        )  # comes ordered by date
        pop_samples = pop_samples.set_index("sample_id")
        pop_samples["expanding"] = "N"

        to_tree_build = trending_details["population_annotations"][population_id].copy()
        to_tree_build = to_tree_build.sort_values(
            "sample_date", ascending=False
        )  # get most recent samples
        to_tree_build["expanding"] = "Y"

        # drop any samples from pop_samples which are in the relevant pcs
        s_pop_samples = set(pop_samples.index.to_list())
        s_to_tree_build = set(to_tree_build.index.to_list())
        to_drop = s_pop_samples.intersection(s_to_tree_build)
        pop_samples = pop_samples.drop(to_drop)

        exp_and_control = pop_samples.head(300).append(to_tree_build.head(300))
        exp_and_control = exp_and_control.fillna("-")  # metadata
        logging.info(
            "Population # {1}: Selected expanding population and control (max 300 samples each) for tree build n = {0}".format(
                len(exp_and_control.index), population_id
            )
        )

        # build MSA including reference seq which we'll use to root.  We will then remove the reference.
        sample_ids = exp_and_control.index.to_list()
        sample_ids.append("--Reference--")
        msa_result = hc.multi_sequence_alignment(sample_ids)

        target_dir = os.path.join(analysis_dir, "tree_with_controls")
        iqr = iq.build(
            msa_result.msa_fasta(), msa_result.fconst, targetdir=target_dir
        )  # note: if directory is hard coded, will generate conflict if two processes run at the same time
        newick_tree = iqr["output"]["newick"]
        rt = ManipulateTree(newick_tree)
        rt.reroot("--Reference--")
        rt.remove_outgroup()
        newick_tree = rt.newick()
        iqr["output"]["newick"] = newick_tree  # store the re-rooted tree back

        # outputfile =  "testdata/ete3/test{0}.nwk".format(population_id)
        # with open(outputfile, "wt") as f:
        #    f.write(newick_tree)

        # outputfile = "testdata/ete3/test{0}.pickle".format(population_id)
        # with open(outputfile, "wb") as f:
        #    pickle.dump(exp_and_control, f)

        title_info = "{0} {1}; {2} {3}".format(
            population_info.level_1_category_type,
            population_info.level_1_category,
            population_info.level_2_category_type,
            population_info.level_2_category,
        )
        mt = DepictTree(
            newick_tree,
            exp_and_control,
            title=[title_info],
            genome_length=len(CONFIG["reference"]),
        )
        export_trees(
            metadata=exp_and_control, 
            target_dir = target_dir, 
            mt= mt, 
            pdm = pdm, 
            has_controls="+",
            population_id = population_id,
            newick_tree = newick_tree)

        # cleanup
        if args.remove_temporary_trees:
            print("Deleting temporary files")
            shutil.rmtree(analysis_dir)  # is this vulnerable to symlink attack? TBD

    logging.info("Build finished.  Results are in database.  Software will terminate shortly, but may take 1-2 minutes to run .tar file validity checking first.")


# startup
if __name__ == "__main__":
    main()
