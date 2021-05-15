""" produces a file contain all edges (distances between guids)

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
import logging
import logging.handlers
import argparse
import progressbar
import csv

# config
from findn.common_utils import ConfigManager

# startup
from findn.mongoStore import fn3persistence

if __name__ == "__main__":

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""dump_edge_list.py - produces a list of stored pairwise distances.
                                     

Example usage: 
============== 

## does not require findNeighbour4_server to be running
Example:
pipenv run python3 dump_edge_list.py ../demos/covid/covid_config_v3.json  --outputdir /backup/edgelist

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
        help="the output directory.  Will try to make the directory if it does not exist",
    )
    args = parser.parse_args()

    ############################ LOAD CONFIG ######################################
    print("dump_edge_list.py .. reading configuration file.")
    print("Warning: take the server offline before doing this.  If the server is online, there is not guarantee that the dump edges operation will be atomic.")
    if len(args.path_to_config_file) > 0:
        config_file = args.path_to_config_file
        debugmode = False
        logging.info(config_file)
    else:
        raise FileNotFoundError("No config file name supplied ; this is required ")

    cfm = ConfigManager(config_file)
    CONFIG = cfm.read_config()

    ##################################################################################################
    # open PERSIST and hybridComparer object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed

    try:
        PERSIST = fn3persistence(
            dbname=CONFIG["SERVERNAME"],
            connString=CONFIG["FNPERSISTENCE_CONNSTRING"],
            debug=0,
        )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception:

        raise

    ########################### create output directory if it does not exist ##########
    if args.outputdir is not None:
        outputdir = args.outputdir
        logging.info("Writing output to {0}".format(outputdir))
    else:
        exit("Need to specify outputdir")
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    ####################################  START ################################################################

    #################################  prepare to iterate #########################################################
    show_bar = True
    all_samples = PERSIST.guids()
    n_samples = len(all_samples)
    n_links = 0
    if show_bar:
        bar = progressbar.ProgressBar(max_value=n_samples)

    with open(os.path.join(outputdir, "sample_list.csv"), "w") as f1:
        with open(os.path.join(outputdir, "edges_list.csv"), "w") as f2:
            sample_writer = csv.writer(f1, quoting=csv.QUOTE_NONNUMERIC)
            sample_writer.writerow(["sample_int_id", "sample_id"])
            edge_writer = csv.writer(f2, quoting=csv.QUOTE_NONNUMERIC)
            edge_writer.writerow(["edge_int_id", "sample_id_1", "sample_id_2", "dist"])

            for num_loaded, guid in enumerate(all_samples):
                if show_bar:
                    bar.update(num_loaded)
                sample_writer.writerow([num_loaded, guid])
                links = PERSIST.guid2neighbours(guid, returned_format=1)["neighbours"]

                for link in links:
                    n_links += 1
                    edge_writer.writerow([n_links, guid, link[0], link[1]])
                    if n_links % 50e6 == 0 and n_links > 0:
                        print(
                            "Loaded samples = {0}; edges = {1}; edges per sample = {2}".format(
                                num_loaded, n_links, int(n_links / num_loaded)
                            )
                        )
