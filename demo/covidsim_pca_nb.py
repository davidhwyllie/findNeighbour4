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
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

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
import progressbar
import pandas as pd
from collections import Counter
from fn4client import fn4Client
import sentry_sdk
from findn.mongoStore import fn3persistence
from findn.common_utils import ConfigManager
from pca.pca import VariantMatrix, PCARunner
from pca.pcadb import PCADatabaseManager
from pca.fittrend import ModelCounts

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
        description="""Fits models to existing pca runs, stored in sqlite databases
                    """
    )

    parser.add_argument(
        "--sqlitepath",
        type=str,
        action="store",
        nargs="?",
        help="the path to the configuration file",
        default="/data/data/pca/sim/sqlite/sim0_*.sqlite"
    )
    args = parser.parse_args()
    focus_on_most_recent_n_days = 120
    compute_slope_over = 120

    for i,sqlite_file in enumerate(sorted(glob.glob(args.sqlitepath))):

        addition_date = sqlite_file[31:41]
        addition_date = datetime.date.fromisoformat(addition_date)

        stem = os.path.basename(sqlite_file)[0:3]

        outputdir = os.path.dirname(sqlite_file)
        print(i, sqlite_file, addition_date, outputdir)
        
        pdm = PCADatabaseManager(
                connection_config="sqlite:///{0}".format(sqlite_file), debug=False, show_bar=False
            )
        pdm.remove_statistical_models()     # remove any old ones
        
        pcas_df = pdm.pca_summary(
            only_pc_cats_less_than_days_old=focus_on_most_recent_n_days,
            today = addition_date
        )
        pcas_df = pcas_df[pcas_df["n_days_observed"] >= 3]
        n_pc_cats = len(pcas_df["pc_cat"].unique())
        logging.info(
            "Fitting  Negative binomial models.  There are {0} models to fit over {1} pc_cats originating in  {2} days to {3}".format(
                len(pcas_df.index),
                    n_pc_cats, 
                    focus_on_most_recent_n_days,
                    addition_date
            )
        )

        bar = progressbar.ProgressBar(max_value=len(pcas_df.index))
        for i, pcas_int_id in enumerate(pcas_df.index):
            bar.update(i)
            #if i> 50:
            #    break

            pcas_obj = pdm.single_pcas_summary(pcas_int_id)

            cntdata = pdm.pcas_count_table(pcas_obj, format = 2)

            # compute the date range over which the slope should be computed.
            # the last date is addition_date
            # the first date is either
            #  focus_on_most_recent_n_days before addition_date, OR
            #  the date the pc_cat was first seen, if that is later.
            analysis_start = addition_date - datetime.timedelta(days =  compute_slope_over)
            
            if cntdata['earliest_date'] > analysis_start:
                analysis_start = cntdata['earliest_date']


            # specify latest date to be modelled
            mc = ModelCounts(
            counts = cntdata['counts'],
            denominators = cntdata['denominators'],
            pcas_int_id = pcas_int_id,
            pc_cat = pcas_obj.pc_cat,
            earliest_date = analysis_start,
            latest_date=addition_date)
            
            res = mc.fit_nb()
            pdm.store_pcas_model_output(res)
            
            res = mc.fit_poisson()
            pdm.store_pcas_model_output(res)

            bar.finish()
            
            features = pdm.feature_association(feature = 'pangolearn:lineage_1')
            destfile = os.path.join(outputdir, stem + "_" + addition_date.isoformat() + "_nb.assocs.txt")
            with open(destfile, "wt") as f:
                features.to_csv(f, date_format="YYYY-MM-DD", index=False, sep=",")
        
    # finished
    logging.info("Finished, terminating program.")
