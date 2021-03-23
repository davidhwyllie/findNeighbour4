""" builds a guide tree - a phylogenetic tree from a representative sample of sequences,
obtained by SNP based sampling.

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

import os
import glob
import datetime
import json
import pandas as pd
import datetime
import argparse
import logging
import logging.handlers
from logging.config import dictConfig

import pathlib
from random import sample as random_sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from collections import Counter
from fn4client import fn4Client
from read_config import ReadConfig

def neighbours_of(sample, snp_cutoff=3, max_neighbourhood_mixpore = 50):
    """ obtains a list of neighbouring samples of *sample* , and an indication of whether *sample* is mixed.

    Parameters
    sample      the identifier of the sample to start with
    cutoff      how many snps from sample neighbours should be identified.  Recommend 1 or 2 for COVID-19
    max_neighbourhood_mixpore:  the maximum number of similar samples (defined at being less than or equal to the cutoff) used for mixpore computations (assessing mixtures of similar strains).  If None, no mixpore computation is done.

    Returns
    dictionary with keys
        mixpore_test_p  : p value indicating whether number of N/M in sites of recent change (differing between the close neighbours of sample) are more common than in the rest of the genome
                            if significant, this is indicative of a mixture of samples, at least in TB
        neighbours:       a list of neighbours of sample
    """

    res = fn4c.guid2neighbours(sample, threshold = snp_cutoff)        # find neighbours within a snv cutoff
    neighbours = []
    for related_sample, distance in res:
            neighbours.append(related_sample)

    if max_neighbourhood_mixpore is None:
        return {'mixpore_test_p': 1, 'neighbours':neighbours}

    # if there are more than 50 samples, randomly downsample.
    if len(neighbours)>50:
        for_msa = random_sample(neighbours, max_neighbourhood_mixpore)
    else:
        for_msa = neighbours
    for_msa.append(sample)          # build the samples into the MSA

    # to just get the MSA
    msa_df = fn4c.msa(for_msa, output_format='json-records', what='N_or_M')
    msa_df.set_index('guid', inplace=True, drop=True) 
    mixpore_test_p = msa_df.at[sample, 'p_value3']
    return {'mixpore_test_p': mixpore_test_p, 'neighbours':neighbours}

if __name__ == '__main__':

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_guidetreeg, a findNeighbour4 component.
                                     

Example usage: 
============== 

## does not require findNeighbour4_server to be running
Minimal example:
python findNeighbour4_guidetree.py ../config/myConfigFile.json  

if a config file is not provided, it will run (as does findNeighbour4_server) is debug mode: it will run once, and then terminate.  This is useful for unit testing.  If a config file is specified, the clustering will  run until terminated.  

Checks for new sequences are conducted once per minute.



""")
    parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                        help='the path to the configuration file', default=''  )
    parser.add_argument('--select_from_pca_output_file', help='selects samples from an sqlite database, which contains details of good quality samples.  If missing, will select based on --sample_quality_threshold', action='store')
    parser.add_argument('--sample_quality_cutoff', help='Only include in the tree samples with >= this proportion of ACTG bases', action='store')
    parser.add_argument('--snp_cutoff', help='Find representatives for tree building > snp_cutoff apart.  If omitted, maximum distance stored in the relatedness server is used', action='store')
    parser.add_argument('--max_neighbourhood_mixpore', help='exclude samples failing mixpore test on at most max_neighbourhood_mixpore closely related samples.  If omitted, mixpore check is not done', action='store')
    args = parser.parse_args()
    
    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 guidetree .. reading configuration file.")

    if len(args.path_to_config_file)>0:
            configFile = args.path_to_config_file
            debugmode = False
            logging.info(configFile)
    else:
            raise FileNotFoundError("No config file name supplied ; this is required ")

    rc = ReadConfig()
    CONFIG = rc.read_config(configFile)

    ########################### SET UP LOGGING #####################################  
    # create a log file if it does not exist.
    print("Starting logging")
    logdir = os.path.dirname(CONFIG['LOGFILE'])
    pathlib.Path(os.path.dirname(CONFIG['LOGFILE'])).mkdir(parents=True, exist_ok=True)

    # set up logger
    logger = logging.getLogger()
    loglevel=logging.INFO
    if 'LOGLEVEL' in CONFIG.keys():
            if CONFIG['LOGLEVEL']=='WARN':
                    loglevel=logging.WARN
            elif CONFIG['LOGLEVEL']=='DEBUG':
                    loglevel=logging.DEBUG

    # configure logging object 
    logger.setLevel(loglevel)       
    logfile = os.path.join(logdir, "clustering-{0}".format(os.path.basename(CONFIG['LOGFILE'])))
    print("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
    formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
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
    if 'SENTRY_URL' in CONFIG.keys():
            logger.info("DISABLED Launching communication with Sentry bug-tracking service")
            #sentry_sdk.init(CONFIG['SENTRY_URL'], integrations=[FlaskIntegration()])

    ########################### read sqlite file containing samples to select, if it exists ##########
    relabel = False  
    if args.select_from_pca_output_file is not None:
        print("SQLITE FILE {0} provided".format(args.select_from_pca_output_file))
    else:
        logging.info("No SQLite file provided")
    

    exit()

    # instantiate client
    fn4c = fn4Client("http://localhost:{0}".format(CONFIG['REST_PORT']))      # expects operation on local host; pass baseurl if somewhere else.

    snp_cutoff = 3                                 # how many snp from a sample should we include neighbours
    sample_quality_cutoff = 0.9                    # only build tree from samples with >= 90% ACTG
    max_neighbourhood_mixpore = None               # set to integer, e.g. 50 to do mixpore check on all samples
    
    print("Loading all samples & annotations.")
    sample_annotations = fn4c.annotations()
    all_samples = sample_annotations.index.to_list()
    high_quality_samples = sample_annotations.loc[sample_annotations['DNAQuality:propACTG']>sample_quality_cutoff].copy()
    
    sampling_population = high_quality_samples.index.to_list()

    all_samples = set(sampling_population)
    sampling_population = set(sampling_population)          # samples now (note: more may be added during computations by other processes, so we pick a set to work on at the start of the computations)

    print("From all {0} eligible samples currently present, studying {1} high quality (ACTG > 0.9) samples ".format(len( sample_annotations.index.to_list()), len(sampling_population)))
    representative_subsample = set()
    remaining_samples = sampling_population.copy()
    representative_size = {}

    ## DEBUG 
    # remaining_samples = set(random_sample(list(remaining_samples), 5) )

    mixed_samples = set()
    singletons = set()
    iteration = 0
    total_delta = 0
    while len(remaining_samples)>0:

        iteration +=1
        test_sample = random_sample(list(remaining_samples),1)[0]        # randomly sample one
        test_sample_neighbour_info = neighbours_of(test_sample, snp_cutoff, max_neighbourhood_mixpore=max_neighbourhood_mixpore)        
        
        # filter these neighbours.   Only include neighbours which are part of all_samples.
        # because samples get added to the server all the time, this is not guaranteed.
        test_sample_neighbours = []
        for item in test_sample_neighbour_info['neighbours']:
            if item in all_samples:
                test_sample_neighbours.append(item)

        # maybe put in place a quality filter (total Ns)
        delta = 1           # one less at a minimum (test_sample is never reconsidered);

        if len(test_sample_neighbours) == 0:
            singletons.add(test_sample)

        elif test_sample_neighbour_info['mixpore_test_p'] is None:
            # cannot tell, not enough neighbours
            # include it
            representative_subsample.add(test_sample)

        elif test_sample_neighbour_info['mixpore_test_p'] > 1e-5:         # not obviously mixed
            representative_subsample.add(test_sample)
            representative_size[test_sample] = len(test_sample_neighbours)

            # remove sample and neighbours from further consideration
            n1 = len(remaining_samples)
            for item in test_sample_neighbours:
                remaining_samples.discard(item)         # discard as item may be low quality or otherwise pre-excluded
            n2 = len(remaining_samples)
            delta = 1+ n1- n2
            total_delta = total_delta + delta   # total number of samples

        elif test_sample_neighbour_info['mixpore_test_p'] <= 0.05:
            mixed_samples.add(test_sample)

        else:
            raise ValueError("Unhandled situation")

        # remove the test sample, if it has not already been removed
        remaining_samples.discard(test_sample)
        if len(representative_subsample)>0:
            ratio_selected = total_delta/float(len(representative_subsample))
        else: 
            ratio_selected = None
        print("{0} | Iteration: {1} | Remaining samples {2} |  Change now/total {9}/{10} | Sample = {3} |  Mixture statistic {4} | Neighbours {5} | Selected representatives {6} | Mixed samples {7} | Singletons {8} | Approx Subsample 1 in {11:.0f}".format(
          datetime.datetime.now().isoformat(),
          iteration,
          len(remaining_samples),
          test_sample,
          test_sample_neighbour_info['mixpore_test_p'],
          len(test_sample_neighbours),
          len(representative_subsample),
          len(mixed_samples),
          len(singletons),
          delta, 
          total_delta,
          ratio_selected
          
        ))

        #if iteration > 1:      # DEBUG
        #    break

    print("Generating fasta output")
    # export masked fasta (note: holds all selected sequences in ram)
    fasta_outputfile =  'tree_selected.fasta' 
    seqs = []
    for test_sample in representative_subsample:
        seq = fn4c.sequence(test_sample)
        sr = SeqRecord(     Seq(seq['masked_dna']), 
                                id= test_sample,
                                description="| within {0} snp of {1} others".format(snp_cutoff, representative_size[test_sample]) 
                                )        
        seqs.append(sr)

    with open(fasta_outputfile, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
    print("Selected samples are written to {0}".format(fasta_outputfile))
    print("Finished")
    
