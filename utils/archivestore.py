""" classes to build and store collections of reference compressed sequences
suitable for testing the performance of pca.

Includes an archiveStore object, which offers a subset of methods present in mongoStore 
(a class accessing the mongodb database) and can be used as a drop in replacement for 
mongoStore when 

- unittesting testing PCA generation 
- benchmarking PCA
- testing the detection of new variants by PCA

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
import sys
import pathlib
import json
import logging
import sqlalchemy
import logging.handlers
import warnings
import datetime
import glob
import sys
import pandas as pd
import numpy as np
import sentry_sdk
import argparse
import progressbar
import time
import subprocess
from io import StringIO
from random import sample as random_sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo

from sentry_sdk import capture_message, capture_exception
from sentry_sdk.integrations.flask import FlaskIntegration

# logging
from logging.config import dictConfig

# startup
from findn.mongoStore import fn3persistence
from findn.hybridComparer import hybridComparer
from findn.read_config import ReadConfig
from trees.manipulate_tree import ManipulateTree

if __name__ == '__main__':


    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""make_referencecompressed_archive, a findNeighbour4 component.
                                     

Example usage: 
============== 

## does not require findNeighbour4_server to be running

python make_referencecompressed_archive.py ../config/myConfigFile.json  

if a config file is not provided, it will not run.

PERSIST.refcompressedsequence_guids(
        refcompressedsequence_read
        config   


""")
    parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                        help='the path to the configuration file', default=''  )
    parser.add_argument('--outputdir', type=str, action='store', nargs='?',
                        help='the output directory.  Will try to make the directory if it does not exist'  )
    parser.add_argument('--max_samples', help='the maximum number of samples to export.  If None, all samples are exported', action='store')
    parser.add_argument('--write_fasta', help='write the masked sequences to fasta', action='store_true')
    args = parser.parse_args()
    
    ############################ LOAD CONFIG ######################################
    print("make_referencecompressed_archive .. reading configuration file.")

    if len(args.path_to_config_file)>0:
            configFile = args.path_to_config_file
            debugmode = False
            logging.info(configFile)
    else:
            raise FileNotFoundError("No config file name supplied ; this is required ")

    rc = ReadConfig()
    CONFIG = rc.read_config(configFile)

    ########################### SET UP LOGGING #####################################  
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
    logfile = os.path.join(logdir, "guidetree_{0}".format(os.path.basename(CONFIG['LOGFILE'])))
    print("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
    formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    ##################################################################################################
     # open PERSIST object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed

    logger.info("Connecting to backend data store")
    try:
            PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
                connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                debug=0
                                   )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception as e:
            logger.exception("Error raised on creating persistence object")
            raise

    logging.info("Loading all samples .")
    guids = PERSIST.guids()
 
    print(guids)
    exit()

    ########################### create output filename if it does not exist ##########
    if args.output_filestem is not None:
        filestem = args.output_filestem
    else:
        filestem = "guidetree-{0}".format(datetime.datetime.now().isoformat().replace(':','-').replace('.','-'))

    ########################### create output directory if it does not exist ##########
    if args.outputdir is not None:
        outputdir = args.outputdir
        logging.info("Writing output to {0}".format(outputdir))
    else:
        outputdir = "/tmp/{0}".format(filestem)
        logging.info("Writing output to temporary directory at {0}.  Set --outputdir to change this".format(outputdir))
    if not os.path.exists(outputdir):
            os.makedirs(outputdir)
    
    ########################### sample quality cutoff ##########
    sample_quality_cutoff = 0.9     # default value
    logging.info("Quality cutoff set to {0}".format(sample_quality_cutoff))

    ########################### decide whether to do a mixpore (base composition) check  #######################
    if args.max_neighbourhood_mixpore is not None:
        if args.max_neighbourhood_mixpore >0 :
            # it is valid
            max_neighbourhood_mixpore =  args.max_neighbourhood_mixpore
        else:
            raise ValueError("Invalid max_neighbourhood_mixpore: must be >0 It was {0}".format(args.max_neighbourhood_mixpore))
    else:
        max_neighbourhood_mixpore = None    # default value
    logging.info("Max neighbourhood mixpore set to {0}".format(max_neighbourhood_mixpore))


    ###################### select high quality samples from which to build the tree ############################
    high_quality_samples = sample_annotations.loc[sample_annotations['DNAQuality:propACTG']>sample_quality_cutoff].copy()
    high_quality_samples = high_quality_samples[high_quality_samples.index.isin(all_samples)]   # subselect those within the PCA if appropriate
    sampling_population = high_quality_samples.index.to_list()
    sampling_population = set(sampling_population)          # need rapid search
    logging.info("From all {0} eligible samples currently present, studying {1} high quality (ACTG > 0.9) samples ".format(len( sample_annotations.index.to_list()), len(sampling_population)))
    representative_subsample = set()
    remaining_samples = sampling_population.copy()
    representative_size = {}

    ################################### DEBUG : ONLY CONSIDER 20 for testing #####################################
    if args.debug:
        remaining_samples = set(random_sample(list(remaining_samples), 20) )
        logging.warning("Running in debug mode.  Only 20 samples will be considered")

    #################################  prepare to iterate #########################################################
    all_samples = set(all_samples)
    mixed_samples = set()
    singletons = set()
    iteration = 0
    total_delta = 0

    representatives = []

    total_n_samples  = len(remaining_samples)
    bar = progressbar.ProgressBar(max_value=len(remaining_samples))
    
    # sample high quality sequences, and their neighbours, until none remain.
    while len(remaining_samples)>0:

        iteration +=1

        if need_to_add_root_sample:
            test_sample = root_sample
            need_to_add_root_sample = False     # added it now
        else:
            test_sample = random_sample(list(remaining_samples),1)[0]        # randomly sample one

        test_sample_neighbour_info = PERSIST.guid2neighbours(test_sample, snp_cutoff, returned_format=1)      #  [guid, distance]
        
        # filter these neighbours.   Only include neighbours which are part of all_samples.
        # because samples get added to the server all the time, this is not guaranteed.
        test_sample_neighbours = []
        for item in test_sample_neighbour_info['neighbours']:
            if item[0] in all_samples:
                test_sample_neighbours.append(item)

        # maybe put in place a quality filter (total Ns)
        delta = 1           # one less at a minimum (test_sample is never reconsidered);

        if len(test_sample_neighbours) == 0:
            singletons.add(test_sample)
            is_singleton = True
        else: 
            representative_subsample.add(test_sample)
            representative_size[test_sample] = len(test_sample_neighbours)
            is_singleton= False
            
            # remove sample and neighbours from further consideration
            n1 = len(remaining_samples)
            for item in test_sample_neighbours:
                remaining_samples.discard(item[0])         # discard as item may be low quality or otherwise pre-excluded
            n2 = len(remaining_samples)
            delta = 1+ n1- n2
            total_delta = total_delta + delta   # total number of samples

        # store the neighbours in a list for export
        for neighbour,dist in test_sample_neighbours:
            representatives.append({'representative':test_sample, 'neighbour':neighbour, 'distance':dist, 'is_singleton':is_singleton})

        # remove the test sample, if it has not already been removed
        remaining_samples.discard(test_sample)
        if len(representative_subsample)>0:
            ratio_selected = total_delta/float(len(representative_subsample))
        else: 
            ratio_selected = None
        bar.update(total_n_samples-len(remaining_samples))

    ########################### write representatives to sqlite file, if it exists ##########
    reps = pd.DataFrame.from_records(representatives)
    if args.select_from_pca_output_file is not None:

        reps.to_sql('guidetree_reps', conn, if_exists='replace')
        logging.info("Wrote guidetree members to sqlite")
        conn.close()


    ########################### write representatives to fasta file ##########
    logging.info("Generating fasta output from {0} samples; ignoring {1} singletons.".format(len(representative_subsample), len(singletons)))
    # export masked fasta (note: holds all selected sequences in ram)

    fasta_outputfile =  os.path.join(outputdir, '{0}.fasta'.format(filestem))
    csv_outputfile = os.path.join(outputdir, '{0}.csv'.format(filestem))
    nwk_outputfile = os.path.join(outputdir, '{0}.nwk'.format(filestem))
    nwkr_outputfile = os.path.join(outputdir, '{0}.rooted.nwk'.format(filestem))
    meta_outputfile = os.path.join(outputdir, '{0}.treeinfo.txt'.format(filestem))
    reps.to_csv(csv_outputfile)

    seqs = []
    for test_sample in representative_subsample:
        seq = hc.uncompress_guid(test_sample)
        seq = seq.replace('N','-')  ## fastTree convention
        sr = SeqRecord(     Seq(seq), 
                                id= test_sample,
                                description="| within {0} snp of {1} others".format(snp_cutoff, representative_size[test_sample]) 
                                )        
        seqs.append(sr)

    with open(fasta_outputfile, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
    logging.info("Selected samples are written to {0}".format(fasta_outputfile))
    logging.info("Selected representatives are written to {0}".format(csv_outputfile))

    ########################## run fasttree #################################

    if os.environ.get("FASTTREE_DIR") is not None:
        single_cmd = "{0}/FastTreeMP -fastest -gtr -nopr -nt {1} > {2} 2> {3}".format(os.environ.get("FASTTREE_DIR"), fasta_outputfile, nwk_outputfile, meta_outputfile)
        subprocess_cmd = [
            "{0}/FastTreeMP".format(os.environ.get("FASTTREE_DIR")),
        "-fastest",
        "-gtr",
        "-nopr",
        "-nt",
        fasta_outputfile]
        logging.info("Running fasttree command equivalent to \n {0}".format(single_cmd))



        ## caution: if one changes the .env file FASTTREE_DIR environment variable to some other software called FastTreeMP, it will get run.  security risk ? how to mitigate
        res = subprocess.run(subprocess_cmd, text=True, capture_output=True, check=True)
        
        with open(nwk_outputfile, 'wt') as f:
            f.write(res.stdout)
        with open(meta_outputfile, 'wt') as f:
            f.write(res.stderr)

        logging.info("Newick as produced by fastTree written to {0}".format(nwk_outputfile))
        logging.info("Metadata written to {0}".format(meta_outputfile))

        # build rooted tree, if root is provided
        if root_sample is not None:
            logging.info("Producing a tree rooted using {0}".format(root_sample))
            tree = Phylo.read(StringIO(res.stdout), 'newick')
            tree.root_with_outgroup({"name":root_sample})
            tree.collapse(root_sample)
            with open(nwkr_outputfile,'w') as f:
                Phylo.write(tree, f, 'newick')
            logging.info("Rooted tree written to {0}".format(nwkr_outputfile))
            

    else:
        logging.error("No tree built as FASTTREE environment variable not found.  put it in .env if using a virtual environment, and point it to the directory containing the fasttreeMP executable.  The fasttreeMP should be complied with double precision flags, see fastTree docs.")
   
    ## option: downsample e.g. pipenv run python3 Treemmer_v0.3.py /backup/bricestudy/pca20210323/guidetree.nwk  -X 5000 2500 1000 (slow) or randomly (faaster)
    ## TODO: store the guidetree in the database, make accessible via front end. ?

