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

### TODO REFACTOR INTO A UNITTESTED CLASS ###

"""

# import libraries
import os
import pathlib
import logging
import sqlalchemy
import logging.handlers
import datetime
import pandas as pd
import sentry_sdk
import argparse
import progressbar
import subprocess
from io import StringIO
from random import sample as random_sample
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo

# config
from findn.common_utils import ConfigManager

# startup
from findn.mongoStore import fn3persistence
from findn.hybridComparer import hybridComparer
from tree.manipulate_tree import ManipulateTree

if __name__ == '__main__':

    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_guidetreeg, a findNeighbour4 component.
                                     

Example usage: 
============== 

## does not require findNeighbour4_server to be running
Minimal example:
python findNeighbour4_guidetree.py ../config/myconfig_file.json  

Realistic example for testing
pipenv run python3 findNeighbour4_guidetree.py ../demos/covid/covid_config_v3.json --select_from_pca_output_file /backup/bricestudy/pca20210323/pca_output.sqlite --output_filestem guidetree --outputdir /backup/bricestudy --debug
pipenv run python3 findNeighbour4_guidetree.py ../demos/covid/covid_config_v3.json --select_from_pca_output_file /backup/bricestudy/pca20210323/pca_output.sqlite --output_filestem guidetree --outputdir /backup/bricestudy --root_sample="--Wuhan-Reference--" --debug

Realistic example for running
pipenv run python3 findNeighbour4_guidetree.py ../demos/covid/covid_config_v3.json --select_from_pca_output_file /backup/bricestudy/pca20210323/pca_output.sqlite --output_filestem guidetree --outputdir /backup/bricestudy --root_sample="--Wuhan-Reference--" 

if a config file is not provided, it will run (as does findNeighbour4_server) is debug mode: it will run once, and then terminate.  This is useful for unit testing.  If a config file is specified, the clustering will  run until terminated.  

Checks for new sequences are conducted once per minute.

""")
    parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                        help='the path to the configuration file', default=''  )
    parser.add_argument('--outputdir', type=str, action='store', nargs='?',
                        help='the output directory.  Will try to make the directory if it does not exist'  )
    parser.add_argument('--select_from_pca_output_file', help='selects samples from an sqlite database, which contains details of good quality samples.  Will write a guidetree table into the database containing samples selected. If missing, will select based on --sample_quality_threshold', action='store')
    parser.add_argument('--sample_quality_cutoff', help='Only include in the tree samples with >= this proportion of ACTG bases', action='store')
    parser.add_argument('--snp_cutoff', help='Find representatives for tree building > snp_cutoff apart.  If omitted, maximum distance stored in the relatedness server is used', action='store')
    parser.add_argument('--max_neighbourhood_mixpore', help='exclude samples failing mixpore test on at most max_neighbourhood_mixpore closely related samples.  If omitted, mixpore check is not done', action='store')
    parser.add_argument('--root_sample', help='regard the closest isolate to this sample as the root.  The root_sample can be ancestral (e.g. Wuhan reference, or the inferred TB Ancester (I Comas 2013)), or it can be an outgroup.  A tree will be generated without this sample.', action='store')
    parser.add_argument('--output_filestem', help='the first part of the output file basename.  e.g. if you want the output file called my_tree.nwk, enter my_tree', action='store')
    parser.add_argument('--debug', help='only analyse the first 20 samples', action='store_true')
    args = parser.parse_args()
    
    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 guidetree .. reading configuration file.")

    if len(args.path_to_config_file)>0:
            config_file = args.path_to_config_file
            debugmode = False
            logging.info(config_file)
    else:
            raise FileNotFoundError("No config file name supplied ; this is required ")

    cfm = ConfigManager(config_file)  
    CONFIG = cfm.read_config()

    ########################### SET UP LOGGING & BUG TRACKING #################################  
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
    logfile = os.path.join(logdir, "guidetree_{0}".format(os.path.basename(CONFIG['LOGFILE'])))
    print("Logging to {0} with rotation".format(logfile))
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
    formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # launch sentry if API key provided
    if 'SENTRY_URL' in CONFIG.keys():
            logger.info("Launching communication with Sentry bug-tracking service")
            sentry_sdk.init(CONFIG['SENTRY_URL'])

    ##################################################################################################
     # open PERSIST and hybridComparer object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed

    logger.info("Connecting to backend data store")
    try:
            PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
                connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                debug=0
                                   )  # if in debug mode wipes all data.  This is not what is wanted here, even if we are using unittesting database

    except Exception:
            logger.exception("Error raised on creating persistence object")
            raise

    ################################# object to do MSA if required #############################################
    # open PERSIST and hybridComparer object used by all samples
    # this is only used for data access and msa.
    # inserts are not allowed
    logger.info("Building hybridComparer object")
    hc = hybridComparer(reference=CONFIG['reference'],
        maxNs=CONFIG['MAXN_STORAGE'],
        snpCeiling=  CONFIG['SNPCEILING'],
        excludePositions=CONFIG['excluded'],
        preComparer_parameters={},
        PERSIST=PERSIST,
        disable_insertion = True)
    


    logging.info("Loading all samples & annotations.")
    sample_annotations = pd.DataFrame.from_dict(PERSIST.guid2items(None,None), orient='index')
 
    server_samples = sample_annotations.index.to_list()
    
    ########################### read sqlite file containing samples to select, if it exists ##########
    if args.select_from_pca_output_file is not None:
        if not os.path.exists(args.select_from_pca_output_file):
            raise FileNotFoundError(args.select_from_pca_output_file)
        if not args.select_from_pca_output_file.endswith('.sqlite'):
            raise FileNotFoundError("The file {0} exists but it does not end with .sqlite".format(args.select_from_pca_output_file))

        # Read sqlite query results into a pandas DataFrame
        engine = sqlalchemy.create_engine("sqlite:///{0}".format(args.select_from_pca_output_file), echo=True)
        conn =  engine.connect()

        df = pd.read_sql_query("SELECT * from built_with_guids", conn)
        all_samples = df['built_with_guids'].to_list()

        extra_samples_in_pca_only = set(all_samples) - set(server_samples)
        if len(extra_samples_in_pca_only)>0:
            raise KeyError("The PCA sqlite file is incompatible with the server data; the PCA contains {0} samples  not present in the server".froamt9len(extra_samples_in_pca_only))

        logging.info("Recovered {0} samples from PCA result in {1}".format(len(all_samples), args.select_from_pca_output_file))
    else:
        logging.info("No SQLite file containing PCA output provided.  The tree will be built from samples in the relatedness server.")
        all_samples = server_samples.copy()

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
    
    ########################### set distance ##########
    if args.snp_cutoff is not None:
        if args.snp_cutoff >=0 and args.snp_cutoff<=CONFIG["SNPCEILING"]:
            # it is valid
            snp_cutoff = args.snp_cutoff 
        else:
            raise ValueError("Invalid snp_cutoff: must be either omitted; or 0 or more; or less than {0}.  It was {1}".format(CONFIG['SNPCEILING'], args.snp_cutoff))
    else:
        snp_cutoff = CONFIG["SNPCEILING"]
    logging.info("SNP cutoff set to {0}".format(snp_cutoff))

    ########################### sample quality cutoff ##########
    if args.sample_quality_cutoff is not None:
        if args.sample_quality_cutoff >=0 and args.sample_quality_cutoff <=1:
            # it is valid
            sample_quality_cutoff =  args.snp_cutoff 
        else:
            raise ValueError("Invalid snp_cutoff: must be either omitted; or 0 to 1 inclusive.  It was {1}".format(CONFIG['SNPCEILING'], args.sample_quality_cutoff))
    else:
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

    ########################### decide whether to force inclusion of the root sample #######################
    if args.root_sample is not None:
        if PERSIST.guid_exists(args.root_sample) :
            # it is valid
            root_sample = args.root_sample
            need_to_add_root_sample = True
        else:
            raise ValueError("Root sample does not exist {0}".format(args.root_sample))
    else:
        root_sample = None    # default value
        need_to_add_root_sample = False     # we don't need to add this
    logging.info("Root sample set to {0}".format(root_sample))

    ####################################  START ################################################################
    ##                      TO CONVERT TO CLASS      ###########################################################
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

