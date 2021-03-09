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


"""

# imports
import os
import glob
import datetime
import pandas as pd
import shutil
import Bio
import logging
import logging.handlers
import argparse
import time
from collections import Counter
from fn4client import fn4Client
import sentry_sdk


if __name__ == '__main__':


   # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour4_server, a service for bacterial relatedness monitoring.
                                     

Example usage: 
============== 
# show command line options 
python updating_covid_load.py --help  

# load into server specified at url
python updating_covid_load.py "http://localhost:5023"

""")
    parser.add_argument('server_url', type=str, action='store', nargs='?',
                        help='the url (and port) of the findNeighbour4 server into which we should insert', default='')
    parser.add_argument('fastadir', type=str, action='store', nargs='?',
                        help='the directory in which fasta files will appear', default='')
    parser.add_argument('fileglob', type=str, action='store', nargs='?',
                        help='a pattern to glob for', default='')
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
    logdir = os.path.join(args.fastadir, 'logs')
    completedir = os.path.join(args.fastadir,'completed')
    for testdir in [logdir, completedir]:
        os.makedirs(testdir, exist_ok = True)

    # launch logger
    logger = logging.Logger('fn4_load_multifasta')
    logger.setLevel(logging.INFO)
    timenow = datetime.datetime.now().isoformat()
    logfile = os.path.join(logdir, "fn4_load_multifasta.log")
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
    formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
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
        logger.info("No error monitoring via sentry.  Set environment variable FN_SENTRY_URL to do so.")

    # instantiate client
    fn4c = fn4Client(args.server_url)      # expects operation on local host; pass baseurl if somewhere else.

    
    existing_guids = set(fn4c.guids())
    clustering_created = False
    logger.info("There are {0} existing guids".format(len(existing_guids)))

    # add the reference sequence as the root if not already present
    ref_guid = '--Wuhan-Reference--'
    ref_guid_present = fn4c.guid_exists(ref_guid)
    if not ref_guid_present:
        logger.info("Adding reference")
        for record in Bio.SeqIO.parse("../reference/nc_045512.fasta", 'fasta'):

            seq = str(record.seq).upper()
            res = fn4c.insert(guid=ref_guid,seq=seq)
    else:
        logger.info("Reference already present")

    for fastafile in glob.glob(os.path.join(args.fastadir, args.fileglob)):
        logger.info("Scanning {0}".format(fastafile))

        nSkipped = 0
        nBad = 0
        nGood = 0
        failed = []
        i = 0
        sr = None
        for record in Bio.SeqIO.parse(fastafile, 'fasta'):

            # build in pause if high storage ratio ('fragmentation')
            if nGood % 50 == 0:
                server_database_usage = fn4c.server_database_usage()
                sr = server_database_usage['latest_stats']['storage_ratio']
                # check whether database is keeping repacked adequately
                logger.info("Examined {0} / skipped {1}.  Database neighbour fragmentation is {2} (target: 1)".format(i,nSkipped,sr))
                while sr > 30:            #   ratio of records containing neighbours to those containing samples - target is 1:1
                    logger.info("Waiting 3 minutes to allow repacking operations.  Will resume when fragmentation, which is now {0}, is < 30.".format(sr))
                    time.sleep(180)    # restart in 3 mins  if below target
                    server_database_usage = fn4c.server_database_usage()
                    sr = server_database_usage['latest_stats']['storage_ratio']
                logger.info("Restarting after  pause. fragmentation is now {0} ".format(sr))  

            i = i + 1
            t1 = datetime.datetime.now()
            guid = record.id
            guid = guid.replace('/','-')
            guid = guid.replace(':','-')

            seq = str(record.seq).upper()
            seq = seq.replace('?','N')
            seq = seq.replace(' ','N')
            counter = Counter(list(seq))
            
            res = {'record_id':record.id,'guid':guid, 'seqlen':len(seq), **counter}
            if not guid in existing_guids:
                res = fn4c.insert(guid=guid,seq=seq)
                if not res.status_code == 200:
                    # failed to add
                    failed.append((i,guid))
                    msg = "** FAILED **"

                else:
                    msg = "Succeeeded"
                    nGood +=1

                t2 = datetime.datetime.now()
                i1 = t2-t1
                s = i1.total_seconds()
                logger.info("Scanned {0} Adding #{1} ({2}) {3} in {4} secs.".format(i, nGood, guid, msg, s, counter))
                
            else:
                nSkipped +=1

            if i % 2500 == 0:
                logger.info("Examined {0} / skipped {1}".format(i,nSkipped))


        logger.info("Complete.  Skipped {0} guids which already exist in the server.  There are {1} bad sequences".format(nSkipped, nBad))

        if len(failed)>0:
            # do not move the file
            logger.warning("findneighbour_load | not all fasta files could be uploaded; {0} failed".format(len(failed)))
            logger.info("Failed samples are: {0}".format(failed))
        else:
            logger.info("Sample load succeeded.  Moving fasta file to /completed directory")
            shutil.move(fastafile, completedir)
        
        # finished
        logging.info('Finished, terminating program.')