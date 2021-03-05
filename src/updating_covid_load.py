""" loads COVID-19 samples into a findNeighbour4 server
from a fasta file (e.g. the COVID-19 alignment produced by COG-UK)

assumes a findNeighbour4 server is running

Paths expected are currently hard-coded

Example usage: 
============== 
# show command line options 
python updating_covid_load.py --help  

python3 updating_covid_load.py [server_url] [directory to look for fastas in] [filename to glob within directory]
# example usage
pipenv run python3 updating_covid_load.py http://localhost:5023 /srv/data/covid COVID_MSA*.fasta


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
    logger = logging.Logger('updating_covid_load')
    logger.setLevel(logging.INFO)
    timenow = datetime.datetime.now().isoformat()
    logfile = os.path.join(logdir, "updating_covid_load.log")
    file_handler = logging.handlers.RotatingFileHandler(logfile, mode = 'a', maxBytes = 1e7, backupCount = 7)
    formatter = logging.Formatter( "%(asctime)s | %(pathname)s:%(lineno)d | %(funcName)s | %(levelname)s | %(message)s ")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.info("Startup with arguments: {0}".format(args))

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
    logging.info("There are {0} existing guids".format(len(existing_guids)))

    # add the reference sequence as the root if not already present
    ref_guid = '--Wuhan-Reference--'
    ref_guid_present = fn4c.guid_exists(ref_guid)
    if not ref_guid_present:
        logging.info("Adding reference")
        for record in Bio.SeqIO.parse("../reference/nc_045512.fasta", 'fasta'):

            seq = str(record.seq).upper()
            res = fn4c.insert(guid=ref_guid,seq=seq)
    else:
        logging.info("Reference already present")

    for fastafile in glob.glob(os.path.join(args.fastadir, args.fileglob)):
        logging.info("Scanning {0}".format(fastafile))

        nSkipped = 0
        nBad = 0
        nGood = 0
        failed = []
        i = 0
        for record in Bio.SeqIO.parse(fastafile, 'fasta'):
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
                logging.info("Scanned {0} Added {1} Sample: {2} {3} in {4} secs.  Composition: {5}".format(i, nGood, guid, msg, s, counter))
                
            else:
                nSkipped +=1

            if i % 1000 == 0:
                logger.info("Examined {0} / skipped {1}".format(i,nSkipped))   
     
        logging.info("Complete.  Skipped {0} guids which already exist in the server.  There are {1} bad sequences".format(nSkipped, nBad))

        logging.info("Failed samples are: {0}".format(failed))
        
        if len(failed)>0:
            # do not move the file
            pass
