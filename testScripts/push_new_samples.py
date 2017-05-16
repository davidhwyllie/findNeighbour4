""" push new fasta files to EW

Conducts a scan of a directory for fasta files, and pushes those that don't exist into the EW server

THIS ASSUMES THERE IS A CONFIG FILE CONTAINING A DICTIONARY CONFIG, with keys 'IP' AND 'PORT' PASSED AS sys.argv[1]
i.e. the first and only argument.

An example of such a file would read:
{
"DESCRIPTION":"A test server operating in ../unittest_tmp, only suitable for testing",
"PORT":8184,
"GLOBPATH":"/home/dwyllie/data/relatednesstest/TB_FASTA/*_v3.fasta" }
"""
import logging
import readline
import time
import json
import datetime
import gzip
import os.path
import xmlrpc.client
import socket
import sys
import os
import glob
import sys
from Bio import SeqIO

# set up logging
# if not within the COMPASS framework, need to decide where to log to
#logging.getLogger()

# check input
json_config_file=sys.argv[1]
with open(json_config_file, 'rt') as f:
	txt=f.read()
	logging.info("Read config file: {0}".format(txt))
	CONFIG=json.loads(txt)

if not type(CONFIG)==dict:
	raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	
if not set(CONFIG.keys()) == set(['DESCRIPTION','IP','PORT','GLOBPATH']):
	raise KeyError("The dictionary must have four components: DESCRIPTION, IP, PORT, and GLOBPATH.  The latter is passed to glob.glob to find the files.  The dictionary actually looks like this: {0}".format(CONFIG))

# try to start the client.  Will throw and error if it cannot connect.
# should wrap with try/catch to log to logfile, if any
logging.info("Trying to make server connection ...")
try:
	client=xmlrpc.client.ServerProxy("http://{0}:{1}".format(CONFIG['IP'],CONFIG['PORT']))
except Exception as e:
	logging.exception(e)
	raise e


logging.info("Connected, checking existing guids ...")
guidlist = json.loads(client.get_all_guids())
guids = set()
for item in guidlist:
	guids.add(item['guid'])
nTested=0
nAdded=0

logging.info("Connected, checking existing guids vs. those found using the glob pattern  ...")
fastaFiles = glob.glob(CONFIG['GLOBPATH'])
for fastaFile in fastaFiles:
	with open(fastaFile, 'rt') as f:
		for seq_record in SeqIO.parse(f, 'fasta'):
			guid = str(os.path.basename(fastaFile)[0:36])
			nTested += 1
			if nTested % 250==0:
				logging.info("{1} checked {0}".format(nTested, datetime.datetime.now()))
			if not guid in guids:				
				if not client.exist_sample(guid):
					seq = str(seq_record.seq)
					result = client.insert(guid,seq)
					nAdded += 1
		#if nTested>500:
		#	break

logging.info("Scan finished.  Added {0}".format(nAdded))
