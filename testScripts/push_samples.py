import readline
import time
import json
import uuid
import datetime
import pdb
import gzip
import os.path
import xmlrpc.client
import socket
import sys
import os
import glob
import sys
from Bio import SeqIO

### push fasta files to findNeighbour
### python push_samples.py localhost R00000039

json_config_file=sys.argv[1]
with open(json_config_file, 'rt') as f:
	txt=f.read()
	print("Read config file: {0}".format(txt))
	CONFIG=json.loads(txt)

if not type(CONFIG)==dict:
	raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	
client=xmlrpc.client.ServerProxy("http://{0}:{1}".format(CONFIG['IP'],CONFIG['PORT']))

# note timings
startParse=datetime.datetime.now()
nTested=0
keyList=[]
# we are going to read the sequences from a test file provided by Trien.
testpath="/home/dwyllie/data/relatednesstest/TB_FASTA/*_v3.fasta" 
nRead=0
fastaFiles=glob.glob(testpath)
for fastaFile in fastaFiles:
	with open(fastaFile, 'rt') as f:
		for seq_record in SeqIO.parse(f, 'fasta'):
			guid=str(os.path.basename(fastaFile)[0:36])
			if True: 	# not client.exist_sample(guid):
				seq = str(seq_record.seq)
				print(datetime.datetime.now(), nTested, guid,len(seq),seq[0:25])
				nTested+=1
				if nTested % 50==0:
					print("{1} Loaded {0}".format(nTested, datetime.datetime.now()))
				result = client.insert(guid,seq)
			else:
				result = 'Already present, as determined by exist_sample'
			print(guid, result)

# note  the end of parse time.      
endParse=datetime.datetime.now()
#print("LOADED to ",this_persistenceStore, startParse, endParse, nTested, "Time per sample=",(endParse-startParse)/nTested)
try:
  perEvent=(endParse-startParse)/nTested
except ZeroDivisionError:
  print("No sequences found")
  exit()
  
print("Loading speed per sample was (seconds) ",perEvent.total_seconds())

