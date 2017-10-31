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
import resource
import psutil
import csv


json_config_file=sys.argv[1]
with open(json_config_file, 'rt') as f:
	txt=f.read()
	print("Read config file: {0}".format(txt))
	CONFIG=json.loads(txt)
 
if not type(CONFIG)==dict:
	raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	
client=xmlrpc.client.ServerProxy("http://{0}:{1}".format(CONFIG['IP'],CONFIG['PORT']))

print("################### query_get_all_values (edgelist)#########################")

retval=client.get_all_values()
edges=json.loads(retval)
nTested=0
outputfile="edges.txt"
with open(outputfile,'wt') as f:
	writer=csv.writer(f, dialect='excel-tab')
	writer.writerow(['guid1','guid2','dist','N1','N2','Nboth'])
	for edge in edges:
		writer.writerow(edge)
		nTested+=1
print("Write finished. {0} rows.".format(nTested))
exit()



