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
import unittest
import logging
import urllib
import requests
import base64
import hmac
import hashlib
import pandas as pd
########################################################


json_config_file=sys.argv[1]
with open(json_config_file, 'rt') as f:
	txt=f.read()
	print("Read config file: {0}".format(txt))
	CONFIG=json.loads(txt)

if not type(CONFIG)==dict:
	raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	
client=xmlrpc.client.ServerProxy("http://{0}:{1}".format(CONFIG['IP'],CONFIG['PORT']))

print("################### query_get_all_guids from ew2 ###########################")
nTested=1
startParse=datetime.datetime.now()

#retval=client.get_all_values()
guids=set()
retval=client.get_all_guids()
retlist=json.loads(retval)
for d in retlist:
	guids.add(d['guid'])
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("get guids",perEvent.total_seconds())


#################################################################################
# compare with EW
print("Comparing with EW")
print("Reading EW2 edges")
inputfile=os.path.join("..","ew1cmp",'ew2_edges.txt')
edges2=pd.read_csv(inputfile,  sep='\t', header=0)
edges2.columns=['guid1','guid2','dist_2','n1','n2','n3']

print("Reading bugmat results")
inputfile=os.path.join("..","ew1cmp","bugmat_snp.txt")
edges1=pd.read_csv(inputfile,  sep='\t', header=None)
edges1.columns=['guid1','guid2','dist']

print("Scanning bugmat file for those guids which should be present in ew2:")
edges1['expectResult']= False
for i in range(len(edges1.index)):
	if edges1.iloc[i,0] in guids and edges1.iloc[i,1] in guids:
		edges1.iloc[i,3] = True
		if i % 10000 ==0:
			print("Checking row {0}".format(i))
edges3=edges1.merge(edges2,
					how='left',
					left_on = ['guid1','guid2'],
					right_on = ['guid1','guid2'])

outputfile=os.path.join("..","ew1cmp","one_vs_two.txt")
edges3.to_csv(outputfile, sep='\t')
