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

### push fasta files to findNeighbour

json_config_file=sys.argv[1]
with open(json_config_file, 'rt') as f:
	txt=f.read()
	print("Read config file: {0}".format(txt))
	CONFIG=json.loads(txt)

if not type(CONFIG)==dict:
	raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	
client=xmlrpc.client.ServerProxy("http://{0}:{1}".format(CONFIG['IP'],CONFIG['PORT']))


nTested=0
keyList=[]
print("Reading guid names from directory")
# we are going to read the sequences from a test file provided by Trien.
testpath="/home/dwyllie/data/relatednesstest/TB_FASTA/*_v3.fasta" 
nRead=0
guids=set()
fastaFiles=glob.glob(testpath)
for fastaFile in fastaFiles:
	nRead+=1
	guid=str(os.path.basename(fastaFile)[0:36])
	guids.add(guid)
	if nRead>200:
		break
print("Loaded {0} guid names for testing".format(nRead))


print("############## checking server connectivity #################")
print("Note: millisecond responses are expected to the first calls")

print("############## checking server time         #################")

startParse=datetime.datetime.now()
retval=client.server_time()
print(retval)
nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
print("Reading server time (whole operation) ",perEvent.total_seconds())

print("############## checking server memory usage   #################")

startParse=datetime.datetime.now()
retval=client.server_memory_usage()
print(retval)

nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading server memory (whole operation) ",perEvent.total_seconds())

print("###############   checking server config   ###################")

startParse=datetime.datetime.now()
retval=client.server_config()
nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading server memory (whole operation) ",perEvent.total_seconds())

print("################### exist_sample ###########################")
nTested=0
startParse=datetime.datetime.now()
for guid in guids:
	retval=client.exist_sample(guid)
	nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading for guid exists (per operation)",perEvent.total_seconds())

print("################### get_all_filtered_guids ###########################")
nTested=0
startParse=datetime.datetime.now()
retVal=client.get_all_filtered_guids(0.85)
nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading all_filtered_guids (total operation) ",perEvent.total_seconds())


print("################### get_all_guids ###########################")
nTested=0
startParse=datetime.datetime.now()
retVal=client.get_all_guids()
nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading all_guids (total operation) ",perEvent.total_seconds())

print("################### get_all_guids_examination_time ###########################")
nTested=0
startParse=datetime.datetime.now()
retVal=client.get_all_guids_examination_time()

nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading get_all_guids_examination_time (total operation) ",perEvent.total_seconds())


print("################### query_get_value_snp_filter ###########################")
nTested=0
startParse=datetime.datetime.now()
for guid in guids:
	print(guid)
	retVal=client.query_get_value_snp_filter(guid, 12)
	nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit()


print("## query_get_value_snp_filter example returned_format=1 vs 2 #######")
nTested=0
startParse=datetime.datetime.now()
for guid in list(guids)[0:50]:
	retVal=client.query_get_value_snp_filter(guid, 12, 0.66, 1)
	print(nTested, "Format 1: ",retVal)
	retVal=client.query_get_value_snp_filter(guid, 12, 0.66, 2)
	print(nTested, "Format 2: ",retVal)	
	nTested+=1

try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit()


print("################### query_get_detail ###########################")
nTested=0
startParse=datetime.datetime.now()
comparator=guids.pop()

for guid in guids:
	retval=client.query_get_detail(comparator, guid)
	nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("query_get_detail (per operation)",perEvent.total_seconds())

exit()


     
print("################### get_all_annotations ###########################")
nTested=0
startParse=datetime.datetime.now()
retVal=client.get_all_annotations()
nTested+=1
try:
	endParse=datetime.datetime.now()
	perEvent=(endParse-startParse)/nTested

except ZeroDivisionError:
  print("No sequences found")
  exit() 
print("Reading get_all_annotations (whole operation) ",perEvent.total_seconds())



