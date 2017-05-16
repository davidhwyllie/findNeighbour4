# ew1 vs ew2
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
from gapi import geneticsAPIInterface
from fastaMasker import fastaMasker

# create a fastaMasker object
excluded_genes = ('rrl','rrs','rpoB', 'Rv2082')
fr=fastaMasker(proscribedGenesFile='../reference/proscribedGenes.txt')
fr.proscribeGenes(excluded_genes)

fastapath = "/mnt/md0/dwyllie/data/relatednesstest/TB_FASTA/{0}_v3.fasta"

json_config_file=sys.argv[1]
with open(json_config_file, 'rt') as f:
	txt=f.read()
	print("Read config file: {0}".format(txt))
	CONFIG=json.loads(txt)

if not type(CONFIG)==dict:
	raise KeyError("CONFIG must be either a dictionary or a JSON string encoding a dictionary.  It is: {0}".format(CONFIG))
	
client=xmlrpc.client.ServerProxy("http://{0}:{1}".format(CONFIG['IP'],CONFIG['PORT']))

print("##################### opening ew1 connection ##############################")
print("Running self test")
gapi = geneticsAPIInterface()
gapi.selfTest_ew()

print("Passed, connection is functional")
print("################### query_get_all_guids from ew2 ###########################")
nTested=1
startParse=datetime.datetime.now()

#retval=client.get_all_values()
ew2_guids=set()
retval=client.get_all_guids()
retlist=json.loads(retval)
for d in retlist:
	ew2_guids.add(d['guid'])

print("################### query_get_all_guids from ew1 ###########################")
retlist = gapi.ewAllGuids()
ew1_guids=set()
for d in retlist:
	ew1_guids.add(d['guid'])

print("################### summary ###########################")
print("There are {0} guids from ew1".format(len(ew1_guids)))
print("There are {0} guids from ew2".format(len(ew2_guids)))
both_guids = ew1_guids.intersection(ew2_guids)
print("There are {0} guids in both".format(len(both_guids)))

# get the guids found for both data sources.
# the guids from ew1 should be identical, or a superset, or those from ew2.

def process_rc(rc):
	# reproduce the gapi procedure for ew2
	#print(rc, type(rc))
	if rc == ['Err','Error in the server side']:
	   raise requests.exceptions.HTTPError("Elephant walk server not working: returned 'Error on the server side' in response to {0}".format(absurl))

	if rc == ['Err', 'missing sample']:
	   # we treat this as an error
	   raise KeyError("Elephant walk server was asked for neighbours of {0} but this guid is not present in the server associated with {1}".format(guid,reference))

	if rc == ['Bad', 'bad sequence']:
	   return(['Bad',[]])
	retVal=json.loads(rc)
	#print(retVal,type(retVal))
	return retVal

def compute_exact_difference(guid1, guid2):
	""" loads and compares sequences of guid1 with guid2, with masking """
	ff1 = fastapath.format(guid1)
	ff2 = fastapath.format(guid2)
	with open(ff1,'rt') as f:
		for r in SeqIO.parse(f,'fasta'):
			ff1seq=str(r.seq)

	with open(ff2,'rt') as f:
		for r in SeqIO.parse(f,'fasta'):
			ff2seq=str(r.seq)
	
	# mask them
	ff1seq = fr.maskString(ff1seq)
	ff2seq = fr.maskString(ff2seq)
	nDiffExact = 0
	for i in range(len(ff1seq)):
		if not ff2seq[i] in ('N','-'):
			if not ff1seq[i] in ['N','-']:
				if not ff1seq[i] == ff2seq[i]:
					nDiffExact +=1
	return(nDiffExact)

print("################### run comparison ###########################")

nTested=0
startParse=datetime.datetime.now()
nDiscrepant = 0

for this_guid in both_guids:
	ew1_result=gapi.ewNeighboursOf(guid=this_guid, SNPCutoff=20)


	# run comparison
	# build a dictionary for each
	ew1_dict={}
	for (guid,dist) in ew1_result:
		if guid in both_guids:		# we can compare the two
			ew1_dict[guid]=dist

	ew2_result=process_rc(client.query_get_value_snp_filter(this_guid, 20))
	ew2_dict={}
	for (guid,dist) in ew2_result[1]:
		if guid in both_guids:
			ew2_dict[guid]=dist
		else:
			raise KeyError("{0} is not present, although it should be".format(guid))

	if not ew1_dict == ew2_dict:
		nDiscrepant+=1

		
		if ew1_dict.keys() == ew2_dict.keys():
			
			for guid in ew1_dict.keys():

				if not ew1_dict[guid]==ew2_dict[guid]:
					nDiffExact = compute_exact_difference(this_guid, guid)
					print("{0}\tDiscordant guid pair ({1},{2})\tEW1={3}\tEW2={4}\tdelta={5}\texact={6}\tpass = {7}".format( nTested,
																						 this_guid,
																						 guid,
																						 ew1_dict[guid],
																						 ew2_dict[guid],
																						 ew1_dict[guid]-ew2_dict[guid],
																						 nDiffExact,
																						 nDiffExact == ew2_dict[guid]
																						))
		else:
			
			s1= set(ew1_dict.keys())
			s2= set(ew2_dict.keys())
			#print("Link in EW1 but not EW2", s1-s2, "Link in EW2 but not EW1", s2-s1)
			for guid in s2-s1:
					  
				nDiffExact = compute_exact_difference(this_guid, guid)																			
				# compute a detailed pairwise comparison
				retval=client.query_get_detail(this_guid, guid)
				resDict = json.loads(retval)
				
				print("{0}\tVariation between guid pair ({1},{2})\tEW1= >cutoff\tEW2= {3}\tEW2 detail= {4} \texact = {5}\tpass = {6}".format(
					nTested,
					this_guid,
					guid,
					ew2_dict[guid],
					resDict['nDiff'],
					nDiffExact,
					ew2_dict[guid] == nDiffExact )
					)

			for guid in s1-s2:

				nDiffExact = compute_exact_difference(this_guid, guid)																			
				# compute a detailed pairwise comparison
				retval=client.query_get_detail(this_guid, guid)
				resDict = json.loads(retval)
				
				print("{0}\tVariation between guid pair ({1},{2})\tEW1={3}\tEW2= >cutoff\tEW2 detail = {4}\texact = {5}\tpass = {6}".format(
					nTested,
					this_guid,
					guid,
					ew1_dict[guid],
					resDict['nDiff'],
					nDiffExact,
					"unknown")
					)
			
	nTested+=1
	if nTested>1000:
		break
	
print("Found {0} / {1} guids discrepant.".format(nDiscrepant, nTested))
exit()
	