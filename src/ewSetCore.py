#!/usr/bin/env python3
import logging
import sqlalchemy
import uuid
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Float, DateTime, MetaData, select, func, or_
from sqlalchemy import create_engine, exc as excb
from sqlalchemy.orm import sessionmaker
from sqlalchemy import ForeignKey
from sqlalchemy.orm import relationship, backref, exc
import datetime
import unittest
import os
import glob
import sys
from Bio import SeqIO
import datetime
import pickle
import hashlib
import collections
import uuid
import json

from seqComparer import seqComparer
from ewsnpstore import db_ewss, ElephantWalkDataSNPaccessor



class ewSetCore():
    """ implements reference compressed, set based in-memory comparisons of reference mapped sequences
    
    The implementation is in python; tested with v. 3.5 and 2.7.
    persistence is to disc as pickled files.
    
    The names of the methods are the same ones used by the class findneighbour in module webservice-server at 
    https://github.com/davidhwyllie/findNeighbour/blob/master/webservice-server.py.
    
    """
    def __init__(self, CONFIG):
        """ instantises the ewSet_core class, which loads sequences from disc into memory,
        and provides methods for compressing, adding, new sequences as well as
        computing, storing and searching for SNP distances between sequences.
        
        A number of parameters are required to do this.
        These need to be passed as a 'configuration dictionary' providing key-value pairs.
    
        __init__ checks that the required key are present:
            INPUTREF:       the path to fasta format reference file.
            PERSISTENCEDIR: the directory where the reference compressed sequences are stored.  Must be writable by the server daemon.
            EXCLUDEFILE:    a file containing the zero-indexed positions in the supplied sequences which should be ignored in all cases.
                            Typically, this is because the software generating the mapped fasta file has elected not to call these regions,
                            in any samples, e.g. because of difficulty mapping to these regions.
                            Such regions can occupy up 5- 20% of the genome and it is important for efficient working of this software
                            that these regions are supplied for exclusion on sequence loading.  Not doing so will slow loading, and markedly increase
                            memory requirements, but will not alter the results produced.
            SNPDIR:         if an SQLite database is used to persist edges, this is the directory in which it will be placed.  Must be writable by the server daemon.
            DEBUGMODE:      False by default.  If true, will stop loading after 500 sequences have been incorporated.
            SERVERNAME:     the name of the SQLite database, if used.
            EDGEDB_CONNSTRING: a valid SQL alchemy connection string.  Will accept one parameter, SERVERNAME, formatted using .format.
                            If <<DEFAULT>> is included in t he string, SNPDIR will be incorporated.  Example:
                            "sqlite:///<<DEFAULT>>/{0}.db"
                            SQL server connection strings have also been tested and work correctly.
            MAXN:           The maximum number of Ns in the sequence <excluding those defined in > EXCLUDEFILE which should be indexed.
                            Other files, e.g. those with all Ns, will be tagged as 'invalid'.  Although a record of their presence in the database
                            is kept, they are not compared with other sequences.
 
            Example settings are below:
           				CONFIG={
						'PORT':8184,
						'IP':'localhost'
						'INPUTREF':os.path.join("/home/compass/ELEPHANTWALK2/","reference","reference.fasta"),
						'PERSISTENCEDIR':"/home/compass/R039/data/",
						'EXCLUDEFILE':os.path.join("/home/compass/R039/","reference","TB.txt"),
						'SNPDIR':"/home/compass/R039/sqlite",
						'DEBUGMODE':0,
						'SERVERNAME':'TBSNP',      
						'EDGEDB_CONNSTRING':"sqlite:///<<DEFAULT>>/{0}.db",
						'MAXN_STORAGE':100000,
	                }
      
        As part of startup, the __init__ method will load all reference compressed sequences from PERSISTENCEDIR into RAM.
        Typical load rate is 100-200 sequences/ second (i.e. between 6,000 and 12,000 sequences per minute).
        Therefore, server restart can take a few minutes.  Doing so is normal.
        
        Exemplar usage follows:
        if __name__=="__main__":
        
            # configuration parameters relevant to the server
            CONFIG={    
                'INPUTREF':os.path.join("..","reference","reference.fasta"),
                'PERSISTENCEDIR':"/home/local/GEL/dwyllie/data/ew2/compressedseqs_test",
                'EXCLUDEFILE':os.path.join("..","reference","TB.txt"),
                'SNPDIR':"/home/local/GEL/dwyllie/data/ew2/snpstore",
                'DEBUGMODE':True,
                'SERVERNAME':'TBSNP',      
                'EDGEDB_CONNSTRING':"sqlite:///<<DEFAULT>>/{0}.db",
                'MAXN':100000
                }     
        
            # instantiate
            ewsc=ewSetCore(CONFIG=CONFIG)
            
            # we are going to read the sequences from a test file provided by Trien.
            testFile="/mnt/microbio/pipeline/TB-FASTA/maskedSamples.txt" 
            nTested=0
            guidsAdded=list()
            with open(testFile, 'rt') as f:
                for line in f:
                    (tag,guid,seq)=line.split('\t')
                    nTested+=1
                    seq=seq.rstrip(' ')
                    seq=seq.rstrip('\n')
                    if CONFIG['DEBUGMODE']==True and nTested>50:
                        break
                    if nTested % 50==0:
                        print("Loaded {0}".format(nTested))
                    guidsAdded.append(guid)                
                    # insert the sequence into the store.
                    ewsc.insert(guid=guid, seq=seq)
        
            # now exemplify recovery;
            print("Testing for guid presence:")
            for guid in guidsAdded:
                print(guid, ewsc.exist_sample(guid))
            
            print("Returning edges < 5")
            for item in ewsc.query_get_values_snp(cutoff=5):
                print(item)
                
            print("Returning neighbours of a single guid with cutoff < 5.")
            for guid in guidsAdded:
                print("Neighbours of {0}".format(guid))
                for item in ewsc.query_get_values(guid=guid, cutoff=5):
                    print(item)
        
            print("Returning neighbours of a pair of guids")
            for guid1 in guidsAdded:
                for guid2 in guidsAdded:
                    print(guid1, guid2)
                    result=ewsc.query_get_value(guid1=guid1, guid2=guid2)
                    print(result)
                                
            print('Finished')
        
        """
        
		# logging
        logging.getLogger()
		
        # validate input
        required_keys=set(['INPUTREF','PERSISTENCEDIR','EXCLUDEFILE','SNPDIR','DEBUGMODE','SERVERNAME','EDGEDB_CONNSTRING','MAXN_STORAGE'])
        missing=required_keys-set(CONFIG.keys())
        if not missing==set([]):
            raise KeyError("Required keys were not found in CONFIG. Missing are {0}".format(missing))
        
        for key in ['INPUTREF','PERSISTENCEDIR','EXCLUDEFILE']:
            if not os.path.exists(CONFIG[key]):
                raise IOError("File/Directory specified by key {0} does not exist: {1}".format(key, CONFIG[key]))
         
        # set properties relevant to configuration.   
        self.CONFIG=CONFIG
        self.connString=self.CONFIG['EDGEDB_CONNSTRING'].format(self.CONFIG['SERVERNAME'])
        
        # read reference
        with open(self.CONFIG['INPUTREF'],'rt') as f:
            for r in SeqIO.parse(f,'fasta'):
                self.reference=str(r.seq)
                
        # create sequence comparer object using reference, and a supply a directory to store the persistent objects;
        # instantiation will load the profiles, but not the full sequences, into memory.
        # defaults are supplied by the constructor about what proportion of bases are to be sampled.
        self.sc=seqComparer(reference=self.reference,
                           persistenceStore='localmemory',
                           persistenceDir=self.CONFIG['PERSISTENCEDIR'],
                           excludeFile=self.CONFIG['EXCLUDEFILE'],
                           startAfresh=False,
						   NCompressionCutoff=self.CONFIG['NCOMPRESSIONCUTOFF'],
						   snpCeiling=self.CONFIG['SNPCEILING'],
						   maxNs=self.CONFIG['MAXN_STORAGE'],
                           debugMode=self.CONFIG['DEBUGMODE'])
        
        # initiate elephantwalk storage system, for persisting SNPs; supply temporary directory and sqlite filename 
        self.ewdir=os.path.join(self.CONFIG['SNPDIR'])
        self.dbname=self.CONFIG['SERVERNAME']
        
        # delete any existing sqlite db if present and we are in DEBUGMODE
        if self.CONFIG['DEBUGMODE']==True:
                logging.warning("Deleting existing database as debugMode=True;")  
                dbfilename=os.path.join(self.CONFIG['SNPDIR'],"{0}.db".format(self.dbname))
                if os.path.exists(dbfilename):
                    os.unlink(dbfilename)
       
        # create EW edge storage object
        self.ewc=ElephantWalkDataSNPaccessor(db=db_ewss,
											 engineName=self.connString,
											 ewdir=self.ewdir,
											 maxDistance=self.CONFIG['SNPCEILING'],
											 refName=self.CONFIG['SERVERNAME'])

        # conduct integrity check.
        inRAM  =self.sc.guidscachedinram()     # this should be a list of guids which have been examined.
        onDisc=self.sc.guidscachedtofile()
        inDb=set(self.ewc.subsetOfNodes(examinationStatus=1))      # the ones which have been examined.
        
        # is there anything in RAM which is not also on disc?
        if not len(inRAM-onDisc)==0:
            raise ValueError("Integrity check #1 failed")
        if not len(inDb - inRAM)==0:        # find things which are inDB but not in RAM.  
            raise ValueError("Integrity check #2 failed")

    def insert(self, guid, seq):
        """ compresses seq, which is a string comprising the sequence, identified by guid, and store it in RAM and on disc."""

        co=None
        if not self.sc.iscachedinram(guid):      # if the guid is not already there
            if not self.sc.iscachedtofile(guid): # we have not previously compressed it;
                co=self.sc.compress(seq)         # so we compress it
            else:
                co=self.sc.load_fromfile(guid)   # otherwise we load the copy we made before;
            self.sc.persist(refCompressedSequence=co, guid=guid, method='localmemory')     # store the parsed object in ram and ensure it is on disc

        # now compute the genetic neighbours of this new sequence.
        links={'guid':guid, 'neighbours':[]}
        for key2 in self.sc.guidscachedinram():
                if not guid==key2:
                    (guid1,guid2,dist,n1,n2,nboth, N1pos, N2pos, Nbothpos)=self.sc.countDifferences_byKey(keyPair=(guid,key2))
                    #print(guid1,guid2,dist,n1,n2,nboth, N1pos, N2pos, Nbothpos)
                    if dist is not None:
                        links['neighbours'].append([guid2,dist, n1,n2,nboth])

        self.ewc.storeRelations(links)      # store them in db.

    def exist_sample(self, guid):
        """ returns True or False depending on whether the guid has been stored.  """
        return(self.ewc.nodeExists(guid))

    def query_get_values_snp(self, cutoff=20):
        """ returns all pairwise distances between samples <= distance cutoff """
        return(self.ewc.edges(cutoff=cutoff))

    def query_get_values(self, guid, cutoff=20, returned_format=1):
        """ returns the pairwise links associated with one guid.

		returned_format=1 returns the contacts in a legacy format delivered by EW (default),
			[otherGuid, dist]
			
		returned_format=2 returns a more extended version:
		    [otherGuid, dist, result.N_just1, result.N_just2, result.N_either]
		"""
        return self.ewc.neighboursOf(guid=guid, cutoff=cutoff, returned_format=returned_format)

    def query_get_value(self, guid1, guid2):
        """ returns the pairwise distance (if stored) associated with one pair of guids"""
        return(self.ewc.pairwiseComparison(guid1=guid1, guid2=guid2))

    def query_get_detail(self, guid1, guid2):
        """ returns the pairwise distance, number and positions of Ns
        provided the sequences are of sufficient quality to be stored at all)
        associated with one pair of guids"""	
        guid1_exists=self.ewc.nodeExists(guid1)
        guid2_exists=self.ewc.nodeExists(guid2)
        retVal=({'guid1_exists':guid1_exists, 'guid2_exists':guid2_exists, 'success':0})	
        if guid1_exists and guid2_exists:
            (key1, key2, nDiff, n1,n2,nboth, N1pos, N2pos, Nbothpos)=self.sc.countDifferences_byKey((guid1,guid2), cutoff=1e10)
            retVal['success']=1
            retVal['nDiff']=nDiff
            retVal['N_just1']=n1
            retVal['N_just2']=n2
            retVal['N_either']=nboth
            if N1pos is not None:
                retVal['N_just1_positions']=list(N1pos)
            if N2pos is not None:
                retVal['N_just2_positions']=list(N2pos)
            if Nbothpos is not None:
                retVal['N_either_positions']=list(Nbothpos)
        return(retVal)

