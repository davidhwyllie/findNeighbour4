#!/usr/bin/env python
""" performs a large scale test of matrix storage.

findNeighbour3 stores a relatedness matrix in a mongodb data store.
This code generates a large number of variants.
Storing these in RAM is very cheap, but the matrix stored will include all pairwise differences.

To run the test, start up a server
python findNeighbour3-server.py ../config/large_matrix_config.json

And then run this script
demo_large_matrix.py

This script determines 
"""
import copy
import os
import random
import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from fn3client import fn3Client

if __name__ == '__main__':
    
    # maximum number to add
    max_sequences = 2000
    
    # H37Rv reference sequence
    print("reading h37rv control sequence")
    inputfile = "../COMPASS_reference/R39/R00000039.fasta"
    with open(inputfile, 'rt') as f:
        for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                seq = str(record.seq)
    seqbases = list(seq)        # a list with one nt per element.

    print("opening connection to fn3 server")
    fn3c = fn3Client()

    # determine all masked positions
    excluded_positions = fn3c.nucleotides_excluded()
    # we can mutate any positions which are not masked
    available_positions = sorted(list(set(range(len(seq)))-set(excluded_positions['excluded_nt'])))
    print("There are {0} available positions to mutate".format(len(available_positions)))
   
    # determine how many samples there are currently in the server.
    nSamples = len(fn3c.guids())
    print("There are {0} existing guids.  Adding more ..".format(nSamples))

    # create output file with header line
    outputfile = os.path.join('..','demos','large_matrix_1', 'output', 'timings_{0}.tsv'.format(nSamples))
    
    with open(outputfile, 'w+t') as f:
        output_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format('nSamples', 's_insert', 'e_insert', 'd_insert', 'r_read', 'e_read', 'd_read')
        f.write(output_line)        
        while nSamples < max_sequences:   
            nSamples +=1
            
            if nSamples>len(available_positions):
                print("All available positions have been mutated")
                exit(0)
    
          
            # mutation nSamples bases, starting at the nSample th position in available_sequences
            new_seqbases = copy.copy(seqbases)
            
            for i in [0]:                # just make one mutation per sample.  range(nSamples):
                current_base  = set(seqbases[i+available_positions[nSamples]])
                non_current_base = set(['A','C','T','G']) - current_base
                can_mutate_to = list(non_current_base)
                random.shuffle(can_mutate_to)
                mutated_base = can_mutate_to[0]
                new_seqbases[i+available_positions[nSamples]]= mutated_base
                
            mutseq = ''.join(new_seqbases)
            guid = 'guid_{0}'.format(nSamples)
            
            # add
            print("Inserting", guid)
            stime1 = datetime.datetime.now()
            fn3c.insert(guid=guid, seq=mutseq)
            etime1 = datetime.datetime.now()
            delta1= etime1-stime1
            
            # recover neighbours of guid1
            print("Recovering", guid)
            stime2 = datetime.datetime.now()
            fn3c.guid2neighbours(guid, threshold=12)
            etime2 = datetime.datetime.now()
            delta2 = etime2 - stime2
            
            output_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(nSamples, stime1, etime1, delta1, stime2, etime2, delta2)
            f.write(output_line)
            f.flush()
            
        print("Have added {0} sequences, stopping.".format(nSamples))
        exit(0)              