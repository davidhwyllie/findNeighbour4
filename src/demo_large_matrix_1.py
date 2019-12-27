#!/usr/bin/env python
""" performs a large scale test of matrix storage.

findNeighbour3 stores a relatedness matrix in a mongodb data store.
This code generates a large number of variants.
Storing these in RAM is very cheap, but the matrix stored will include all pairwise differences.

To run the test, start up a server, e.g.
1 python3 findNeighbour3-server.py ../config/large_matrix_config.json 
Then run the test; this configuration does not pause during insertion.
2 python3 demo_large_matrix_1.py 2000 ../demos/large_matrix_1/output 0 0 
Then analyse the output
3 Rscript demo_depict_timings.R ../demos/large_matrix_1/output


"""
import copy
import os
import random
import datetime
import argparse
import pathlib
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from fn4client import fn4Client

if __name__ == '__main__':
    
    # maximum number to add
    parser = argparse.ArgumentParser(description='Generate a large number of similar sequences')
    parser.add_argument('max_sequences', type=int, nargs=1,
                        help='sequences will be added until max_sequences exist in the server.')
    parser.add_argument('outputdir', type=str, nargs=1,
                        help='output will be written to the outputdir')
    parser.add_argument('pause_after', type=int, nargs=1,
                        help='insertion will pause after adding pause_after sequences (set to zero to never pause)')
    parser.add_argument('pause_duration', type=int, nargs=1,
                        help='pause_duration in seconds')

    args = parser.parse_args()
    max_sequences = args.max_sequences[0]
    pause_after = args.pause_after[0]
    pause_duration = args.pause_duration[0]
    outputdir = os.path.abspath(args.outputdir[0])

    # make the directory
    p = pathlib.Path(outputdir)
    p.mkdir(parents=True, exist_ok=True)
    
    # H37Rv reference sequence
    print("reading h37rv control sequence")
    inputfile = "../COMPASS_reference/R39/R00000039.fasta"
    with open(inputfile, 'rt') as f:
        for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                seq = str(record.seq)
    seqbases = list(seq)        # a list with one nt per element.

    print("opening connection to fn3 server")
    fn4c = fn4Client(baseurl = "http://127.0.0.1:5020")

    # determine all masked positions
    excluded_positions = fn4c.nucleotides_excluded()
    # we can mutate any positions which are not masked
    available_positions = sorted(list(set(range(len(seq)))-set(excluded_positions['excluded_nt'])))
    print("There are {0} available positions to mutate".format(len(available_positions)))
   
    # determine how many samples there are currently in the server.
    nSamples = len(fn4c.guids())
    print("There are {0} existing guids.  Adding more ..".format(nSamples))

    # create output file with header line
    outputfile = os.path.join(outputdir, 'timings_{0}.tsv'.format(nSamples))
    nAdded_this_batch = 0
    with open(outputfile, 'w+t') as f:
        output_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format('nSamples', 's_insert', 'e_insert', 'd_insert', 's_read', 'e_read', 'd_read')
        f.write(output_line) 
        while nSamples < max_sequences:   
            nSamples +=1  
            nAdded_this_batch +=1       
        
            if pause_after >0:
                if nAdded_this_batch % pause_after==0:
                     print("Insert paused; will resume after {0} seconds.".format(pause_duration))
                     time.sleep(pause_duration)
        
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
            print("Inserting", guid, "(samples in this batch = ",nAdded_this_batch,"); will pause every ",pause_after)
            stime1 = datetime.datetime.now()
            resp = fn4c.insert(guid=guid, seq=mutseq)
            print("Insert yielded status code {0}".format(resp.status_code))

            etime1 = datetime.datetime.now()
            delta1= etime1-stime1
            
            # recover neighbours of guid
            stime2 = datetime.datetime.now()
            # check it exists
            if not fn4c.guid_exists(guid):
                print("Guid {0} was not inserted".format(guid))   
            else:
                neighbours = fn4c.guid2neighbours(guid, threshold=10000000)
                etime2 = datetime.datetime.now()
                delta2 = etime2 - stime2
                print("Recovered {1} neighbours of {0}".format(guid, len(neighbours)))            
                output_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(nSamples, stime1, etime1, delta1, stime2, etime2, delta2)
                f.write(output_line)
                f.flush()
            
        print("Have added {0} sequences, stopping.".format(nSamples))
        exit(0)              
