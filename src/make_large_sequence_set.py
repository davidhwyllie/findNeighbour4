#!/usr/bin/env python3
""" produces a simulated sequences, with ancestors and their descendents, storing output to outputdir
Suitable for generating large sets of data similar to that produced by TB sequencing, for stress testing
relatedness servers.

Algorithm:
==========
for each of x samples, from a reference sequence
-- it introduces k1 variations at random positions, generating an ancestral sequence.
-- for each ancestral sequence
----- it generates 3 children, each with a small number of variations at random positions.
      'small number' is drawn from a Poisson distribution with mean s;
      children are then mutated as above, until up to c children of the ancestor are produced.
for each output sequence, it generates Ns (unknown bases) at random positions in the consensus at frequency p.
write output to outputdir

"""

import ete3
import random
import pandas as pd
import os
import numpy as np
import argparse
import json
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

def introduce_n(seq, p):
    """ introduce Ns at random into a sequence, seq.
    The probability that each base is N is p, where 0<=p<=1
    """

    # compute the number of Ns given the length of seq and p
    seqlen = len(seq)
    n = np.random.binomial(seqlen, p)
    
    # compute bases to turn to N
    to_change = random.sample(range(seqlen),n)
    seqlist = list(seq)
    for pos in to_change:
        seqlist[pos] = 'N'
    
    return ''.join(seqlist)
         
def mutate(seq, n):
    """ introduce n mutations into seq """

    # compute bases to turn to something different
    seqlen = len(seq)
    to_change = random.sample(range(seqlen),n)
    seqlist = list(seq)
    for pos in to_change:
        if seqlist[pos]=='G':
            seqlist[pos] = 'T'
        else:
            seqlist[pos] = 'G'
    
    return ''.join(seqlist)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Produces a simulated sequences, with ancestors and their descendents,
storing output to disc.
Suitable for generating large sets of data similar to that produced by TB sequencing, for stress testing
relatedness servers.

Algorithm:
==========
for each of n samples, from a reference sequence
-- it introduces k1 variations at random positions, generating an ancestral sequence.
-- for each ancestral sequence
----- it generates r0 children, each with a small number of variations at random positions.
      'small number' is drawn from a Poisson distribution with mean s;
      children are then mutated as above, until up to max_c children of the ancestor are produced.
for each output sequence, it generates Ns (unknown bases) at random positions in the consensus at frequency p.
write output to outputdir

Example usage:

python make_large_sequence_set.py 10 1000 3 50 1e-8 ../output/lss_tb

""")
    parser.add_argument('n_sequences', type=int, nargs=1,
                        help='number of ancestors sequences to simulate. Suitable value = 10')
    parser.add_argument('k1', type=int, nargs=1, 
                        help='number of nucleotides by which ancestor sequences differ from reference. Suitable value = 1000')
    parser.add_argument('s', type=int, nargs=1,
                        help='the mean number of SNV by which children differ from their parent.  Suitable value = 2')
    parser.add_argument('max_c', type=int, nargs=1,
                        help='The maximum number of sequences produced from each ancestor.  The total number of sequences generates is n_sequences*max_c.  Suitable value = 50')
    parser.add_argument('p', type=float, nargs=1,
                        help='probability that a single base is N.  Suitable value= 0 (if no Ns are to be introduced)')    
    parser.add_argument('outputdir', type=str, nargs=1,
                        help='output will be written to the outputdir')    
    args = parser.parse_args()
    n_sequences = args.n_sequences[0]
    outputdir = os.path.abspath(args.outputdir[0])
    k1 = args.k1[0]
    s= args.s[0]
    max_c = args.max_c[0]
    p = args.p[0]
    # if the output directory does not exist, create it
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    print("reading h37rv control sequence")
    inputfile = "../COMPASS_reference/R39/R00000039.fasta"
    with open(inputfile, 'rt') as f:
        for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                refseq = str(record.seq)
             
    for i in range(n_sequences):
        seqs = {}
        print("Building tree")  # we only use the topology
        t= ete3.Tree()
        t.populate(max_c)
        node_id = 0
        for node in t.traverse():
            node_id+=1
            node.name = "N{0}_{1}".format(i, node_id)
            if node.is_root():
                seqs[node.name]= mutate(refseq, k1)
                print("Mutating ancestor # {0}".format(i))
            else:
                dist =  np.random.poisson(s)
    
                seqs[node.name] = mutate(seqs[node.up.name], dist)
               
                if node_id >= max_c:
                    break
                
                if node_id % 10 == 0:
                    print("generated {0} children ..".format(node_id))
        
        print("Exporting {0} sequences, with Ns added".format(node_id))
        for seq_name in seqs.keys():
            with open(os.path.join(outputdir, "{0}.fasta".format(seq_name)),'wt') as f:
                seq = introduce_n(seqs[seq_name],p)
                f.write(">{0}\n{1}\n".format(seq_name,seq))
    
    print("Export complete.")