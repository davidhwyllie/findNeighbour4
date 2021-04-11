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


A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.


"""

import ete3
import random
import os
import numpy as np
import argparse
import itertools
import glob

from string import ascii_uppercase
from Bio import SeqIO

class LargeSequenceSetGenerator():
    """ generates sets of closely related samples.
    
    produces a simulated sequences, with ancestors and their descendents, storing output to outputdir
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
    """
    def __init__(self, refSeq):
        """ makes sequences derived from refSeq
        
        Inputs:
            refSeq  the reference sequence to start with.  should be in string format.
        
        Returns:
            None
        """
        self.refSeq = refSeq
        
    def make_sequences(self, prefix, n_sequences, k1, s, max_c, p):
        """ makes a series of closely related samples.
        
        Inputs:
            prefix  a prefix to apply to the name of the sequence returned.  Allows the calling function to enforce unique naming.
            n_sequences number of ancestors sequences to simulate. Suitable value = 10
            k1  number of nucleotides by which ancestor sequences differ from reference. Suitable value = 1000
            s   the mean number of SNV by which children differ from their parent.  Suitable value = 2
            max_c   The maximum number of sequences produced from each ancestor.  The total number of sequences generates is n_sequences*max_c.  Suitable value = 50')
            p       probability that a single base is N.  Suitable value= 0 (if no Ns are to be introduced)    
  
        Returns:
            a generator providing (sequence_id, sequence) pairs

        """
        for i in range(n_sequences):
            seqs = {}
            
            # build random tree; we only use the topology
            t= ete3.Tree()
            t.populate(max_c)
            node_id = 0
            for node in t.traverse():
                node_id+=1
                node.name = "{0}_{1}".format(i, node_id)
                if node.is_root():
                    seqs[node.name]= self.mutate(self.refSeq, k1)
                    
                else:
                    dist =  np.random.poisson(s)
        
                    seqs[node.name] = self.mutate(seqs[node.up.name], dist)
                   
                    if node_id >= max_c:
                        break         
           
            for seq_name in seqs.keys():
                yield("{1}_{0}".format(seq_name, prefix),  self.introduce_n(seqs[seq_name],p))

        
    def possible_prefixes(self):
        for size in itertools.count(1):
            for s in itertools.product(ascii_uppercase, repeat=size):
                yield "".join(s)
    
    def introduce_n(self,seq, p):
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
             
    def mutate(self, seq, n):
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
    # command line interface; writes data to disc.
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
    nWrit = 0
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    print("reading h37rv control sequence")
    inputfile = "../COMPASS_reference/R39/R00000039.fasta"
    with open(inputfile, 'rt') as f:
        for record in SeqIO.parse(f,'fasta'):               
                refSeq = str(record.seq)
 
    # find unique prefix
    lsg = LargeSequenceSetGenerator(refSeq)
    for prefix in lsg.possible_prefixes():
        globpath = os.path.join(outputdir, "{0}_*.fasta".format(prefix))
        outputfiles = glob.glob(globpath)
        if len(outputfiles)==0:     # the prefix is unique
                break
    print("Using file prefix", prefix)

    for guid, seq in lsg.make_sequences(prefix, n_sequences, k1, s, max_c, p):
        outputfile = os.path.join(outputdir, "{0}_{1}.fasta".format(prefix,guid))
        with open(outputfile,'wt') as f:
            f.write(">{2}_{0}\n{1}\n".format(guid, seq, prefix))
            nWrit +=1
        if nWrit % 50 == 0:
            print("Sequences written", nWrit)

    print("Export complete.")