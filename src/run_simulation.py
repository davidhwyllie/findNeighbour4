#!/usr/bin/env python3
""" enters sequences produced by make_simulation.py into a findNeighbour3 server,
    and reports on the results.
    
"""

import pandas as pd
import os
import glob
import numpy as np
import json
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

from fn3client import fn3Client

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=""" Loads data produced by make_simulation.py into a running findNeighbour3 server.
Example usage:

# first, a server must be running
python findNeighbour3-server.py ../demos/simulation/config/config.json
# then
python run_simulation.py  ../output/simulation_set_1""")
    parser.add_argument('inputdir', type=str, nargs=1, help='data will be read from the inputdir')    
    args = parser.parse_args()
    basedir = os.path.abspath(args.inputdir[0])
    
    
    # connect to server
    fn3c = fn3Client()
    
    # iterate over simulated data
    for inputdir in glob.glob(os.path.join(basedir,'*')):
        print(inputdir)

        # define filenames
        fasta_filename = os.path.join(inputdir,'phylogeny.fasta')
        sequence_filename = os.path.join(inputdir,'phylogeny.txt')
        observed_filename = os.path.join(inputdir,'observed.txt')
        tree_filename = os.path.join(inputdir,'phylogeny.nwk')
        ref_filename = os.path.join(inputdir,'reference.fasta')
        treepic_filename = os.path.join(inputdir,'{0}.png'.format('tree_image'))
        annotated_treepic_filename = os.path.join(inputdir,'{0}.png'.format('annotated_tree_image')) 
        distmat_filename = os.path.join(inputdir,'distmat.txt')
        obs_distmat_filename = os.path.join(inputdir,'obs_distmat.txt')
        params_filename = os.path.join(inputdir,'params.json')
        
        # read reference sequence
        print("reading h37rv control sequence")
        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                    seq = str(record.seq)
    

  
        # read data
        with open(params_filename,'rt') as f:
            params = json.load(f)
        with open(observed_filename,'rt') as f:
            observed = pd.read_csv(f, sep='\t')

        
        ref_filename = os.path.join("..",'demos','simulation','config','ref.fasta')
        l = params['len_simulated_alignment']
        print("Making reference of length", l)
        miniseq=seq[0:l]
        with open(ref_filename,'wt') as f:
            f.write(">reference\n{0}\n".format(miniseq))
  
        fn3c.reset()        # restart server with no data and new reference sequence
        
        inserted = []
        for ix in observed.index:
            to_insert = observed.loc[ix,'seq']
            guid = observed.loc[ix,'node']
            fn3c.insert(guid, to_insert)
            print("Loaded ",guid)
            inserted.append(guid)
            if not len(inserted) == len(fn3c.guids()):
                raise KeyError("Guids were not inserted correctly")

        # examine each sample.  is it mixed?
        for sample in inserted:
            neighbours = fn3c.guid2neighbours(guid=sample, threshold = 20)
            nneighbours = len(neighbours)
            guids = [sample]
            for item in neighbours:
                guids.append(item[0])
            if nneighbours < 2:
                status = 'Not assessable'
            else:

                msa = fn3c.msa(guids)
                p_value_cutoff = 1e-5/len(msa.index)
                ismixed = msa.query("p_value1 < @p_value_cutoff")
                status = 'Not mixed'
                if len(ismixed.index)==1:
                    if sample in ismixed.index:
                        status = 'Sample mixed'
                if len(ismixed.index)>1:
                    status = 'multiple samples mixed'
                    outfile = os.path.join(inputdir,'msa_{0}.xlsx'.format(sample))
                    msa.to_excel(outfile)
    
                
            print(sample, len(neighbours), status, 'True mix = ',params['observed_sequence_mixed'])
        exit()
            # find its neighbours, possibly with downsampling
        # make msa with sample + neighbours - this is what is proposed 
        # is it mixed?
        msa = fn3c.msa(inserted)
        
        print(ismixed)
        print(params)

    