#!/usr/bin/env python3
""" enters sequences produced by make_simulation.py into a findNeighbour3 server,
    and reports on the results.
    
Loads data produced by make_simulation.py into a running findNeighbour3 server.
Example usage:

# first, a server must be running
python findNeighbour3-server.py ../demos/simulation/config/config.json
# then simulations must be generated (e.g. with run_simulation)
python run_simulation.py  ../output/simulation_set_1
"""

import pandas as pd
import os
import glob
import numpy as np
import json
import argparse
import random
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

from fn4client import fn4Client

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=""" Loads data produced by make_simulation.py into a running findNeighbour3 server.
Example usage:

# first, a server must be running
python findNeighbour3-server.py ../demos/simulation/config/config.json
# then simulations must be generated (e.g. with run_simulation)
python run_simulation.py  ../output/simulation_set_1""")
    
    
    parser.add_argument('inputdir', type=str, nargs=1, help='data will be read from the inputdir')    
    args = parser.parse_args()
    basedir = os.path.abspath(args.inputdir[0])
    
    
    # connect to server
    fn4c = fn4Client("http://localhost:5020")
    
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
  
        print("Loading new data")
        fn4c.reset()        # restart server with no data and new reference sequence
        
        inserted = []
        for ix in observed.index:
            to_insert = observed.loc[ix,'seq']
            guid = observed.loc[ix,'node']
            fn4c.insert(guid, to_insert)
            #print("Loaded ",guid)
            inserted.append(guid)
            if not len(inserted) == len(fn4c.guids()):
                raise KeyError("Guids were not inserted correctly, {0} {1}".format(inserted, fn4c.guids()))

        # examine all samples.  are they mixed?
        print("Conducting per-sample analysis")
        per_sample = []
        for sample in  fn4c.guids():

            neighbours = fn4c.guid2neighbours(guid=sample, threshold = 20)
            nneighbours = len(neighbours)
            guids = [sample]
            for item in neighbours:
                guids.append(item[0])
            if nneighbours < 2:
                status = 'Not assessable'
            else:

                msa = fn4c.msa(guids)
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
    
            isCorrect = True
            if status == 'multiple samples mixed':
                isCorrect = False
            if status == 'Not mixed' and (params['observed_sequence_mixed']==sample):
                isCorrect = False
            res = {'sim': observed_filename, 'guid':sample, 'nneighbours':len(neighbours),
                   'mix_ascertainment':status, 'truly_mixed':params['observed_sequence_mixed']==sample,
                   'isCorrect':isCorrect}
            per_sample.append(res)

        outfile = os.path.join(inputdir,'persample.xlsx')
        if len(per_sample)>0:
            per_sample_results = pd.DataFrame.from_records(per_sample)
            per_sample_results.to_excel(outfile)

            # two outcomes
            # proportion of samples correctly classified (all or none)
            # proportion of samples which cannot be assessed

            xtb = pd.crosstab(per_sample_results['isCorrect'], per_sample_results['mix_ascertainment'] )
            xtb['simulation']=inputdir
            xtb['nSamples']=len(per_sample)
            outfile = os.path.join(inputdir,'persample_xtb.xlsx')
            xtb.to_excel(outfile)
        else:
            print("Caution: no samples found ? path to simulated data is incorrect")
        
        # second analysis
        # start with truly mixed sample
        # find neighbours
        # sample neighbours n=2 to n = ?10.  Up to 10 replicates per set of neighbours.
        # are the samples detected as mixed?
        # report proportion detected in the simulations
        
        # find the neighbours of the mixed sample
        for snv_cutoff in [5]:
            neighbours = fn4c.guid2neighbours(guid=params['observed_sequence_mixed'], threshold = snv_cutoff)
            nneighbours = len(neighbours)
            neighbouring_guids = []
            print("The mixed sample {0} has {1} neighbours at {2} snp; resampling".format(params['observed_sequence_mixed'], nneighbours, snv_cutoff))
            
            # the acid test from a clustering perspective
            # is whether the falsely low relationship between two unmixed neighbours of mixed (M) (A,B)
            # is avoided by mixture detection
            for item in neighbours:
                neighbouring_guids.append(item[0])
            if nneighbours >= 2:
                # we sample the neighbours
                max_neighbours = nneighbours
                if max_neighbours > 20:
                    max_neighbours = 20
                sampled = []
                for s in range(2,nneighbours):
                    max_samples = math.factorial(s)     # max number of samples possible
                    if max_samples > 5:
                        max_samples = 5
                    for i in range(max_samples):
                        guid_sample = random.sample(neighbouring_guids, s)
                        guid_sample.append(params['observed_sequence_mixed'])
                        msa = fn4c.msa(guid_sample)
                        #print(msa)
                        p_value_cutoff = 1e-5/len(msa.index)
                        ismixed = msa.query("p_value1 < @p_value_cutoff")
                        status = 0
                        if len(ismixed.index)==1:
                            if params['observed_sequence_mixed'] in ismixed.index:
                                status = 1
                        sampling_result = {'dataset':inputdir, 'snv_cutoff':snv_cutoff, 'sample': i, 'sample_size':s, 'assessed_mixed': status, 'distance_between_mixed':params['distance_between_mixed']}
                        sampled.append(sampling_result)
                        #exit()
            
                outfile = os.path.join(inputdir,'resampling.xlsx')
                if len(sampled)>0:
                    sampled_df = pd.DataFrame.from_records(sampled)
                    sampled_df.to_excel(outfile)
                    print(sampled_df)
                    
