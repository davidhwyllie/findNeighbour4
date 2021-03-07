""" illustrates use of findNeighbour4 with covid samples
assumes a findNeighbour4 server is running
"""

if __name__ == '__main__':
    import os
    import glob
    import datetime
    import json
    import pandas as pd
    from random import sample as random_sample
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from collections import Counter
    from fn4client import fn4Client


    def neighbours_of(sample, snp_cutoff=1, max_neighbourhood_mixpore = 50):
        """ obtains a list of neighbouring samples of *sample* , and an indication of whether *sample* is mixed.

        Parameters
        sample      the identifier of the sample to start with
        cutoff      how many snps from sample neighbours should be identified.  Recommend 1 or 2 for COVID-19
        max_neighbourhood_mixpore:  the maximum number of similar samples (defined at being less than or equal to the cutoff) used for mixpore computations (assessing mixtures of similar strains)

        Returns
        dictionary with keys
            mixpore_test_p  : p value indicating whether number of N/M in sites of recent change (differing between the close neighbours of sample) are more common than in the rest of the genome
                              if significant, this is indicative of a mixture of samples, at least in TB
            neighbours:       a list of neighbours of sample
        """

        res = fn4c.guid2neighbours(sample, threshold = snp_cutoff)        # find neighbours within a snv cutoff
        neighbours = []
        for related_sample, distance in res['neighbours']:
                neighbours.append(related_sample)

        # if there are more than 50 samples, randomly downsample.
        if len(neighbours)>50:
            for_msa = random_sample(neighbours, max_neighbourhood_mixpore)
        else:
            for_msa = neighbours
        for_msa.append(sample)          # build the samples into the MSA

        # to just get the MSA
        msa_df = fn4c.msa(for_msa, output_format='json-records', what='N_or_M')
        msa_df.set_index('guid', inplace=True, drop=True) 
        mixpore_test_p = msa_df.at[sample, 'p_value3']
        return {'mixpore_test_p': mixpore_test_p, 'neighbours':neighbours}


    # instantiate client
    fn4c = fn4Client("http://findneighbours04.unix.phe.gov.uk:5023")      # expects operation on local host; pass baseurl if somewhere else.

    snp_cutoff = 5                                 # how many snp from a sample should we include neighbours
    sample_quality_cutoff = 0.9                    # only build tree from samples with >= 90% ACTG

    print("Loading all samples & annotations.")
    sample_annotations = fn4c.annotations()
    all_samples = sample_annotations.index.to_list()
    high_quality_samples = sample_annotations.loc[sample_annotations['DNAQuality:propACTG']>0.9].copy()
    
    sampling_population = high_quality_samples.index.to_list()

    all_samples = set(all_samples)
    sampling_population = set(sampling_population)          # samples now (note: more may be added during computations by other processes, so we pick a set to work on at the start of the computations)

    print("From all {0} samples currently present, studying {1} high quality (ACTG > 0.9) samples ".format(len(all_samples), len(sampling_population)))
    representative_subsample = set()
    remaining_samples = sampling_population.copy()
    representative_size = {}

    ## DEBUG 
    # remaining_samples = set(random_sample(list(remaining_samples), 5) )

    mixed_samples = set()
    singletons = set()
    iteration = 0
    total_delta = 0
    while len(remaining_samples)>0:

        iteration +=1
        test_sample = random_sample(list(remaining_samples),1)[0]        # randomly sample one
        test_sample_neighbour_info = neighbours_of(test_sample, snp_cutoff)
        
        # filter these neighbours.   Only include neighbours which are part of all_samples.
        # because samples get added to the server all the time, this is not guaranteed.
        test_sample_neighbours = []
        for item in test_sample_neighbour_info['neighbours']:
            if item in all_samples:
                test_sample_neighbours.append(item)

        # maybe put in place a quality filter (total Ns)
        delta = 1           # one less at a minimum (test_sample is never reconsidered);

        if len(test_sample_neighbours) == 0:
            singletons.add(test_sample)

        elif test_sample_neighbour_info['mixpore_test_p'] is None:
            # can tell, not enough neighbours
            # include it
            representative_subsample.add(test_sample)

        elif test_sample_neighbour_info['mixpore_test_p'] > 1e-5:         # not obviously mixed
            representative_subsample.add(test_sample)
            representative_size[test_sample] = len(test_sample_neighbours)

            # remove sample and neighbours from further consideration
            n1 = len(remaining_samples)
            for item in test_sample_neighbours:
                remaining_samples.discard(item)         # discard as item may be low quality or otherwise pre-excluded
            n2 = len(remaining_samples)
            delta = 1+ n1- n2
            total_delta = total_delta + delta   # total number of samples

        elif test_sample_neighbour_info['mixpore_test_p'] <= 0.05:
            mixed_samples.add(test_sample)

        else:
            raise ValueError("Unhandled situation")

        # remove the test sample, if it has not already been removed
        remaining_samples.discard(test_sample)

        print("{0} | Iteration: {1} | Remaining samples {2} |  Change now/total {9}/{10} | Sample = {3} |  Mixture statistic {4} | Neighbours {5} | Selected representatives {6} | Mixed samples {7} | Singletons {8} | Approx Subsample 1 in {11:.0f}".format(
          datetime.datetime.now().isoformat(),
          iteration,
          len(remaining_samples),
          test_sample,
          test_sample_neighbour_info['mixpore_test_p'],
          len(test_sample_neighbours),
          len(representative_subsample),
          len(mixed_samples),
          len(singletons),
          delta, 
          total_delta,
          total_delta/float(len(representative_subsample))
          
        ))

        #if iteration > 1:      # DEBUG
        #    break

    # export masked fasta (note: holds all selected sequences in ram)
    fasta_outputfile =  'tree_selected.fasta' 
    seqs = []
    for test_sample in representative_subsample:
        seq = fn4c.sequence(test_sample)
        sr = SeqRecord(     Seq(seq['masked_dna']), 
                                id= test_sample,
                                description="| within {0} snp of {1} others".format(snp_cutoff, representative_size[test_sample]) 
                                )        
        seqs.append(sr)

    with open(fasta_outputfile, 'w') as f:
        SeqIO.write(seqs, f, "fasta")
    print("Selected samples are written to {0}".format(fasta_outputfile))
    print("Finished")
    