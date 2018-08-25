#!/usr/bin/env python3
""" produces simulated sequences, storing output to disc.

simulates phylogenetic trees,
sequences reflecting the members of the phylogeny,
sequences with Ns added randomly derived from the phylogeny,
and mixtures of sequences from within the phylogeny.
Also simulates the presence of a risk factor in mixed and unmixed sequences.

"""

import ete3
import random
import pyvolve
import pandas as pd
import os
import numpy as np
import argparse
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide

def delete_single_child_internal(t):
    """Utility function that removes internal nodes
    with a single child from tree"""

    for node in t.traverse("postorder"):
        if(not node.is_leaf() and len(node.get_children()) < 2):
            node.delete()

    if len(t.get_children()) == 1:
        t = t.children[0]
        t.up = None
def birth_only_tree(birth, nsize=10, max_time=None):
    """Generates a uniform-rate birth-only tree.
    Arguments:
        - ``birth`` : birth rate
        - ``nsize``  :  desired number of leaves
        - ``max_time`` : maximum allowed time for evolution
        
    code is discussed
    https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/
    and taken from
    https://gist.githubusercontent.com/maclandrol/09f2b2f47713da9d465a1a7f72078fe7/raw/9874305b9996515fea5d94a90a976c705440fd2e/tree_simul.py

    """

    done = False
    total_time = 0

    # create initial root node and set length to 0
    tree = ete3.ete3.Tree()
    tree.dist = 0.0
    # check if a stopping condition is provided
    if not (nsize or max_time):
        raise ValueError('A stopping criterion is required')

    while True:

        # get the current list of extant species
        leaf_nodes = tree.get_leaves()

        # get required waiting time before next speciation
        wtime = random.expovariate(len(leaf_nodes) / birth)

        # check if stopping criterion is reached then
        # stop loop if yes
        if len(leaf_nodes) >= nsize or (max_time and total_time + wtime >= max_time):
            done = True

        # update the length of all current leaves
        # while not exceeding the maximum evolution time
        max_limited_time = min(
            wtime, (max_time or total_time + wtime) - total_time)
        for leaf in leaf_nodes:
            leaf.dist += max_limited_time

        if done:
            break

        # update total time of evolution
        total_time += max_limited_time

        # add a new node to a randomly chosen leaf.
        node = random.choice(leaf_nodes)
        c1 = ete3.Tree()
        c2 = ete3.Tree()
        node.add_child(c1)
        node.add_child(c2)
        c1.dist = 0.0
        c2.dist = 0.0

    # Label leaves here and also update
    # the last branch lengths

    leaf_nodes = tree.get_leaves()
    leaf_compteur = 1
    for (ind, node) in enumerate(leaf_nodes):
        node.name = 'T%d' % leaf_compteur
        leaf_compteur += 1

    return tree
def birth_death_tree(birth, death, nsize=10, max_time=None, remlosses=True, r=True):
    """Generates a birth-death tree.
    Arguments:
        - ``birth`` : birth rate
        - ``death`` : death rate
        - ``nsize`` : desired number of leaves
        - ``max_time`` : maximum time of evolution
        - ``remlosses`` : whether lost leaves (extinct taxa) should be pruned from tree
        - ``r`` : repeat until success
        
    code is discussed
    https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/
    and taken from
    https://gist.githubusercontent.com/maclandrol/09f2b2f47713da9d465a1a7f72078fe7/raw/9874305b9996515fea5d94a90a976c705440fd2e/tree_simul.py

    """
    # initialize tree with root node
    tree = ete3.Tree()
    tree.add_features(extinct=False)
    tree.dist = 0.0
    done = False
 
    # get current list of leaves
    leaf_nodes = tree.get_leaves()
    curr_num_leaves = len(leaf_nodes)

    total_time = 0
    died = set([])

    # total event rate to compute waiting time
    event_rate = float(birth + death)

    while True:
        # waiting time based on event_rate
        wtime = random.expovariate(event_rate)
        total_time += wtime
        for leaf in leaf_nodes:
            # extinct leaves cannot update their branches length
            if not leaf.extinct:
                leaf.dist += wtime

        if curr_num_leaves >= nsize:
            done = True

        if done:
            break

        # if event occurs within time constraints
        if max_time is None or total_time <= max_time:

            # select node at random, then find chance it died or give birth
            # (speciation)
            node = random.choice(leaf_nodes)
            eprob = random.random()
            leaf_nodes.remove(node)
            curr_num_leaves -= 1

            # birth event (speciation)
            if eprob < birth / event_rate:
                child1 = ete3.Tree()
                child1.dist = 0
                child1.add_features(extinct=False)
                child2 = ete3.Tree()
                child2.dist = 0
                child2.add_features(extinct=False)
                node.add_child(child1)
                node.add_child(child2)
                leaf_nodes.append(child1)
                leaf_nodes.append(child2)
                # update add two new leave
                # (remember that parent was removed)
                curr_num_leaves += 2

            else:
                # death of the chosen node
                if curr_num_leaves > 0:
                    node.extinct = True
                    died.add(node)
                else:
                    if not r:
                        raise ValueError(
                            "All lineage went extinct, please retry")
                    # Restart the simulation because the tree has gone
                    # extinct
                    tree = ete3.Tree()
                    leaf_nodes = tree.get_leaves()
                    tree.add_features(extinct=False)
                    tree.dist = 0.0
                    curr_num_leaves = 1
                    died = set([])
                    total_time = 0

            # this should always hold true
            assert curr_num_leaves == len(leaf_nodes)

    if remlosses:
        # prune lost leaves from tree
        leaves = set(tree.get_leaves()) - died
        tree.prune(leaves)
        # remove all non binary nodes
        delete_single_child_internal(tree)

    leaf_nodes = tree.get_leaves()
    leaf_compteur = 1
    for ind, node in enumerate(leaf_nodes):
        # label only extant leaves
        if not node.extinct:
            # node.dist += wtime
            node.name = "T%d" % leaf_compteur
            leaf_compteur += 1
    return tree
def mix_sequences(seq1, seq2):
    """ simulates a mixture of sequences. 
    converts any nucleotides differing between seq1 and seq2 to N  
    seq1 and seq2 are strings.
    returns a string representing the output of consensus base calling of a mix of seq1 and seq2."""
    
    if not isinstance(seq1,str) and isinstance(seq2, str):
        raise TypeError("Seq1 and seq2 must both be strings")
    if not len(seq1)==len(seq2):
        raise ValueError("Strings seq1 and seq2 must be the same length")
    
    retVal = []
    for i in range(len(seq1)):
        if seq1[i] == seq2[i]:
            retVal.append(seq1[i])
        else:
            retVal.append('N')
    return ''.join(retVal)
def introduce_n(seq, p):
    """ introduce Ns at random into a sequence, seq.
    The probability that each base is N is p, where 0<=p<=1
    """

    # compute the number of Ns given the length of seq and p
    seqlen = len(seq)
    n = np.random.binomial(seqlen, p)
    
    # compute bases to turn to N
    to_change = random.sample(range(n),n)
    seqlist = list(seq)
    for pos in to_change:
        seqlist[pos] = 'N'
    
    return ''.join(seqlist)
def snv_distance(seq1, seq2):
    """ returns the number of SNV different between two sequences """
    ls1 = list(seq1)
    ls2 = list(seq2)
    if not len(ls1)==len(ls2):
        raise ValueError("The two lists must be the same length")
    dist = 0
    for i in range(len(ls1)):
        if ls1[i] == 'N' or ls2[i]== 'N':
            pass
        else:
            if not ls1[i] == ls2[i]:
                dist+=1
    return(dist)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Produce simulated data for testing software.  Produces simulated sequences, storing output to disc.
simulates phylogenetic trees,
sequences reflecting the members of the phylogeny,
sequences with Ns added randomly derived from the phylogeny,
and mixtures of sequences from within the phylogeny.
Also simulates the presence of a risk factor in mixed and unmixed sequences.
Example usage: python3 make_simulation.py 10 15 0.01 0.05 0.1 ../output/t1""")
    parser.add_argument('n_sequences', type=int, nargs=1,
                        help='number of sequences to simulate. Suitable value = 10')
    parser.add_argument('len_simulated_alignment', type=int, nargs=1, 
                        help='length of the simulated sequences to generate. Suitable value = 15')
    parser.add_argument('prop_random_N', type=float, nargs=1,
                        help='proportion of nucleotides which are called N by chance.  Suitable value = 0.01')
    parser.add_argument('prop_risk_not_mixed', type=float, nargs=1,
                        help='proportion of non-mixed samples with a risk factor. Suitable value = 0.05')
    parser.add_argument('prop_risk_mixed', type=float, nargs=1,
                        help='proportion of mixed samples with a risk factor. Suitable value = 0.10')    
    parser.add_argument('outputdir', type=str, nargs=1,
                        help='output will be written to the outputdir')

    args = parser.parse_args()
    n_sequences = args.n_sequences[0]
    outputdir = os.path.abspath(args.outputdir[0])
    l = args.len_simulated_alignment[0]
    prop_random_N = args.prop_random_N[0]
    prop_risk_mixed = args.prop_risk_mixed[0]
    prop_risk_not_mixed = args.prop_risk_not_mixed[0]
        
    # names of output files
    fasta_filename = os.path.join(outputdir,'phylogeny.fasta')
    sequence_filename = os.path.join(outputdir,'phylogeny.txt')
    observed_filename = os.path.join(outputdir,'observed.txt')
    tree_filename = os.path.join(outputdir,'phylogeny.nwk')
    ref_filename = os.path.join(outputdir,'reference.fasta')
    treepic_filename = os.path.join(outputdir,'{0}.png'.format('tree_image'))
    distmat_filename = os.path.join(outputdir,'distmat.txt')
    params_filename = os.path.join(outputdir,'params.json')
      
    # if the output directory does not exist, create it
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
        
    print("reading h37rv control sequence")
    inputfile = "../COMPASS_reference/R39/R00000039.fasta"
    with open(inputfile, 'rt') as f:
        for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                seq = str(record.seq)

    print("Making mini genomes of length", l)
    miniseq=seq[1:l]
    with open(ref_filename,'wt') as f:
        f.write(">reference\n{0}\n".format(miniseq))
        
    # generate a birth death tree
    print("generating birth-death tree")
    bdtree = birth_death_tree(10,0, nsize=n_sequences, remlosses=True)

    # name internal nodes and store distance to root for all nodes
    counter = 0
    t2n = []
    node_names = []
    for node in bdtree.traverse(strategy='postorder'):
        counter +=1
        new_node_name= "0000"+str(counter)
        new_node_name= new_node_name[len(new_node_name)-4:]
        new_node_name = "N"+new_node_name
        node.name = new_node_name
        node_names.append(node.name)
    for node_name in node_names:           
        dtr = 0
        this_node = node&node_name
        while not this_node.is_root():
            dtr = dtr + this_node.dist
            this_node = this_node.up
        t2n.append({'name':node_name,'dist_to_root':dtr})
    nprops = pd.DataFrame.from_records(t2n, index='name')
    nprops.sort_values('dist_to_root', inplace=True)
    nprops['order'] = range(len(nprops.index))
    ts = ete3.TreeStyle()

    # simulate evolution from the H37Rv reference
    bdtree.write(outfile = tree_filename, format=3)
    t = pyvolve.read_tree(tree=bdtree.write(format=3))  # scale_tree = 0.005
    m = pyvolve.Model("nucleotide" )
    p = pyvolve.Partition(models=m, root_sequence = miniseq)
    e = pyvolve.Evolver(partitions = p, tree=t)
    
    print("Running evolution")
    e()
    simulated_sequences = e.get_sequences()
    with open(sequence_filename,'wt') as f:
        for key in sorted(simulated_sequences.keys()):
            f.write("{0}\t{1}\n".format(key,simulated_sequences[key]))
    with open(fasta_filename,'wt') as f:
        for key in sorted(simulated_sequences.keys()):
            f.write(">{0}\n{1}\n".format(key,simulated_sequences[key]))

    print("Computing pairwise distance matrix")
    dist_list = []
    for i in simulated_sequences.keys():
        for j in simulated_sequences.keys():
            if i<j:
                dist_list.append({'node1':i,'node2':j,'snv':snv_distance(simulated_sequences[i], simulated_sequences[j])})

    distmat = pd.DataFrame.from_records(dist_list)
    distmat.to_csv(distmat_filename)
    
    # what pairwise distances are present?
    dists = sorted(distmat['snv'].unique())
    #print("Pairwise SNV are as follows:", dists)

    #  depict tree
    for i in range(len(nprops.index)):
        for node in bdtree.traverse():
            
            nstyle = ete3.NodeStyle()
            node_order = nprops.loc[node.name, 'order']
            nstyle['size']=0
            nstyle['fgcolor']='black'
            nstyle['bgcolor']='white'
            #node.add_face(ete3.TextFace('-'), column=1, position = 'branch-right')
            if node.is_leaf():
                if node_order < i:
                    #node.add_face(ete3.TextFace(str(i)), column=1, position='branch-right')
                    nstyle['size']=1
                    nstyle['fgcolor']='blue'
                    nstyle['bgcolor']='white'                        

            node.set_style(nstyle)
    bdtree.render(treepic_filename, tree_style= ts,units='px', w=928, h=928) 

        
    print("Generating mixture")
    # select two nodes at random from the available distances
    # this is stratified sampling, ensuring that all parts of the snv distribution are sampled
    distneeded = random.sample(dists, 1)
    pairwise_available = distmat.query("snv=={0}".format(distneeded))

    pair_selected_index = random.sample(list(pairwise_available.index), 1)[0]

    mixed =      pairwise_available.at[pair_selected_index, 'node1']
    mixed_with = pairwise_available.at[pair_selected_index, 'node2']
    

    common_ancestor = bdtree.get_common_ancestor([mixed, mixed_with]).name
    d1 = bdtree.get_distance(mixed, common_ancestor)
    d2 = bdtree.get_distance(mixed_with, common_ancestor)
    dist = snv_distance(simulated_sequences[mixed], simulated_sequences[mixed_with])
    observed_sequences = simulated_sequences
    for key in observed_sequences:
        observed_sequences[key] = introduce_n(observed_sequences[key], prop_random_N)

    print("Exporting observed sequences")        
    observed_sequences[mixed] =  mix_sequences(observed_sequences[mixed], observed_sequences[mixed_with])
    with open(observed_filename,'wt') as f:
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format('node','numN', 'isMixed','mixed_with','common_ancestor','dist_to_mixed','risk_factor_present','seq'))
        for key in sorted(observed_sequences.keys()):
            ns = list(observed_sequences[key]).count('N')
            if key == mixed:
                risk_factor = np.random.binomial(1, prop_risk_mixed)
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(key, ns, key==mixed, mixed_with, common_ancestor, dist,risk_factor,observed_sequences[key]))
            else:
                risk_factor = np.random.binomial(1, prop_risk_not_mixed)
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(key, ns, key==mixed, '', '', 0,risk_factor, observed_sequences[key]))

    # create dictionary with key parameters, to accompany the output
    params = {'n_sequences':n_sequences,
              'len_simulated_alignment':l,
              'prop_random_N':prop_random_N,
              'prop_risk_not_mixed':prop_risk_not_mixed,
              'prop_risk_mixed':prop_risk_mixed,
              'observed_sequence_mixed':mixed,
              'observed_sequence_mixed_with':mixed_with,
              'distance_between_mixed':dist,
              'common_ancestor':common_ancestor}
    with open(params_filename,'wt') as f:
        json.dump(params,f)
    
