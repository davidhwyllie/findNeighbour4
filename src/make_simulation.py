#!/usr/bin/env python3
""" produces simulated sequences, storing output to disc.

simulates phylogenetic trees,
sequences reflecting the members of the phylogeny,
sequences with Ns added randomly derived from the phylogeny,
and mixtures of sequences from within the phylogeny.
Also:
* simulates the presence of a risk factor in mixed and unmixed sequences.
* performs clustering based on observed SNV with and without the mixed sample.

"""

import ete3
import random
import pyvolve
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

def delete_single_child_internal(t):
    """Utility function that removes internal nodes
    with a single child from tree"""

    for node in t.traverse("postorder"):
        if(not node.is_leaf() and len(node.get_children()) < 2):
            node.delete()

    if len(t.get_children()) == 1:
        t = t.children[0]
        t.up = None
def birth_death_tree(birth, death, nsize=10, max_time=None, remlosses=True, r=True):
    """Generates a birth-death tree.
    Arguments:
        - ``birth`` : birth rate
        - ``death`` : death rate
        - ``nsize`` : desired number of leaves
        - ``max_time`` : maximum time of evolution
        - ``remlosses`` : whether lost leaves (extinct taxa) should be pruned from tree
        - ``r`` : repeat until success
        
    Approach is discussed here:  https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/
    Code is modified from:
    https://gist.githubusercontent.com/maclandrol/09f2b2f47713da9d465a1a7f72078fe7/raw/9874305b9996515fea5d94a90a976c705440fd2e/tree_simul.py

    """
    # initialize tree with root node
    tree = ete3.Tree()
    tree.add_features(extinct=False)
    tree.dist = 0.0

    total_time = 0
    died = set([])

    # total event rate to compute waiting time
    event_rate = float(birth + death)

    while True:
        # waiting time based on event_rate
        wtime = random.expovariate(event_rate)
        total_time += wtime
        
        leaf_nodes = tree.get_leaves()
        curr_num_leaves = len(leaf_nodes)
        for leaf in leaf_nodes:
            # extinct leaves cannot update their branches length
            if not leaf.extinct:
                leaf.dist += wtime

        if curr_num_leaves >= nsize:
            break

        # if event occurs within time constraints
        if max_time is None or total_time <= max_time:

            # select node at random, then find chance it died or gave birth
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
    to_change = random.sample(range(seqlen),n)
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
Example usage:

python make_simulation.py 50 10000 0.02 0.01 0.05 0.1 20 1 ../output/simulation_set_1

""")
    parser.add_argument('n_sequences', type=int, nargs=1,
                        help='number of sequences to simulate. Suitable value = 10')
    parser.add_argument('len_simulated_alignment', type=int, nargs=1, 
                        help='length of the simulated sequences to generate. Suitable value = 15')
    parser.add_argument('prop_bases_mutated', type=float, nargs=1,
                        help='proportion of nucleotides which are mutated.  Must be low to be relevant to TB.  Suitable value = 0.01')
    parser.add_argument('prop_random_N', type=float, nargs=1,
                        help='proportion of nucleotides which are called N by chance.  Suitable value = 0.01')
    parser.add_argument('prop_risk_not_mixed', type=float, nargs=1,
                        help='proportion of non-mixed samples with a risk factor. Suitable value = 0.05')
    parser.add_argument('prop_risk_mixed', type=float, nargs=1,
                        help='proportion of mixed samples with a risk factor. Suitable value = 0.10')    
    parser.add_argument('snp_cutoff', type=int, nargs=1,
                        help='SNV clustering cutoff, for illustrating impact of mixed samples.  Suitable value = 5-10')
    parser.add_argument('n_simulations', type=int, nargs=1,
                        help='number of simulations performed.  Suitable value = 5')
    parser.add_argument('outputdir', type=str, nargs=1,
                        help='output will be written to the outputdir')    
    args = parser.parse_args()
    n_sequences = args.n_sequences[0]
    basedir = os.path.abspath(args.outputdir[0])
    l = args.len_simulated_alignment[0]
    prop_bases_mutated = args.prop_bases_mutated[0]
    prop_random_N = args.prop_random_N[0]
    prop_risk_mixed = args.prop_risk_mixed[0]
    prop_risk_not_mixed = args.prop_risk_not_mixed[0]
    n_simulations = args.n_simulations[0]
    snp_cutoff = args.snp_cutoff[0]
    
    for sim_id in range(n_simulations):
        outputdir = os.path.join(basedir, str(sim_id))
        print("Producing simulation {0} of {1}".format(sim_id, n_simulations))
        # names of output files
        fasta_filename = os.path.join(outputdir,'phylogeny.fasta')
        sequence_filename = os.path.join(outputdir,'phylogeny.txt')
        observed_filename = os.path.join(outputdir,'observed.txt')
        tree_filename = os.path.join(outputdir,'phylogeny.nwk')
        ref_filename = os.path.join(outputdir,'reference.fasta')
        treepic_filename = os.path.join(outputdir,'{0}.png'.format('tree_image'))
        annotated_treepic_filename = os.path.join(outputdir,'{0}.png'.format('annotated_tree_image'))
        
        distmat_filename = os.path.join(outputdir,'distmat.txt')
        obs_distmat_filename = os.path.join(outputdir,'obs_distmat.txt')
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
        miniseq=seq[0:l]
        with open(ref_filename,'wt') as f:
            f.write(">reference\n{0}\n".format(miniseq))
            
        # generate a birth death tree
        print("generating birth-death tree")
        bdtree = birth_death_tree(1,0, nsize=n_sequences, remlosses=False)
    
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
        ts.show_scale = False
    
        ## simulate evolution from the H37Rv reference
        
        # read tree and determine root to tip distance
        max_rtt = max([x['dist_to_root'] for x in t2n])
        
        # we want at most 1-3% of our sequence to vary.  This is because we want to
        # minimise the chance of back mutations, as these are very rare in our application
        scaling_factor = prop_bases_mutated/max_rtt
        
        for node in bdtree.traverse():           
            node.dist = node.dist * scaling_factor
        bdtree.write(outfile = tree_filename, format=3)

        t = pyvolve.read_tree(tree=bdtree.write(format=3))  
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
    
        print("Computing pairwise distance matrix for true data")
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
                
                if node.is_leaf():
                    if node_order < i:
                        
                        nstyle['size']=1
                        nstyle['fgcolor']='blue'
                        nstyle['bgcolor']='white'                        
    
                node.set_style(nstyle)
        bdtree.render(treepic_filename, tree_style= ts,units='px', w=500, h=1024) 
    
            
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
    

        print("Computing pairwise distance matrix for observed data")
        obs_dist_list = []
        obs_dist_list_nomix = []
        for i in observed_sequences.keys():
            for j in observed_sequences.keys():
                if i<j:
                   edge = {'node1':i,'node2':j,'snv':snv_distance(observed_sequences[i], observed_sequences[j])}
                   obs_dist_list.append(edge)
                   if not (i == mixed or j == mixed):
                       obs_dist_list_nomix.append(edge)
        
        obs_distmat = pd.DataFrame.from_records(obs_dist_list)
        obs_distmat.to_csv(obs_distmat_filename)
        
        # perform SNV based clustering using mixed and unmixed samples
        G1= nx.Graph()      # mixed included
        G2= nx.Graph()      # mixed not included
        nEdges = 0
        for key in sorted(observed_sequences.keys()):
            G1.add_node(key)
            G2.add_node(key)
        for item in obs_dist_list:
            if item['snv'] <= snp_cutoff:
                nEdges += 1
                G1.add_edge(item['node1'], item['node2'])
        print("Identified {0} edges < {1} ".format(nEdges,snp_cutoff))
        for item in obs_dist_list_nomix:
            if item['snv'] <= snp_cutoff:
                G2.add_edge(item['node1'], item['node2'])
        g1_components = list(nx.connected_components(G1))
        g2_components = list(nx.connected_components(G2))
        
        print("Exporting observed sequences and associated tree")        
        observed_sequences[mixed] =  mix_sequences(observed_sequences[mixed], observed_sequences[mixed_with])
        with open(observed_filename,'wt') as f:
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format('node','numN', 'isMixed','mixed_with','common_ancestor','dist_to_mixed','risk_factor_present','seq','cluster_all','cluster_no_mix'))
            for key in sorted(observed_sequences.keys()):
                
                # determine what clustering groups it is in
                cl_number_1 =0
                for cl in g1_components:
                    cl_number_1 +=1
                    if key in cl:
                        break
                cl_number_2 =0
                if not key == mixed:
                    for cl in g2_components:
                        cl_number_2 +=1
                        if key in cl:
                            break
                
                # annotate tree
                subtree = bdtree&key
                subtree.add_face(ete3.TextFace(cl_number_1,ftype='Courier New', fsize=6), column=2, position='aligned')
                subtree.add_face(ete3.TextFace(cl_number_2,ftype='Courier New', fsize=6), column=3, position='aligned')
                #subtree.add_face(ete3.TextFace(observed_sequences[key], ftype='Courier New', fsize=6), column=0, position='aligned')
                #subtree.add_face(ete3.TextFace('-', ftype='Courier New', fsize=6), column=0, position='aligned')
                
                ns = list(observed_sequences[key]).count('N')
                subtree.add_face(ete3.TextFace(ns, ftype='Courier New', fsize=6), column=0, position='aligned')
                if key == mixed:
                    subtree.add_face(ete3.TextFace('c. {0}'.format(mixed_with), ftype='Courier New', fsize=6), column=2, position='aligned')
                    risk_factor = np.random.binomial(1, prop_risk_mixed)
                    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(key, ns, key==mixed, mixed_with, common_ancestor, dist,risk_factor,observed_sequences[key], cl_number_1, cl_number_2))
                else:
                    subtree.add_face(ete3.TextFace('No',ftype='Courier New', fsize=6), column=1, position='aligned')
                    risk_factor = np.random.binomial(1, prop_risk_not_mixed)
                    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(key, ns, key==mixed, '', '', 0,risk_factor, observed_sequences[key], cl_number_1, cl_number_2))
    
        # export tree depicting clustering with mixed and unmixed samples
        bdtree.render(annotated_treepic_filename, tree_style= ts,units='px', w=1024, h=1024) 
    
        # create dictionary with key parameters, to accompany the output
        print("Exporting parameters")
        params = {'n_sequences':n_sequences,
                  'len_simulated_alignment':l,
                  'prop_bases_mutated':prop_bases_mutated,
                  'prop_random_N':prop_random_N,
                  'prop_risk_not_mixed':prop_risk_not_mixed,
                  'prop_risk_mixed':prop_risk_mixed,
                  'observed_sequence_mixed':mixed,
                  'observed_sequence_mixed_with':mixed_with,
                  'distance_between_mixed':dist,
                  'common_ancestor':common_ancestor}
        with open(params_filename,'wt') as f:
            json.dump(params,f)
        
        # produce little SNV matrices for each cluster identified.
        print("Exporting distance matrices")
        obs = pd.read_csv(observed_filename, sep='\t')
        for clustering_method in ['cluster_all','cluster_no_mix']:
            clusters = obs[clustering_method].unique()
            for this_cluster in clusters:

                cl_data = obs.query("{0}=={1}".format(clustering_method,this_cluster))
                node_names = cl_data['node'].unique()
                mat = obs_distmat.query("node1 in @node_names")
                if len(mat.index)>0:
                    mat =     mat.query("node2 in @node_names")
                if len(mat.index)>0:
                    piv = mat.pivot(index='node1', columns='node2', values='snv')
                    clusterdir = os.path.join(outputdir,'clusters')
                    if not os.path.exists(clusterdir):
                         os.makedirs(clusterdir)
                    pivot_output_file =  os.path.join(outputdir,'clusters','{0}_{1}.xlsx'.format(clustering_method, this_cluster))
                    piv.to_excel(pivot_output_file)
                    