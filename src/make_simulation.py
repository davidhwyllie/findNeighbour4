#!/usr/bin/env python3
""" produces simulated sequences, storing output to disc.

simulates
* phylogenetic trees, based on a birth death process
* true sequences reflecting the members of the phylogeny,
* sequences with Ns added randomly derived from true sequences (i.e. missing at random, missing data)
* sequences representing mixtures of true sequences.

Also:
* simulates the presence of a risk factor in mixed and unmixed sequences.
* performs single-linkage type clustering based on observed SNV with and without the mixed sample.

"""

import ete3
import random
import pyvolve
import pandas as pd
import os
import numpy as np
import argparse
import json
import copy
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    converts any nucleotides differing between seq1 and seq2 to M  
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
            retVal.append('M')
    return ''.join(retVal)
def introduce_nm(what, seq, p):
    """ introduce character what ('N' or 'M') at random into a sequence, seq.
    The probability that each base is N is p, where 0<=p<=1
    """

    # compute the number of Ns given the length of seq and p
    seqlen = len(seq)
    n = np.random.binomial(seqlen, p)
    
    # compute bases to turn to N
    to_change = random.sample(range(seqlen),n)
    seqlist = list(seq)
    for pos in to_change:
        seqlist[pos] = what
    
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
def seq2distmat(seqdict, output_filename, mixed):
    """ computes a pairwise distance matrix for a dictionary contain sequences
    {'id1':'AAAA', 'id2':'AAAT', ...}
    Write output to output_filename
    Mixed is the identity of a mixed sample.  Pass None if none exists.
    Returns a list containing the distance matrix.
    """
    obs_dist_list = []
    for i in seqdict.keys():
        for j in seqdict.keys():
            if i<j:
               edge = {'node1':i,'node2':j,'snv':snv_distance(seqdict[i], seqdict[j])}
               if not mixed in [i,j]:       # don't link mixed samples.
                   obs_dist_list.append(edge)
  
    df = pd.DataFrame.from_records(obs_dist_list)
    df.to_csv(output_filename)
    return obs_dist_list

def edgelist2graph(node_list, el, snp_cutoff):
    """ converts an edgelist to a networkx graph; returns connected components. """           
    # perform SNV based clustering using mixed and unmixed samples
    G1= nx.Graph()      # mixed included
    nEdges = 0
    node_names= set()
    for item in node_list:
        G1.add_node(item)

    for item in el:
        if item['snv'] <= snp_cutoff:
            nEdges += 1
            G1.add_edge(item['node1'], item['node2'])
    print("Identified {0} edges < {1} ".format(nEdges,snp_cutoff))
    conn_comps = list(nx.connected_components(G1))
    node2cluster = {}
    for i,item in enumerate(conn_comps):
        for node_name in item:
            node2cluster[node_name] = i

    # now make a comparison matrix
    cmpmat = []
    for node_name1 in node2cluster.keys():
        for node_name2 in node2cluster.keys():
            if node_name1 < node_name2:
                item = {'id':"{0}-{1}".format(node_name1, node_name2), 'same_cluster':node2cluster[node_name1]==node2cluster[node_name2]}
                cmpmat.append(item)
    df = pd.DataFrame.from_records(cmpmat, index='id')
    return(df)

def compare_clustering(df1,df2):
    """ compare whether pairs of items in two dataframes are in the same clusters.
    Only pairs of items present in df1 are considered """
    
    df = df1.merge(df2, how='inner', left_index=True, right_index=True)
    def que(x):
        if x['same_cluster_x']==x['same_cluster_y']:
            return 1
        else:
            return 0
    
    return sum(df.apply(que, axis=1))
  

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""
Produce simulated data for testing algorithms and software, storing output to disc.

simulates
* phylogenetic trees, based on a birth death process
* true sequences reflecting the members of the phylogeny,
* sequences with Ns added randomly derived from true sequences (i.e. missing at random, missing data)
* sequences representing mixtures of true sequences.

Also:
* simulates the presence of a risk factor in mixed and unmixed sequences.
* performs single-linkage type clustering based on observed SNV with and without the mixed sample.

Example usage:

python make_simulation.py 50 10000 0.02 0.01 0.001 0.05 0.1 20 1 ../output/simulation_set_1

""")
    
    ####################################################################################
    #  Read inputs from command line.
    ####################################################################################
    
    parser.add_argument('n_sequences', type=int, nargs=1,
                        help='number of sequences to simulate. Suitable value = 10')
    parser.add_argument('len_simulated_alignment', type=int, nargs=1, 
                        help='length of the simulated sequences to generate. Suitable value = 10000')
    parser.add_argument('prop_bases_mutated', type=float, nargs=1,
                        help='proportion of nucleotides which are mutated within the phylogeny.  Must be low to be relevant to TB.  Suitable value = 0.01')
    parser.add_argument('prop_random_N', type=float, nargs=1,
                        help='proportion of nucleotides which are called N by chance.  Suitable value = 0.01')
    parser.add_argument('prop_random_M', type=float, nargs=1,
                        help='proportion of nucleotides which are called M by chance.  Suitable value = 0.001')
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
    prop_random_M = args.prop_random_M[0]
    prop_risk_mixed = args.prop_risk_mixed[0]
    prop_risk_not_mixed = args.prop_risk_not_mixed[0]
    n_simulations = args.n_simulations[0]
    snp_cutoff = args.snp_cutoff[0]


    ####################################################################################
    #  Produce however many simulations are required.
    ####################################################################################
       
    for sim_id in range(n_simulations):
        
        ####################################################################################
        #  Produce one simulation
        ####################################################################################
   
        outputdir = os.path.join(basedir, str(sim_id))
        print("Producing simulation {0} of {1}".format(sim_id, n_simulations))
        
        # set names of output files
        fasta_filename = os.path.join(outputdir,'phylogeny.fasta')
        sequence_filename = os.path.join(outputdir,'phylogeny.txt')
        observed_filename = os.path.join(outputdir,'observed.txt')
        tree_filename = os.path.join(outputdir,'phylogeny.nwk')
        ref_filename = os.path.join(outputdir,'reference.fasta')
        treepic_filename = os.path.join(outputdir,'{0}.png'.format('tree_image'))
        annotated_treepic_filename = os.path.join(outputdir,'{0}.png'.format('annotated_tree_image'))    
        distmat_filename = os.path.join(outputdir,'distmat.txt')
        obs_distmat_pre_filename = os.path.join(outputdir,'obs_distmat_pre.txt')
        obs_distmat_nomix_filename = os.path.join(outputdir,'obs_distmat_nomix.txt')
        obs_distmat_mix_filename = os.path.join(outputdir,'obs_distmat_mix.txt')
        obs_distmat_mix_marked_filename = os.path.join(outputdir,'obs_distmat_mix_marked.txt')
                              
        params_filename = os.path.join(outputdir,'params.json')
          
        # if the output directory does not exist, create it
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        
        # read H37RV (TB) sequence from disc;    
        print("reading h37rv control sequence")
        inputfile = "../COMPASS_reference/R39/R00000039.fasta"
        with open(inputfile, 'rt') as f:
            for record in SeqIO.parse(f,'fasta'):               
                    seq = str(record.seq)
    
        # we only use a subpart of this for our simulations, in part for reasons of speed.
        # these 'minigenomes' comprise the first l bases of the H37RV reference genome.
        print("Making mini genomes of length", l)
        miniseq=seq[0:l]
        with open(ref_filename,'wt') as f:
            f.write(">reference\n{0}\n".format(miniseq))
            
        # generate a birth death tree
        # in the configuration below, there is no death
        print("generating birth-death tree")
        bdtree = birth_death_tree(1,0, nsize=n_sequences, remlosses=False)
    
        # name internal nodes 
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
            
        # store distance to root for all nodes
        for node_name in node_names:           
            dtr = 0
            this_node = node&node_name
            while not this_node.is_root():
                dtr = dtr + this_node.dist
                this_node = this_node.up
            t2n.append({'name':node_name,'dist_to_root':dtr})
        
        # store the distance to root for each node in a dataframe;
        nprops = pd.DataFrame.from_records(t2n, index='name')
        nprops.sort_values('dist_to_root', inplace=True)
        nprops['order'] = range(len(nprops.index))
        
        # depict the tree generated.
        ts = ete3.TreeStyle()
        ts.show_scale = False
    
        ## now the phylogeny is known, we simulate sequences conformant to that phylogeny
        # starting from the H37Rv reference minigenome
        # we want at most 1-3% of our sequence to vary.  This is because we want to
        # have a relevant (and low) changce of back mutations (i.e. two changes at the same site)
        # as these are very rare for the organism studied.
        
        # read tree and determine root to tip distance
        max_rtt = max([x['dist_to_root'] for x in t2n])
        scaling_factor = prop_bases_mutated/max_rtt      
        for node in bdtree.traverse():           
            node.dist = node.dist * scaling_factor
        bdtree.write(outfile = tree_filename, format=3)

        # now we use the pyvolve module, see http://sjspielman.org/pyvolve/
        # Spielman, SJ and Wilke, CO. 2015.
        # Pyvolve: A flexible Python module for simulating sequences along phylogenies. PLOS ONE. 10(9): e0139047.
        t = pyvolve.read_tree(tree=bdtree.write(format=3))  
        m = pyvolve.Model("nucleotide" )
        p = pyvolve.Partition(models=m, root_sequence = miniseq)
        e = pyvolve.Evolver(partitions = p, tree=t)
     
        # Run evolution
        e()
        
        # Recover sequences from the evolution;
        # write output to file.
        simulated_sequences = e.get_sequences()
        with open(sequence_filename,'wt') as f:
            for key in sorted(simulated_sequences.keys()):
                f.write("{0}\t{1}\n".format(key,simulated_sequences[key]))
        with open(fasta_filename,'wt') as f:
            for key in sorted(simulated_sequences.keys()):
                f.write(">{0}\n{1}\n".format(key,simulated_sequences[key]))
    
        # compute pairwise matrix for sequences observed in the absence of mixtures or Ns
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
    
        #  depict tree comprising the observed sequences
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
    
        
        # now generate a simulated mixture    
        print("Generating mixture")
        # select two nodes at random from the available distances
        # this uses stratified sampling, ensuring that all parts of the snv distribution are sampled
        distneeded = random.sample(dists, 1)
        pairwise_available = distmat.query("snv=={0}".format(distneeded))
    
        pair_selected_index = random.sample(list(pairwise_available.index), 1)[0]
    
        mixed =      pairwise_available.at[pair_selected_index, 'node1']
        mixed_with = pairwise_available.at[pair_selected_index, 'node2']
        
        # compute the distance to common ancestor between the component sequences of the mixture
        common_ancestor = bdtree.get_common_ancestor([mixed, mixed_with]).name
        d1 = bdtree.get_distance(mixed, common_ancestor)
        d2 = bdtree.get_distance(mixed_with, common_ancestor)
        mixed_dist = snv_distance(simulated_sequences[mixed], simulated_sequences[mixed_with])
        
        # compute observed sequences, including random Ns to reflect sequencing error
        observed_sequences_nomix = copy.copy(simulated_sequences)
        for key in observed_sequences_nomix.keys():
            observed_sequences_nomix[key] = introduce_nm('N', observed_sequences_nomix[key], prop_random_N)
            observed_sequences_nomix[key] = introduce_nm('M', observed_sequences_nomix[key], prop_random_M)
        observed_sequences_mix = copy.copy(observed_sequences_nomix)
        observed_sequences_mix[mixed] =  mix_sequences(observed_sequences_nomix[mixed], observed_sequences_nomix[mixed_with])
       
        print("Computing pairwise distance matrix for observed data with mixed  and unmixed data")
        observed_sequences_without_mix = copy.copy(observed_sequences_nomix)
        del observed_sequences_without_mix[mixed]
        obs_distmat_preadd      = seq2distmat(observed_sequences_without_mix, None, obs_distmat_pre_filename) # with addition of unmixed      
        obs_distmat_nomix      = seq2distmat(observed_sequences_nomix, None, obs_distmat_nomix_filename) # with addition of unmixed
        obs_distmat_mix        = seq2distmat(observed_sequences_mix, None, obs_distmat_mix_filename) # with addition of mixed
        obs_distmat_mix_marked = seq2distmat(observed_sequences_mix, mixed, obs_distmat_mix_marked_filename) # with poss. mixed ignored
    
        connections = {}
        connections['0_preadd'] =     edgelist2graph(observed_sequences_without_mix.keys(), obs_distmat_preadd, snp_cutoff)
        connections['1_nomix'] =      edgelist2graph(observed_sequences_nomix.keys(), obs_distmat_nomix, snp_cutoff)
        connections['2_mix'] =        edgelist2graph(observed_sequences_mix.keys(), obs_distmat_mix, snp_cutoff)
        connections['3_mix_marked'] = edgelist2graph(observed_sequences_mix.keys(), obs_distmat_mix_marked, snp_cutoff)

        comparison = {
            'sim_id':sim_id,
            'nEdges':len(connections['0_preadd'].index),
            'dist2ancestor':(d1+d2) * l,
            'mix_components_snv':mixed_dist,
        }
        
        for element in ['0_preadd','1_nomix','2_mix','3_mix_marked']:
            comparison[element] = compare_clustering(connections['0_preadd'], connections[element])
            

        print("Exporting observed sequences and associated tree")        
        with open(observed_filename,'wt') as f:
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format('node','numN', 'isMixed','mixed_with','common_ancestor','dist_to_mixed','risk_factor_present','seq','cluster_all','cluster_no_mix'))
            for key in sorted(observed_sequences_nomix.keys()):
                
                # determine what clustering groups it is in
                ## TODO
                #for element in   ['0_preadd','1_nomix','2_mix','3_mix_marked']:
                #    if key in comparison[element].index:
                #        cln = comparison[element].loc[key,]
                
                # annotate tree
                subtree = bdtree&key
                #subtree.add_face(ete3.TextFace(cl_number_1,ftype='Courier New', fsize=6), column=2, position='aligned')
                #subtree.add_face(ete3.TextFace(cl_number_2,ftype='Courier New', fsize=6), column=3, position='aligned')
                #subtree.add_face(ete3.TextFace(observed_sequences[key], ftype='Courier New', fsize=6), column=0, position='aligned')
                #subtree.add_face(ete3.TextFace('-', ftype='Courier New', fsize=6), column=0, position='aligned')
                
                ns = list(observed_sequences_mix[key]).count('N')
                ms = list(observed_sequences_mix[key]).count('M')
                subtree.add_face(ete3.TextFace(ns, ftype='Courier New', fsize=6), column=0, position='aligned')
                subtree.add_face(ete3.TextFace(ms, ftype='Courier New', fsize=6), column=1, position='aligned')
                if key == mixed:
                    subtree.add_face(ete3.TextFace('Y c. {0}'.format(mixed_with), ftype='Courier New', fsize=6), column=2, position='aligned')
                    risk_factor = np.random.binomial(1, prop_risk_mixed)
                    #f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(key, ns, key==mixed, mixed_with, common_ancestor, dist,risk_factor,observed_sequences_mix[key], cl_number_1, cl_number_2))
                elif not key == mixed_with:         # we don't export mixed_with ??
                    # optionally don't export a proportion ( ? 30%)
                    subtree.add_face(ete3.TextFace('N',ftype='Courier New', fsize=6), column=1, position='aligned')
                    risk_factor = np.random.binomial(1, prop_risk_not_mixed)
                    #f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(key, ns, key==mixed, '', '', 0,risk_factor, observed_sequences_mix[key], cl_number_1, cl_number_2))
    
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
                  'distance_between_mixed':mixed_dist,
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
                    
                    