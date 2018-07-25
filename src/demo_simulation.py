#!/usr/bin/env python
""" simulates phylogenetic trees.

code is discussed
https://mrnoutahi.com/2017/12/05/How-to-simulate-a-tree/
and taken from
https://gist.githubusercontent.com/maclandrol/09f2b2f47713da9d465a1a7f72078fe7/raw/9874305b9996515fea5d94a90a976c705440fd2e/tree_simul.py

"""

import ete3
import random
import pyvolve
import pandas as pd
import os
import imageio
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









if __name__ == '__main__':

    # H37Rv reference sequence
    print("reading h37rv control sequence")
    inputfile = "../COMPASS_reference/R39/R00000039.fasta"
    with open(inputfile, 'rt') as f:
        for record in SeqIO.parse(f,'fasta', alphabet=generic_nucleotide):               
                seq = str(record.seq)
    miniseq=seq[1:100]
    
    # generate a birth death tree
    print("generating birth-death tree")
    bdtree = birth_death_tree(1,0.7, nsize=30, remlosses=False)

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
    images = []
    ts = ete3.TreeStyle()
    ts.title.add_face(ete3.TextFace("Addition of samples from simulated phylogeny", fsize=8), column=0)
    ts.title.add_face(ete3.TextFace("Gray: not added; Blue: added, not mixed; Red: added& mixed", fsize=8), column=0)
    ts.show_leaf_name = False

    print("Writing tree cluster depiction to AVI")
    avi_filename = os.path.join("..",'output','movie.avi')
    
    with imageio.get_writer(avi_filename, mode='I', fps=3, quality=10) as f:

        for i in range(len(nprops.index)):
            for node in bdtree.traverse():
                
                nstyle = ete3.NodeStyle()
                node_order = nprops.loc[node.name, 'order']
                nstyle['size']=0
                nstyle['fgcolor']='black'
                nstyle['bgcolor']='white'
                if node.is_leaf():
                    nstyle['size']=1
                    nstyle['fgcolor']='gray'
                    nstyle['bgcolor']='white'
                    
                    if node_order < i:
                        nstyle['size']=1
                        nstyle['fgcolor']='blue'
                        nstyle['bgcolor']='white'
    
                node.set_style(nstyle)
            filename = os.path.join("..",'output','{0}.png'.format(i))
            bdtree.render(filename, tree_style= ts, w=928, h=928, units='px')
            img = imageio.imread(filename)
            f.append_data(img)
            os.unlink(filename)
        
    f.close()  

    # simulate evolution from the H37Rv reference
    t = pyvolve.read_tree(tree=bdtree.write(format=3),scale_tree = 0.005)
    m = pyvolve.Model("nucleotide" )
    p = pyvolve.Partition(models=m, root_sequence = miniseq)
    e = pyvolve.Evolver(partitions = p, tree=t)
    
    print("Running evolution")
    e()
    simulated_sequences = e.get_sequences()
    print("Complete")
    for key in sorted(simulated_sequences.keys()):
        print(key,simulated_sequences[key])

    ## TODO: depict cluster assignation next to the tree using a Face
    ## link to a demo fn3 instance
    ## demonstrate detection of mixtures
    ## 
