""" manipulate phylogenetic trees """

import ete3 

class ManipulateTree():
    """ manipulates a newick phylogenetic tree """
    def __init__(self):
        """ initiator.  does nothing """
        pass
    def reroot(self, treestring, root_sample):
        """ sets the root of the newick tree in string treestring
        using root_sample as an outgroup; removes the outgroup

        Returns:
            a phylogenetic tree, rerooted and with root_sample removed """

        tree = ete3.Tree(treestring)
        tree.set_outgroup(root_sample)
        r = tree.search_nodes(name=root_sample)[0]
        r.delete()
    
        return tree

