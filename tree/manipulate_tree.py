""" manipulate phylogenetic trees """

import ete3


class ManipulateTree:
    """manipulates a newick phylogenetic tree

    provides an interface to functions largely provided by ete3  http://etetoolkit.org/
    self.tree is an ete3 object, and ete3 methods can be called against this"""

    def __init__(self, treestring):
        """constructs a ManipulateTree object starting with a Newick string

        parameters:
        treestring: a newick tree

        returns: none"""
        self.tree = ete3.Tree(treestring)
        self.outgroup = None

    def nodes(self):
        """returns the nodes of the tree in order of a postorder traverse"""
        node_list = []
        for x in self.tree.traverse("postorder"):
            node_list.append(x)
        return x

    def remove_outgroup(self):
        """removes the outgroup from the tree"""
        if self.outgroup is not None:
            (self.tree & self.outgroup).detach()
            self.outgroup = None
            # r=self.tree.get_tree_root()
            # for rootchild in r.children:
            #        rootchild.dist=0

    def rtd(self):
        """compute root tip distance"""
        # find the root
        r = self.tree.get_tree_root()
        # and the farthest distance from the root
        farthest, dist = r.get_farthest_node()
        return dist

    def reroot(self, root_sample):
        """sets the root of the tree
        using root_sample as an outgroup

        parameters:
        root_sample: the root node

        Returns:
            a phylogenetic tree, rerooted and with root_sample removed"""

        self.tree.set_outgroup(root_sample)
        self.outgroup = root_sample
        return self.tree
