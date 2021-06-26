""" manipulate phylogenetic trees """
import os
import ete3
import logging
import subprocess
import shutil
import datetime
import copy

# import pandas as pd


class IQTree:
    """builds a tree with iqTree

    Iqtree: http://www.iqtree.org/

    """

    def __init__(self, genome_length, use_bootstrap=False):
        """checks whether iqTree is present"""
        if os.getenv("IQTREE_DIR") is not None:
            if not os.path.exists(os.environ["IQTREE_DIR"]):
                raise FileNotFoundError(
                    "no directory at the location specified by environment variable IQTREE_DIR: {0}".format(
                        os.environ["IQTREE_DIR"]
                    )
                )

            self.iqTreeLoc = os.environ["IQTREE_DIR"]
            logging.info("Using iqTree version at {0}")
            self.iqTreeExe = "iqtree2"

            if use_bootstrap is True:
                self.bootstrapOption = "-bb 1000"
            else:
                self.bootstrapOption = ""

            self.min_branch_length = 0.1 * (1 / genome_length)
        else:
            self.iqTreeExe = None

    def build(self, fasta_file_string, fconst, targetdir):
        """builds a maximal likelihood tree using iqtree2

        Parameters:
        fasta_file_string:  the fasta file
        fconst: a Counter object or dictionary containing the counts of the constant bases
        """

        if self.iqTreeExe is None:
            return None
            
        start_time = datetime.datetime.now()

        # check the fconst string.
        if not isinstance(fconst, dict):
            raise TypeError("fconst must be a dictionary, not {0}".format(type(fconst)))

        essential_keys = ["A", "C", "G", "T"]
        for essential_key in essential_keys:
            if essential_key not in fconst.keys():
                fconst[essential_key] = 0
        base_freq_string = "{0},{1},{2},{3}".format(
            fconst["A"], fconst["C"], fconst["G"], fconst["T"]
        )

        # delete the directory if it exists
        if os.path.exists(targetdir):
            shutil.rmtree(targetdir)

        # make the target dir if it does not exist
        os.makedirs(targetdir, exist_ok=True)
        fastafile = os.path.join(targetdir, "alignment.fa")
        with open(fastafile, "wt") as f:
            f.write(fasta_file_string)

        cmd = " ".join(
            [
                os.path.join(self.iqTreeLoc, self.iqTreeExe),
                " -nt 8 -st DNA -s {0}".format(self.bootstrapOption),
                os.path.join(os.path.abspath(targetdir), "alignment.fa"),
                "-m GTR+I",
                "-blmin {0}".format(self.min_branch_length),
                "-fconst",
                base_freq_string,
            ]
        )

        logging.info("iqTree command: {0}".format(cmd))

        try:
            subprocess.check_output(
                cmd,
                shell=True,
                cwd=os.path.abspath(targetdir),
                stderr=subprocess.STDOUT,
            )
            success = True
        except subprocess.CalledProcessError as exc:

            # it crashed
            logging.error("iqtree execution failed with cmd = " + cmd)
            logging.warning(exc.output.decode("utf-8"))
            with os.path.join(
                os.path.abspath(targetdir), "iqTree_error_message.err"
            ) as f:
                f.write(exc.output)
            success = False

        end_time = datetime.datetime.now()

        expected_results = {
            "alignment": "alignment.fa",
            "report": "alignment.fa.iqtree",
            "log": "alignment.fa.log",
            "newick": "alignment.fa.treefile",
        }

        retVal = {
            "success": success,
            "start_time": start_time,
            "end_time": end_time,
            "output": {"build_command": cmd},
        }
        for retVal_item in expected_results.keys():
            targetfile = os.path.join(targetdir, expected_results[retVal_item])
            if os.path.exists(targetfile):
                with open(targetfile, "rt") as f:
                    content = f.read()
                    retVal["output"][retVal_item] = content
        return retVal

        return (success, targetdir)


class DepictTree:
    """depicts a newick phylogenetic tree, attaching metadata"""

    def __init__(
        self, treestring, metadata, title=["Default title"], genome_length=29903
    ):
        """constructs a DepictTree object starting with a Newick string and metadata

        parameters:
        treestring: a newick tree in a character string
        metadata: a pandas dataframe, where the index are the nodes and the columns are for depiction.

        returns: none
        """

        self._refresh()
        self.genome_length = genome_length
        self.tree = ete3.Tree(treestring)

        # zero small branch lengths (< 0.5 SNP)
        self.zero_small_branch_lengths()

        # only display node names & dates for expanding clones - the rest is just a random recent sample for context
        self.expanding = metadata["expanding"]
        displayed_sample_ids = []
        displayed_sample_dates = []
        for ix in metadata.index:
            if metadata.at[ix, "expanding"] == "Y":
                displayed_sample_ids.append(ix)
                displayed_sample_dates.append(metadata.at[ix, "sample_date"])

            else:
                displayed_sample_ids.append(" ")
                displayed_sample_dates.append(" ")
        metadata = metadata.drop(columns=["expanding"])
        cols = metadata.columns.to_list()
        metadata["sample_id"] = metadata.index
        cols.insert(0, "sample_id")
        self.metadata = metadata[cols]
        self.metadata["sample_id"] = displayed_sample_ids
        self.metadata["sample_date"] = displayed_sample_dates

        self.title = title
        # check: are the nodes and the index of self.metadata the same?
        nodes = set(self.nodes())

        metadata_ix = set(self.metadata.index.to_list())
        if not nodes == metadata_ix:
            raise ValueError(
                "Nodes and metadata index do not match: only in nodes{0} vs only in metadata index {1}".format(
                    nodes - metadata_ix, metadata_ix - nodes
                )
            )

        # add all attributes in metadata to the nodes.
        for node in self.tree.traverse("postorder"):
            # fc = ete3.FaceContainer()
            for i, col in enumerate(self.metadata.columns.to_list()):
                # node.add_feature(col,self.metadata.loc[node.name, col])        # add the whole dictionary

                if len(node.name) > 0:
                    nstyle = ete3.NodeStyle()
                    nstyle["shape"] = "circle"
                    nstyle["size"] = 5
                    nstyle["fgcolor"] = "gray"
                    if self.expanding[node.name] == "Y":
                        nstyle["fgcolor"] = "red"

                    node.set_style(nstyle)

                    f = ete3.TextFace(self.metadata.loc[node.name, col], fsize=6)
                    node.add_face(f, position="aligned", column=i)

    def _refresh(self):
        """get ready for next tree"""
        self.genome_length = None
        self.tree = None
        self.branchSupport = None
        self._annotations = None
        self.outgroup = None
        self.has_title = False
        self.rtd_displayed = False
        self.distance_unit = "SNV"
        self.title = []
        self.show_leaf_name = False
        self.show_branch_length = False
        self.show_branch_support = False
        self._headers = {}  # tree column headings
        self._tss = ete3.TreeStyle()  # global configuration parameters for the tree
        self.page_units = "mm"
        self.page_height = 297
        self.page_width = 210
        self.dpi = 300

    def _prepareTreeStyle(self):
        """sets options for tree setting"""
        self._tss = ete3.TreeStyle()
        self._tss.show_leaf_name = self.show_leaf_name
        self._tss.show_branch_length = self.show_branch_length
        self._tss.show_branch_support = self.show_branch_support
        scaleInfo = ""
        self._tss.show_scale = True
        if self.show_branch_support is True or self.show_branch_length is True:
            self._tss.show_scale = False

        ## setting scale.  Scale is pixels per branch length unit.
        # at 300 dpi, and a tree of 6 inches across, that's ~ 2000 pixels.
        # empirically, about 1000 works better
        # however, screen display can run at 72dpi
        # the screen is 15 inches across
        # that's about 1000 px
        # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html
        # suppose we want 1000 pixels to be the width of the tree
        # rtd is the maximal root-tip distance
        # so the scale is set as below.

        dist = self.rtd()
        if dist == 0:
            self._tss.show_scale = False
        if dist <= 5:
            dist = 5
        if dist > 400:
            dist = 400
        scale_unit = float(int(400 / dist))
        # self._tss.mode = 'c'
        self._tss.scale = scale_unit
        self._tss.optimal_scale_level = "full"
        self._tss.mode = "r"
        root_to_tip = self.rtd()
        if self.genome_length is None:
            self.distance_unit = "Subs/ Kb"
            root_to_tip = root_to_tip * 1000
        else:
            self.distance_unit = "SNV"
            root_to_tip = root_to_tip * self.genome_length

        scaleInfo = "Root-to-tip distance is {0:.1f} {1}.".format(
            root_to_tip, self.distance_unit
        )
        pca_info = "For expanding branches, the following are shown: {0}".format(
            ";".join(self.metadata.columns.to_list())
        )
        if self.title is not None:
            nrows = 0
            self._tss.title = ete3.treeview.main.FaceContainer()
            for titlerow in self.title:
                nrows += 1
                self.hasTitle = True
                if nrows == 1:
                    fs = 14
                else:
                    fs = 8
                # write the information provided, row by row
                self._tss.title.add_face(ete3.TextFace(titlerow, fsize=fs), column=0)
            # write the scale informtion
            self._tss.title.add_face(ete3.TextFace(scaleInfo, fsize=6), column=0)
            # add information about the columns
            self._tss.title.add_face(ete3.TextFace(pca_info, fsize=6), column=0)
            control_info = "Gray samples represent a random selection of recent samples which do not belong to the expanding categories"
            self._tss.title.add_face(ete3.TextFace(control_info, fsize=6), column=0)

    def zero_small_branch_lengths(self):
        """zeros small branch lengths less than 0.5 SNP equivalent"""
        for node in self.tree.traverse("postorder"):
            if node.dist < (0.5 / self.genome_length):
                node.dist = 0

    def rtd(self):
        """compute root tip distance"""
        # find the root
        r = self.tree.get_tree_root()
        # and the farthest distance from the root
        farthest, dist = r.get_farthest_node()
        return dist

    def nodes(self):
        """returns the nodes of the tree in order of a postorder traverse"""
        node_list = []
        for x in self.tree.traverse("postorder"):
            if len(x.name) > 0:
                node_list.append(x.name)
        return node_list

    def render(self, outputfile, mode="r"):
        """render the tree to file"""
        self._tss.mode = mode
        self._prepareTreeStyle()
        if os.path.exists(outputfile):
            os.unlink(outputfile)  # get rid of it if exists
        self.tree.render(outputfile, tree_style=self._tss)
        return outputfile


class ManipulateTree:
    """manipulates  a newick phylogenetic tree

    provides an interface to functions largely provided by ete3  http://etetoolkit.org/
    self.tree is an ete3 object, and ete3 methods can be called against this

    """

    def __init__(self, treestring, genome_length=29903):
        """constructs a ManipulateTree object starting with a Newick string

        parameters:
        treestring: a newick tree
        genome_length: the length of the refernece genome to which the sequences making up the tree were mapped

        returns: none"""
        self.tree = ete3.Tree(treestring)
        self.outgroup = None
        self.genome_length = genome_length
        self.stree = {}

    def nodes(self):
        """returns the nodes of the tree in order of a postorder traverse"""
        node_list = []
        for x in self.tree.traverse("postorder"):
            if len(x.name) > 0:
                node_list.append(x.name)
        return node_list

    def remove_outgroup(self):
        """removes the outgroup from the tree"""
        if self.outgroup is not None:
            (self.tree & self.outgroup).detach()
            self.outgroup = None

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
        r = self.tree.get_tree_root()
        r.dist = 0  # branch 'up from' root is length zero
        return self.tree

    def newick(self, format=2):
        """return the tree in newick format

        Parameters:
        format: format to export in, see http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#writing-newick-trees"""
        return self.tree.write(format=format)

    def _find_proximal_home(self, t, node_name, forbidden=set()):
        """moves towards the root of tree t from a node until it finds a node
        with a non-zero distance to its parent, or the root.
        it ignores any nodes with names in the set forbidden.
        It does this so that the 'new parent' of any node which
        is about to be removed using a prune() method can be determined.
        An internal function called during blobogram construction."""
        # find the node referred to

        this_node = t & node_name

        while this_node.dist == 0:  # proceed unless it's the root
            # logging.info("At {0}".format(this_node.name))
            if this_node.is_root():
                return this_node.name  # the root

            # logging.info("Moving up from", this_node.name)
            this_node = this_node.up  # the parent look further up
        return this_node.name

    def zero_small_branch_lengths(self):
        """zeros small branch lengths less than 0.5 SNP equivalent"""
        for node in self.tree.traverse("postorder"):
            if node.dist < (0.5 / self.genome_length):
                node.dist = 0

    def generate_simplified_tree(self):
        """collapses any zero length branches creating a new tree data structure, self.stree, which is visualisable using Bokeh
        with functions below"""

        self.stree = {
            "originalTree": copy.copy(self.tree)
        }  # a data structure which will contain inputs and outputs

        #  the deepcopy is critical

        t = copy.deepcopy(self.tree)  
        # label the internal nodes, mainly so we can monitor what happens next
        initial_leaves = 0
        nInternal = 0
        for node in t.traverse():

            if not node.is_leaf():
                nInternal += 1
                node.name = "_{0}".format(nInternal)
            else:
                initial_leaves += 1

        print(
            "Now simplifying the tree, removing nodes with zero distance to their parent"
        )

        ninitial = 0
        for node in t.traverse():
            if not node.is_root():
                if node.dist == 0:
                    ninitial += 1
        
        # prune
        to_prune = set()
        initial_nodes = set()
        for node in t.traverse():
            initial_nodes.add(node.name)
            if node.dist == 0 and not node.is_root():
                to_prune.add(node.name)

        # iterate over the nodes.  if the distance (up, towards the parent, which is node.dist in ete3 speak
        # zero, we want to remove this and merge its children with node's surviving parent.
        # intially, we use a very fast ete3 method which results in a dichtomous tree.
        node_parents = dict()
        for node in t.traverse():
            if node.is_leaf():  # we are only interested in the leaves
                node_parents[node.name] = self._find_proximal_home(
                    t, node.name, forbidden=to_prune
                )

        remaining = initial_nodes - to_prune
        forbidden = set()
        try:
            t.prune(list(remaining), preserve_branch_length=True)
            for node in t.traverse():
                if node.dist > 0 and not node.is_root():
                    forbidden.add(node.name)
        except ete3.coretype.tree.TreeError:  # raise when the tree is very small and all would be pruned
            logging.warning("Tree prune failed")

        # now we iterate again, removing any samples where a->b->c and a->b = 0 and b->c=0

        for node in t.traverse():
            if node.dist == 0 and not node.is_root():
                node.delete()

        # assemble a dictionary -> list structure describing the names of the nodes in the old tree to which the new tree refers
        logging.info(
            "Building dictionary mapping new node names to a set old node names"
        )
        node_contents = (
            {}
        )  # keys refer to the name of the new tree; contents to the names of the old tree nodes
        # logging.info("Building node contents")
        for node in t.traverse():
            node_contents[node.name] = []  # an empty list

        # node_parents is derived from the original tree.
        # the keys of node_parents are the original node names.
        # the values are the new node names
        final_leaves = 0
        for original_node_name in node_parents.keys():
            original_node_parent_name = node_parents[original_node_name]
            if (
                original_node_parent_name in node_contents.keys()
            ):  # if it still exists in the new tree
                node_contents[original_node_parent_name].append(original_node_name)
                final_leaves += 1
            else:
                logging.warn(
                    "Node no longer exists: {0}".format(original_node_parent_name)
                )

        # sanity check: ensure we have not lost any leaves in our simplification
        if not (initial_leaves == final_leaves):
            logging.critical(
                "Lost leaves during the transformation: {0}; {1}; {2}".format(
                    initial_leaves, final_leaves, node_contents
                )
            )
        assert initial_leaves == final_leaves

        # store results in the self.stree dictionary
        self.stree["node_contents"] = node_contents
        self.stree["simplifiedTree"] = t
        self.stree["rnn"] = node_contents.keys()

        #for node_name in node_contents.keys():
        #    print(node_name, len(node_contents[node_name]), node_contents[node_name])
        return ()
