""" decorate tree with cluster identifier

# note:
tree is generated on a subset of samples only
subsetting with:
pipenv run python3 demo_covid_testfile.py

tree generation command:
./iqtree2 -s /srv/data/mixfiles/covid/milk_micro.fas -T 6 -m GTR+I --redo

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""

from Bio import Phylo
import pandas as pd
import io


if __name__ == "__main__":
    # load cluster designations
    filename = "v14_eigenvalues.csv"
    ev = pd.read_csv(filename)
    ev.set_index("guid", inplace=True)
    guids = set(ev.index)

    # load tree
    filename = "/srv/data/mixfiles/covid/milk_micro.fas.treefile"
    with (open(filename, "r")) as f:
        Tree = Phylo.read(f, format="newick")

    for i, node in enumerate(Tree.get_terminals()):
        if len(node.name) > 0:
            if node.name in guids:
                node.name = "{0} # {1}".format(node.name, ev.at[node.name, "cluster"])

    with io.StringIO(initial_value="") as f:
        Phylo.write(Tree, f, format="newick")
        retVal = f.getvalue()

    print(retVal)

# tree.write(format=0, outfile='v14_named_tree.nwk')
