#!/usr/bin/env python
""" returns a hash on a list.  Used for identifying clusters

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 
 """

import hashlib


class IdentifySequenceSet:
    """generates a unique key for a given set of sequence identifiers"""

    def __init__(self):
        """create the object"""
        self.permitted_collection_types = set(
            ["cluster", "connections", "msa", "distmat", "njtree", "iqtree", "raxml"]
        )
        self.has_outgroup2label = {True: "og", False: "no_og"}
        self.permitted_what = set(["M", "N", "N_or_M", "-"])

    def _hashComponents(self, x: list) -> str:
        """returns an sha1 hash on a list, x"""
        x = sorted(x)

        to_hash = ";".join([str(item) for item in x])
        return hashlib.sha1(to_hash.encode("utf-8")).hexdigest()

    def make_identifier(self, collection_type, what, has_outgroup, x):
        """returns id for a representation of a collection of samples.
        collection_type can be one of the following:
            msa, distmat, njtree, iqtree, raxml [see self.permitted_collection_types]
        what: the uncertain base type: N,M, or N_or_M
        has_outgroup: bool, whether there is an outgroup
        x: a collection of sample names, excluding any outgroup"""

        if collection_type not in self.permitted_collection_types:
            raise ValueError(
                "collection type {0} is not allowed".format(collection_type)
            )
        if what not in self.permitted_what:
            raise ValueError("what {0} is not allowed".format(what))
        if not isinstance(has_outgroup, bool):
            raise TypeError("has_outgroup must be boolean")

        return "{0}|{1}|{2}|{3}".format(
            collection_type,
            what,
            self.has_outgroup2label[has_outgroup],
            self._hashComponents(x),
        )
