#!/usr/bin/env python3
""" produces cluster numbers 

Clusters numbers are codes  of the type AB123
and are generated from an integer. 
System used by PHE/UKHSA to name TB clusters.

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 

"""


import pandas as pd
from collections import Counter


class ClusterNameAssigner:
    """ assigns cluster names to samples """

    def __init__(self, CLUSTERNOMENCLATURE):
        """creates a ClusterNameAssigner

        Parameters:
        CLUSTERNOMENCLATURE: a ClusterNomenclature object"""
        if not isinstance(CLUSTERNOMENCLATURE, ClusterNomenclature):
            raise TypeError("Must be passed a ClusterNomenclature object")
        self.CLUSTERNOMENCLATURE = CLUSTERNOMENCLATURE
        return

    def assign_new_clusternames(self, clusterid2guid, previous_guid2cluster_label={}):
        """assigns a cluster identifier, as provided by self.CLUSTERNOMENCLATURE, to each cluster.

        Parameters:
        cluster_id2guid: a dictionary mapping integer cluster_ids to their contents {cluster_id: {'latest_change_id':latest_change_id, 'guids':members, 'member_hash': member_hash} ,... where:
        previous_guid2cluster_label: a dictionary containing the *previous* cluster names assigned to sample names (usually guids).  Ignored if set to empty {}
        """

        # for each cluster, determine the majority previous call,
        recommended_cluster_label = {}
        for cluster_id in clusterid2guid:

            legacy_cluster_labels = []
            if len(clusterid2guid[cluster_id]["guids"]) > 1:  # we don't note clusterids for clusters of size 1

                # for each guid in the existing cluster
                for guid in clusterid2guid[cluster_id]["guids"]:
                    try:
                        # if there are legacy cluster label(s) associated with this guid
                        for legacy_cluster_label in previous_guid2cluster_label[guid]:  #
                            legacy_cluster_labels.append(legacy_cluster_label)  # note them

                    except KeyError:  # nothing provided
                        pass

                ctr = Counter(legacy_cluster_labels)  # determine the frequency of the legacy cluster labels in the cluster
                most_common = [
                    x[0] for x in ctr.most_common()
                ]  # if multiple cluster_ids have same frequency , there will be more than 1
                # the counter delivers these as tuples in descending order
                n_most_common = len(most_common)
                if n_most_common >= 1:
                    most_common = min(
                        most_common
                    )  # if there is more than one common example, we take the smallest/ earliest sample

                elif n_most_common == 0:  # none provided
                    most_common = self.CLUSTERNOMENCLATURE.new_label()

                recommended_cluster_label[cluster_id] = {"cluster_label": most_common}

        # now check whether the recommended cluster labels are unique to the clusters: in SQL:
        # select cluster_id, count(distinct recommended_cluster_label) group by cluster_id having count(distinct recommended_cluster_labels)>1

        df = pd.DataFrame.from_dict(recommended_cluster_label, orient="index")
        if len(df.index) == 0:
            return {}

        df["cluster_id"] = df.index.to_list()
        cnts = df["cluster_id"].groupby(df["cluster_label"]).agg(["count"])
        df = df.merge(cnts, left_on="cluster_label", right_on="cluster_label")

        # we sort this by cluster_label, count.  The intention is that the earliest cluster labels, which are the smallest alphabetically,
        # and the largest clusters, will be encountered first and will not change.  By contrast, smaller clusters may be renamed if necessary.
        df = df.sort_values(["cluster_label", "count"], ascending=[True, False])
        # print(df)
        already_used = set()
        for ix in df.index:
            # print(ix)
            if df.at[ix, "count"] > 1:
                # print("Count >1", df.at[ix,'cluster_label'])
                if df.at[ix, "cluster_label"] in already_used:  # we have used that label

                    new_label = self.CLUSTERNOMENCLATURE.new_label(from_existing=df.at[ix, "cluster_label"])
                    df.at[ix, "cluster_label"] = new_label
                    # print("Assigned",new_label)
                else:

                    # no change needed, but if we encounter this cluster label again, then we will need to reassign this label as it's not unique.
                    pass
            already_used.add(df.at[ix, "cluster_label"])
            recommended_cluster_label[df.at[ix, "cluster_id"]] = {"cluster_label": df.at[ix, "cluster_label"]}

        return recommended_cluster_label


class ClusterNomenclature:
    """provides names for clusters.
    It provides four functions:

    - keeps track of which labels are used
    - provides a new unused label4
    - provides a new unused label from existing labels
    - marks a label as in use.

    Exactly how this is is performed may vary between ClusterNomenclature objects, but
    the following methods are exposed:

    __init__ (existing_labels, deserialise_from)
    is_used(label)
    used_labels()
    serialise()
    new_label(from, is_used)
    use_label(label)
    """

    def __init__(self, cluster_nomenclature_method="integer", existing_labels=None, deserialise_from=None, **kwargs):
        """creates a ClusterNomenclature object

        Parameters;
        method: the method used to generate cluster names.
                'integer' - a new integer is assigned.
                Alternatives is
                'TB': identifiers of the type AA123 are generated, with AA123.1, AA123.2 etc generated if subdivisions are needed
                Nothing else is currently implemented.
        existing_labels: a list of existing labels, or None
        deserialise_from: a json string serialisation of the object, or None
        **kwargs: any other parameters

        Returns:
        None

        if existing_labels and deserialise_from are both None, a new empty object is created
        """

        self.moduli = [10, 10, 10, 24, 24]
        self.letters = [
            "A",
            "B",
            "C",
            "D",
            "E",
            "F",
            "G",
            "H",
            "J",
            "K",
            "L",
            "M",
            "N",
            "P",
            "Q",
            "R",
            "S",
            "T",
            "U",
            "V",
            "W",
            "X",
            "Y",
            "Z",
        ]
        self.numbers = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
        if deserialise_from is not None and existing_labels is not None:
            raise ValueError("You can either deserialise from a string, or supply existing_labels, but not both")
        implemented_cluster_nomenclature_methods = ["integer", "WBS", "TB"]
        if cluster_nomenclature_method not in implemented_cluster_nomenclature_methods:
            raise ValueError(
                "implemented ClusterNomenclature methods are {0}; {1} is not".format(
                    implemented_cluster_nomenclature_methods, cluster_nomenclature_method
                )
            )

        self.cluster_nomenclature_method = cluster_nomenclature_method
        if deserialise_from is not None:
            self._deserialise(deserialise_from)
        else:
            if existing_labels is None:
                self._labels = set()

            else:
                self._labels = set(existing_labels)

    def is_used(self, label):
        """tests whether a label is in use

        Parameters:
        label: a cluster label

        Returns:
        bool

        Raises:
        ValueError if label is None
        """
        if label is None:
            raise ValueError("label cannot be None")
        else:
            return label in self._labels

    def new_label(self, from_existing=None):
        """produce a new cluster label, optionally from an existing one

        Parameters:
        optionally, from_existing, which is an existing cluster name

        Returns:
        new cluster label
        """
        if from_existing is None:

            new_label = self._new_label()
        else:
            new_label = self._new_label_from_existing(from_existing=from_existing)
        self._labels.add(new_label)
        return new_label

    def _new_label_from_existing(self, from_existing):
        """produces a new label from an existing one
        Parameters:
        from_existing, an existing label

        Returns:
        a new cluster label, as specified by cluster_nomenclature_method.

        If cluster_nomenclature_method is 'integer', returns a new integer
        If cluster_nomenclature_method is 'TB', returns a cluster of the form AB123-3"""

        if self.cluster_nomenclature_method == "integer":
            return self._new_label()
        else:
            first_five_chars = from_existing[:5]  # compute a new label based on the existing one
            similar_labels = [
                x for x in self._labels if x[:5] == first_five_chars
            ]  # and find any existing labels beginning with these
            similar_labels = set(similar_labels)

            if len(similar_labels) == 0:  # should be impossible
                # error
                raise KeyError("internal error: {0} not found in existing labels".format(first_five_chars))
            else:
                versions = []
                for item in similar_labels:
                    if "v" in item:  # it's of the form AA123v6
                        versions.append(int(item.split("v")[1]))
                if len(versions) == 0:
                    new_label_id = 0
                else:
                    new_label_id = max(versions)
            return "{0}v{1}".format(first_five_chars, new_label_id + 1)

    def _new_label(self):
        """produce a new cluster label ab initio
        Parameters:
        None

        Returns:
        a new cluster label, as specified by cluster_nomenclature_method"""
        try:
            existing_max_label = max(self._labels)
        except ValueError:
            if self.cluster_nomenclature_method == "TB":
                existing_max_label = "AA000"
            else:
                existing_max_label = 0

        # if ClusterNomenclature is 'tb', then deliver a new TB type identifier
        if self.cluster_nomenclature_method == "TB":
            i = self._decode_TB(existing_max_label)
            new_label = self._encode_integer(i + 1)
        else:
            new_label = existing_max_label + 1

        return new_label

    def existing_labels(self):
        """returns a sorted list of used labels.

        Parameters:
        None

        Returns:
        sorted list of used labels"""

        return sorted(self._labels)

    def _encode_integer(self, x):
        """ encodes the integer x using two letters and three numbers, as performed with TB clusters """

        # check input
        if not isinstance(x, int):
            raise TypeError("x must be an integer")
        if not x >= 0 and x < 584000:
            raise ValueError("x must be between 0 and 584000; it is {0}".format(x))

        integer_representation = x
        digits = [None, None, None, None, None]

        for i in range(5):
            this_modulus = self.moduli[i]
            remainder = integer_representation % this_modulus
            integer_representation = integer_representation - remainder
            integer_representation = int(integer_representation / this_modulus)
            digits[i] = remainder

        digits = digits[::-1]  # reverse it
        cluster_code = "{0}{1}{2}{3}{4}".format(
            self.letters[digits[0]], self.letters[digits[1]], digits[2], digits[3], digits[4]
        )
        return cluster_code

    def _decode_TB(self, x):
        """ decodes the two letters and three numbers, as performed with TB clusters, into an integer """

        # check input

        if not isinstance(x, str):
            raise TypeError("x must be a string, not {0}.  x is {1}".format(type(x), x))

        if not len(x) >= 5:
            raise ValueError("must be of length >=5, but {0} is not".format(x))

        x = x[:5]  # chop anything > 5
        digits = list(x)
        if not (
            digits[0] in self.letters
            and digits[1] in self.letters
            and digits[2] in self.numbers
            and digits[3] in self.numbers
            and digits[4] in self.numbers
        ):
            raise ValueError("invalid format: need two letters and three numbers")
        digits_represent = [
            [i for i, n in enumerate(self.letters) if n == digits[0]][0] * 24000,
            [i for i, n in enumerate(self.letters) if n == digits[1]][0] * 1000,
            [i for i, n in enumerate(self.numbers) if n == digits[2]][0] * 100,
            [i for i, n in enumerate(self.numbers) if n == digits[3]][0] * 10,
            [i for i, n in enumerate(self.numbers) if n == digits[4]][0],
        ]
        return sum(digits_represent)

    def serialise(self):
        """  a json serialisable dictionary including the object's data """
        return {"cluster_nomenclature_method": self.cluster_nomenclature_method, "labels": self.existing_labels()}

    def _deserialise(self, deserialise_from):
        self._labels = set(deserialise_from["labels"])
        self.cluster_nomenclature_method = deserialise_from["cluster_nomenclature_method"]
