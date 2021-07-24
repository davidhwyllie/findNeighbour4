#!/usr/bin/env python
""" performs fast single-linkage and mixture-aware clustering using SNV distances

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 
"""
import networkit as nk
import uuid
import random
import datetime
import scipy.stats
import statistics
import pandas as pd
import progressbar
import logging
from findn.identify_sequence_set import IdentifySequenceSet


class MockPersistence:
    """simulates the fnPersistence class which provides access to findNeighbour stored data;
    the objective is to allow testing of the Clustering class, which requires a Persistence object to be provided to allow it to access findNeighbour SNP distances.
    therefore, only methods relevant to SNP distance provision are provided.
    The MockPersistence class generates data resembling a collection of samples which are related, their SNP distances,
    These methods are:
    guids()    lists the names of sequences present
    guid2neighbours  links between guids - returns type 1 output

    it also supports an isMixed() method, which is used for simulating whether a sample is mixed or not.
    """

    def cluster_delete_legacy(self, name):
        """ delete any legacy data in the mock persistence store """
        pass
        return

    def __init__(self, n_guids: int):
        """ starts up a MockPersistence object; generates a set of pairwise links compatible with n_guids being part of a set of similar sequences. """
        self.latest_version = 0
        self.latest_version_behaviour = "increment"

        self.node2name = {}
        self.name2node = {}
        self.name2clusterid = {}
        self.node2clusterid = {}
        self._guid2neighbours = {}
        self.store = {}
        self.g = nk.generators.ClusteredRandomGraphGenerator(n_guids, int(n_guids / 2), 1, 0).generate()
        for x in self.g.iterNodes():
            new_guid = str(uuid.uuid4())
            self.node2name[x] = new_guid
            self.name2node[new_guid] = x
            self._guid2neighbours[new_guid] = []

        # determine connected components ('single linkage clusters')
        cc = nk.components.ConnectedComponents(self.g)
        cc.run()
        for clusterid, component in enumerate(cc.getComponents()):
            for node in component:
                self.name2clusterid[self.node2name[node]] = clusterid
                self.node2clusterid[node] = clusterid

        # build a dictionary containing edges (SNV distances) in format 1
        for (x, y) in self.g.iterEdges():
            guid1 = self.node2name[x]
            guid2 = self.node2name[y]
            snv = random.sample(range(7), 1)[0]  # distances drawn randomly from 0-6
            self._guid2neighbours[guid1].append([guid2, snv])
            self._guid2neighbours[guid2].append([guid1, snv])

    def cluster_latest_version(self, clustering_version):
        """ returns fake version information; increments if latest_version_behaviour is 'increment' """
        if self.latest_version_behaviour == "increment":
            self.latest_version += 1
        return self.latest_version

    def guids(self):
        """ returns all guids (sequence identifiers) in the network """
        return set(self.name2node.keys())

    def guids_valid(self):
        """ returns all guids (sequence identifiers) in the network , which are assumed to be valid"""
        return set(self.name2node.keys())

    def guid2neighbours(self, guid: str, returned_format: int = 1) -> dict:
        """returns neighbours of a guid in type 1 format [[guid1, snvdist1], [guid2, snvdist2], ...]

        note: attempts to change the returned_format result in a NotImplementedError

        Parameters:
        guid: the identifier of the sequence whose neighbours are sought
        returned_type: a placeholder, only a value of 1 is accepted
        """
        if not returned_format == 1:
            raise NotImplementedError(
                "the MockPersistence.guid2neighbours() method always returns type 1 output, but type {0} was requested".format(
                    returned_format
                )
            )
        return {"guid": guid, "neighbours": self._guid2neighbours[guid]}

    def cluster_store(self, key, serialisation):
        """ stores serialisation in a dictionary using key """
        self.store[key] = serialisation
        return None

    def cluster_read(self, key):
        try:
            return self.store[key]
        except KeyError:
            return None


class MixtureChecker:
    """ abstract base class for mixture checkers.  do not use """

    def __init__(self, **kwargs):
        """ does nothing """
        self.info = "This is the arbitrary base class; do not use "

    def is_mixed(self, guid):
        """ does not check for mixtures """
        return {"mix_check_method": self.info, "is_mixed": False}


class NoMixtureChecker(MixtureChecker):
    """ a class which implements a MixtureChecker which does not check for mixtures """

    def __init__(self):
        """ does nothing """
        self.info = "No check"


class MixtureCheckerTest(MixtureChecker):
    """ a class which implements a MixtureChecker which sets guids as mixed if they begin with a number """

    def __init__(self):
        """ does nothing """
        self.info = "Test check"

    def is_mixed(self, guid):
        """ does not check for mixtures """
        return {"mix_check_method": self.info, "is_mixed": guid[0] in ["0", "1", "2"]}


class MixPOREMixtureChecker(MixtureChecker):
    """ a class which implements a MixtureChecker based on the mixPORE approach """

    def __init__(
        self,
        HYBRIDCOMPARER,
        snv_threshold,
        mixture_criterion,
        cutoff,
        uncertain_base_type,
        max_seqs_mixpore,
        mixture_snv_cutoff,
        **kwargs,
    ):
        """sets up the MixtureChecker.

        Input:
                a hybridComparer object, configured to analyse stored sequences in a findNeighbour4 persistence store.
                snv_threshold - the clustering threshold; samples are joined if <= snv from each other.

                'exclude':  the clusters do not include guids with 'is_mixed'=True properties.
                            mixed samples exist only as single-element clusters.
                'include':  samples are included in clusters to which they are
                            similar.
                            One guid can belong to more than one cluster.

                uncertain_base_type: dictates which bases are considered in computations about mixtures;
                             valid values are 'N' 'M' 'N_or_M'.

                mixture_criterion: pvalue_n, where p = 1,2,3 or 4.  refers to the nature of the statistical test peformed; see  hybridComparer.msa() function for details.  pvalue_2 is recommended.
                cutoff:  the p value at which the result is considered mixed.  not corrected for multiple comparisons.  one comparison is performed per sample analysed. consider 1e-8
                snv_cutoff:  if >0, the result will not be considered mixed if align{what} < snv_cutoff.  For example, if less than 4 mixed bases in the alignment do not compromise clustering, snv_cutoff could be set to 4.
                max_seqs_mixpore:  how many sequences the test sequence will be compared with to perform mixpore.  2 minimum.  recommended: 10

                **kwargs: other args are ignored
                Note that in all the below documentation, 'guid' refers to a node identified by a guid,
                    and 'guids' to multiple such nodes.
        """
        self.info = "MixPORE"
        self.hc = HYBRIDCOMPARER
        self.snv_threshold = snv_threshold
        self.mixture_criterion = mixture_criterion
        self.p_value_cutoff = cutoff
        self.snv_cutoff = mixture_snv_cutoff
        self.uncertain_base_type = uncertain_base_type
        self.max_seqs_mixpore = max_seqs_mixpore
        self.sequence_ms = {}
        # we need to load enough samples into the hybridComparer to do reliable composition based computations.  The hybridcomparer doesn't contain all the sequences, and doesn't need to, but it does need to contain enough to compute composition data used by msa().  Important: the hybridComparer must have disable_insertion = True, to avoid possible conflicts with different catwalks.

        # if not self.hc.disable_insertion:
        #    raise NotImplementedError("Cannot use a hybridComparer with insertion enabled ")

        self.ensure_composition()

    def ensure_composition(self):
        """ loads enough enough samples into the hybridComparer to do reliable composition based computations.  The hybridcomparer attached to this object doesn't contain all the sequences, and doesn't need to, but it does need to contain enough to compute composition data used by msa(). """
        n_current_sequences = len(self.hc.pc.guids())

        if n_current_sequences < 100:  # the target
            self.hc.repopulate_sample(100)  # will do up to 100

    def is_mixed(self, guid):
        """ check for mixtures using mixpore method."""
        self.ensure_composition()
        neighbours = self.hc.PERSIST.guid2neighbours(guid, returned_format=1)["neighbours"]
        neighbours = sorted(neighbours, key=lambda x: int(x[1]))
        neighbours = [x for x in neighbours if x[1] <= self.snv_threshold]

        # if we have more neighbours than we are allowed to analyse (self.max_seqs_mixpore) then we preferentially analyse the ones with high mixed base numbers
        if len(neighbours) > self.max_seqs_mixpore:
            guid2mix = []
            for neighbour in neighbours:
                if not neighbour[0] in self.sequence_ms.keys():
                    nMixed = 0
                    try:
                        nMixed = self.hc.PERSIST.guid_annotation(neighbour[0])[neighbour[0]]["DNAQuality:mixed"]
                    except KeyError:
                        pass  # not recorded
                    self.sequence_ms[neighbour[0]] = nMixed
                guid2mix.append([neighbour[0], self.sequence_ms[neighbour[0]]])
            neighbours = sorted(guid2mix, key=lambda x: -int(x[1]))

        if len(neighbours) > 1:  # need 2 or more to apply mixpore
            neighbours_analysed = neighbours[: self.max_seqs_mixpore]  # take up to specific number
            to_msa = neighbours_analysed.copy()
            to_msa.append([guid, 0])

            guids_to_msa = [x[0] for x in to_msa]
            msa_result = self.hc.multi_sequence_alignment(guids_to_msa, uncertain_base_type=self.uncertain_base_type)
            res = msa_result.df.loc[guid].to_dict()
            del res["aligned_seq"]  # not useful, wastes ram
            del res["aligned_mseq"]  # not useful, wastes ram

            enough_bases_mixed = None
            res["is_mixed"] = False
            if res[self.mixture_criterion] is not None:
                align_mixed = "align{0}".format(self.uncertain_base_type)
                if align_mixed in res.keys():
                    if res[align_mixed] > self.snv_cutoff:
                        enough_bases_mixed = True
                    else:
                        enough_bases_mixed = False
                else:
                    enough_bases_mixed = False
                if res[self.mixture_criterion] < self.p_value_cutoff and enough_bases_mixed:
                    res.update({"is_mixed": True})
            res["enough_bases_mixed"] = enough_bases_mixed

            ## return result
            res.update({"mix_check_method": self.info})

            return res
        else:
            # not assessed
            return {"mix_check_method": self.info, "is_mixed": None}


class MixtureAwareLinkageResult:
    """provides access to stored clustering from MixtureAwareLinkage

    It exposes the following methods:
    guid2clustermeta(after_change_id)
    guid2clusters(guid)
    .is_mixed(guid)
    .refresh()

    It exposes the following properties:
    .parameters : a dictionary, which may include
        snv_threshold
        uncertain_base_type
        mixed_sample_management
        mixture_criterion
        cutoff
        cluster_nomenclature_method

    Parameters also include the following:
    .clustering_time
    .refresh_time
    .guid2cluster (whole dictionary)
    .cluster2guid
    .change_id

    It is a lightweight; all it does is load relevant data from the clustering collection in mongo, perform limited & rapid rearrangement (to allow fast indexed access) and exposes it.
    The class does not itself do any clustering; it uses data written into the database by findNeighbour4_clustering in order to provide clustering information to front end services.
    """

    def __init__(self, serialisation=None, PERSIST=None, name="Test"):
        """loads data from serialisation, which is a dictionary generate by the
        MixtureAwareLinkage.serialise_output() method.
        parameters:
            serialisation: a dictionary generated by the MixtureAwareLinkage.serialise_output() method.  May be None; if none, will attempt to load a serialisation from the PERSIST object (see below)
            PERSIST:  a findNeighbour persistence object, which interfaces with Mongo.  If both this an serialisation are None, raises an error.
            name: the name of the clustering process saving the data
        returns:
            None if there is not clustering data; or a changeId reflecting the clustering version.
        """
        self.loaded_version = None
        self.PERSIST = PERSIST
        self.name = name
        self.storage_key = "{0}-{1}".format(self.name, "output")

        self.refresh(serialisation=serialisation)

    def refresh(self, serialisation=None):
        """reloads data from persistence store, if necessary.
        > If serialisation is a dictionary, loads from that
        > If serialisation is None and self.PERSIST is not None:
            A data source exists; load a new version from it if it exists.
        > If serialisation is None and self.PERSIST is None:
            This is a first run situation.
        """
        if isinstance(serialisation, dict):
            self._deserialise_from_dict(serialisation)
            self.loaded_version_load_time = datetime.datetime.now()
            self.loaded_version = None  # not linked to a stored version

            return

        latest_version = self.PERSIST.cluster_latest_version(self.storage_key)
        if serialisation is None and self.PERSIST is not None:  # try to recover data
            if self.loaded_version is not None and self.loaded_version == latest_version:  # no update needed
                return

            elif self.loaded_version is None or (
                not self.loaded_version == latest_version
            ):  # nothing stored in ram,  or latest version not loaded;

                serialisation = self._recover_serialisation()
                if serialisation is not None:
                    self._deserialise_from_dict(serialisation)
                    self.loaded_version = latest_version

                    self.loaded_version_load_time = datetime.datetime.now()
                    return

        # otherwise
        # First run
        self._first_run()

    def _recover_serialisation(self):
        """reads a serialisation of  the clustering output from the Mongo db.
        cluster
        """

        if self.PERSIST is None:
            raise ValueError("Must define a PERSIST Object to persist to disc")
        return self.PERSIST.cluster_read(self.storage_key)

    def _deserialise_from_dict(self, serialisation):
        """ recovers the clustering data from a serialisation dictionary """
        expected_keys = set(["parameters", "clustering_time", "guid2clustermeta", "clusterid2clusterlabel", "name"])
        if not set(serialisation.keys()) == expected_keys:
            raise KeyError(
                "Expected keys are not present in serialisation output: {0} vs {1}".format(
                    set(serialisation.keys()), expected_keys
                )
            )

        self.parameters = serialisation["parameters"]
        self.refresh_time = datetime.datetime.now().isoformat()
        self.guid2cluster = serialisation["guid2clustermeta"]
        self._clusterid2clusterlabel = self._dictkey2int(serialisation["clusterid2clusterlabel"])
        self.uncertain_base_type = self.parameters["uncertain_base_type"]
        self.snv_threshold = self.parameters["snv_threshold"]
        # compute change_id

        try:
            self.change_id = max([x["add_change_id"] for x in self.guid2cluster.values()])
        except ValueError:  # no data
            self.change_id = -1

        # build cluster2guid lookup
        self.cluster2guid = {}
        for guid in self.guid2cluster.keys():
            cluster_ids = self.guid2cluster[guid]["cluster_id"]
            for cluster_id in cluster_ids:
                if cluster_id not in self.cluster2guid.keys():
                    self.cluster2guid[cluster_id] = []
                self.cluster2guid[cluster_id].append(guid)
        self.current_version_load_time = datetime.datetime.now()

    def _first_run(self):
        """ sets up an empty clustering entry """
        self.parameters = {"note": "No data found"}
        self.refresh_time = datetime.datetime.now().isoformat()
        self.guid2cluster = {}
        self._clusterid2clusterlabel = {}
        self.uncertain_base_type = "?"
        self.snv_threshold = None
        self.current_version_load_time = datetime.datetime.now()

        # compute change_id
        self.change_id = 0

        # build cluster2guid lookup
        self.cluster2guid = {}

    def guids(self):
        """ returns the clustered guids """
        return set(self.guid2cluster.keys())

    def clusters2guid(self):
        """ returns a cluster2guid lookup """
        return self.cluster2guid

    def guid2clustermeta(self):
        """" returns guid2 cluster lookup """
        return self.guid2cluster

    def guid2clusters(self, guid):
        """ returns which clusters the guid belongs to """
        try:
            cluster_ids = self.guid2cluster[guid]["cluster_id"]
        except KeyError:
            return None
        retVal = []
        for cluster_id in cluster_ids:
            try:
                clusterlabel = self._clusterid2clusterlabel[cluster_id]["cluster_label"]
            except KeyError:
                clusterlabel = "-"
            retVal.append(
                {
                    "clustering_algorithm": self.name,
                    "guid": guid,
                    "cluster_id": cluster_id,
                    "cluster_label": clusterlabel,
                    "cluster_loadtime": self.current_version_load_time.isoformat(),
                }
            )
        return retVal

    def is_mixed(self, guid, reportUnknownAsFalse=True):
        """returns mixture status.
        internally, findNeighbour4 scores this as three values: True, False, and None (=not assessable)
        if reportUnknownAsFalse is True, reports none (unknown) values as False (this is FindNeighbour3's behaviour)
        reports None if guids does not exist"""
        try:
            mix = self.guid2cluster[guid]["is_mixed"]
        except KeyError:
            return None
        if mix is None and reportUnknownAsFalse:
            return False
        else:
            return mix

    def _dictkey2string(self, outputdict):
        """ converts the keys of a dictionary from int to string """
        retVal = {}
        for key in outputdict:
            retVal[str(key)] = outputdict[key]
        return retVal

    def _dictkey2int(self, outputdict):
        """ converts the keys of a dictionary from  string to int """
        retVal = {}
        for key in outputdict:
            retVal[int(key)] = outputdict[key]
        return retVal

    def clusters2guidmeta(self, after_change_id=None):
        """ returns a cluster -> guid mapping """

        retVal = []
        for guid in sorted(self.guids()):

            # set changeid when cluster was updated or when mixture was detected, whichever is later; in practice, they are likley to be the same.
            change_id = self.guid2cluster[guid]["latest_change_id"]  # refers to cluste
            is_mixed = self.is_mixed(guid)
            if is_mixed:
                if self.guid2cluster[guid]["mix_detected_at_change_id"] > change_id:
                    change_id = self.guid2cluster[guid]["mix_detected_at_change_id"]

            for cluster_id in self.guid2cluster[guid]["cluster_id"]:
                try:
                    clusterlabel = self._clusterid2clusterlabel[cluster_id]["cluster_label"]
                except KeyError:
                    clusterlabel = "-"
                if (after_change_id is None) or (change_id > after_change_id):
                    retVal.append(
                        {
                            "guid": guid,
                            "cluster_id": cluster_id,
                            "cluster_label": clusterlabel,
                            "change_id": change_id,
                            "is_mixed": is_mixed,
                        }
                    )
        return retVal


class MixtureAwareLinkage:
    """joins samples (identified by guids) which are less than some SNP cutoff apart, but handles mixtures.
    The clusters generated are identified by integer cluster_ids; these cluster_ids are not guaranteed to be stable.
    This class does not generate stable identifiers (which we call cluster_labels); however, it allows these to be set using the set_cluster_labels
    method, and does persist them with the cluster in an atomic manner.

    Mixtures are not computed by this class; they are computed by an external class which is provided to this class to use.

    The implementation initially (optionally) ignores mixed samples and uses an in-memory graph.  An undirected graph is constructed where an edge represents a pairwise distance less than cutoff.  A connected component represents a cluster.  Optionally, mixed samples are subsequently added to these clusters.

    This implementation has very high capacity- it has been tested with up to 100 million samples and about 3 billion edges.
    However, to reduce the number of edges, the class can reduce edges if edge numbers become large.   The simplify() method does this.  This isn't called automatically at present.

    It preserves the connectedness of samples in a cluster, but reduces the number of edges markedly (for tb data, by about 95%).  After simplification, statistics such as betweenness, centrality and degree cannot be meaningfully applied to the resulting graph, but we don't need to do this.
    """

    def __init__(
        self,
        PERSIST,
        MIXCHECK=None,
        snv_threshold: int = 20,
        serialisation=None,
        mixed_sample_management="ignore",
        parameters={},
        name="NoName",
    ):
        """creates a MixtureAwareLinkage object

        parameters:
                snv_threshold:  sample pairs with SNP distances large than this will be ignored when adding
                serialisation: a python representation of a MixtureAwareLinkage object, as generated by self.serialise(), which is a way of serialising the object for storage as json.
                                If None, a fresh (empty) graph is returned.
                                If not None, snv_threshold is ignored, with the snv_threshold taken from from_dict.
                PERSIST: a findNeighbour persistence object.
                         The guids(), guids2neighbours(), cluster_read and cluster_store methods are used.  For unittesting, an instance of MockPersistence can be used.
                MIXCHECK: an object whose is_mixed() method, called with a guid, returns a dictionary including an is_mixed key, a boolean indicating whether the sample is mixed or not.  Mixed samples are not used in the primary clustering.  If not supplied or None, a NoMixtureChecker object is used, which does not check for mixtures.


                mixed_sample_management: dictates how mixed samples are dealt with
                'ignore': the clustering ignores the 'is_mixed' property.  This is the behaviour of standard 'snp address' and related approaches.

                'exclude':  the clusters do not include guids with 'is_mixed'=True properties.
                            mixed samples exist only as single-element clusters.
                'include':  samples are included in clusters to which they are
                            similar.
                            One guid can belong to more than one cluster.

                'parameters': a dictionary of other parameters, which are not used but are made available to the  MixtureAwareLinkageResults class.
                'name' :    the name of the clustering.  Must be compatible with being a linux file name.  example: SNV12_Ignore
        returns:
                None
        """

        self.PERSIST = PERSIST  # store the persistence in the MixtureAwareLinkage object
        self.MIXCHECK = MIXCHECK  # a mixture checker object

        if self.MIXCHECK is None:
            self.MIXCHECK = NoMixtureChecker()

        # check it is the right class
        if not isinstance(self.MIXCHECK, MixtureChecker):
            raise TypeError("MIXCHECK must be a Mixture Checker")
        self.iss = IdentifySequenceSet()

        self._node2name = {}
        self._name2node = {}
        self.clustered_at_change_id = None
        self._name2meta = {}
        self.cluster2names = {}
        self.name2cluster = {}
        self.is_simplified = False
        self.mixed_sample_management = mixed_sample_management
        self.parameters = parameters  #
        self.name = name
        self._clusterid2clusterlabel = {}  # a dictionary of the type {1:{'cluster_label':'AA0041'}}

        start_afresh = False
        if serialisation is None:
            if self.PERSIST is None:
                start_afresh = True
            else:
                serialisation = self._recover_serialisation()

                if serialisation is None:
                    start_afresh = True
            if start_afresh:
                # create a new graph
                self.snv_threshold = snv_threshold
                self.g = nk.graph.Graph(weighted=False, directed=False, n=0)  # empty graph
                self.cc = nk.components.ConnectedComponents(self.g)

            else:
                self._deserialise_from_dict(serialisation)

        else:
            self._deserialise_from_dict(serialisation)
        self.dc = nk.centrality.DegreeCentrality(self.g, normalized=False, ignoreSelfLoops=True)

    def name2meta(self):
        """returns the guid to metadata (including mixture and if appropriate clustering data) information as a pandas dataframe.
        there is one row per guid.
        clusters are presented as a list; if guid aa1234 is in cluster 1, clusters are [1]"""
        res = pd.DataFrame.from_dict(self._name2meta, orient="index")

        # add in clustering if available
        if len(self.name2cluster) > 0:
            # there is clustering data
            res = res.merge(
                pd.DataFrame.from_dict(self.name2cluster, orient="index"), left_index=True, right_index=True, how="left"
            )

        return res

    def guid2cluster_labels(self):
        """List the labelled clusters in which guids belong.  If no labels for the cluster have
        been assigned, none will be listed.  Clusters will only be included if at least two
        cluster members exist.

        Parameters:
        None

        Returns:
        For guids in clusters of size at least 2, returns a dictionary of the type
        {'guid1':['AA001','AA002'], 'guid2':['AA003','AA004']}"""
        retVal = {}
        for cl in self.cluster2names.keys():
            if cl in self._clusterid2clusterlabel.keys():
                members = self.cluster2names[cl]["guids"]
                if len(members) >= 2:
                    for guid in members:
                        try:
                            retVal[guid] = set()
                        except KeyError:
                            pass
                        retVal[guid].add(self._clusterid2clusterlabel[cl]["cluster_label"])
            for guid in retVal:
                retVal[guid] = list(retVal[guid])
        return retVal

    def existing_labels(self):
        """ returns a list of the existing cluster labels """
        existing_labels = set()
        for labels in self._clusterid2clusterlabel.values():
            existing_labels.add(labels["cluster_label"])
        return list(sorted(existing_labels))

    def guid2clustermeta(self):
        """Parameters:
                    None

                    Returns
        a guid to metadata lookup (including mixture and if appropriate clustering data) information as a pandas dataframe.
                    one guid can exist in more than one cluster.
                    if guid aa1234 is in clusters 1 and 2, there are two rows added"""
        res = pd.DataFrame.from_dict(self._name2meta, orient="index")

        # add in clustering if available
        if len(self.name2cluster) > 0:
            # there is clustering data
            res = res.merge(
                pd.DataFrame.from_dict(self.name2cluster, orient="index"), left_index=True, right_index=True, how="left"
            )

        return res

    def serialise(self):
        """serialises the graph to a dictionary.
        The following are serialised:
            self.snv_threshold # integer
            self.g  # serialised as a string in EdgeList format
            self.node2name # dictionary
            self.name
            The following are not deserialised and are regenerated on load:
                self._name2node = {}

                self.cc = nk.components.ConnectedComponents(self.g) # this is an object, and is not serialisable

                # the below, which are set by cluster()
                self.clustered_at_change_id = None
                self.cluster2names= {}
                self.name2cluster = {}
        """

        edges = []
        for (u, v) in self.g.iterEdges():
            edges.append([u, v])

        # write it back to file
        retVal = {
            "snv_threshold": self.snv_threshold,
            "_edges": edges,
            "name": self.name,
            "is_simplified": self.is_simplified,
            "_node2name": self._node2name,
            "_name2meta": self._name2meta,
            "mixed_sample_management": self.mixed_sample_management,
            "parameters": self.parameters,
            "clusterid2clusterlabel": self._dictkey2string(self._clusterid2clusterlabel),
        }

        return retVal

    def _dictkey2string(self, outputdict):
        """ converts the keys of a dictionary from int to string """
        retVal = {}
        for key in outputdict:
            retVal[str(key)] = outputdict[key]
        return retVal

    def _dictkey2int(self, outputdict):
        """ converts the keys of a dictionary from  string to int """
        retVal = {}
        for key in outputdict:
            retVal[int(key)] = outputdict[key]
        return retVal

    def persist(self, what):
        """stores a serialisation of either the clustering graph itself, or the output in the Mongo db.

        what is one of 'graph' or 'output'
        Will store the serialisation in the clusters collection, using the key
        {name}-{what}

        Note: will not remove old versions
        call .remove_legacy() to do this.

        """
        if self.PERSIST is None:
            raise ValueError("Must define a PERSIST Object to persist to disc")
        storage_key = "{0}-{1}".format(self.name, what)
        if what == "graph":
            serialisation = self.serialise()
        elif what == "output":
            serialisation = self.serialise_output()
        else:
            raise ValueError("what must be one of graph, output")

        self.PERSIST.cluster_store(storage_key, serialisation)

        return

    def remove_legacy(self):
        """removes any legacy versions of clustering.
        This is a good idea because otherwise every time a clustering is performed, a new version will be stored (which happens routinely) and old ones retained (which is not normally relevant): while an audit trail may be helpful during development, this in general undesirable as very large amounts of disc space (hundreds of gigabytes) can readily be consumed"""
        self.PERSIST.cluster_delete_legacy(self.name)

    def _recover_serialisation(self):
        """reads a serialisation of  the clustering graph from the Mongo db.


        Will recover the serialisation from the clusters collection, using the key
        {name}-graph

        """
        what = "graph"
        if self.PERSIST is None:
            raise ValueError("Must define a PERSIST Object to persist to disc")
        storage_key = "{0}-{1}".format(self.name, what)
        return self.PERSIST.cluster_read(storage_key)

    def to_dict(self):
        """ serialises the object to a dictionary.  synonym of serialise() """
        return self.serialise()

    def _deserialise_from_dict(self, sdict):
        """deserialises the graph from dictionary sdict.
        The following are recovered from the dictionary:
            self.snv_threshold # integer
            self.g  # regenerated from nodes and edges
            self._node2name # dictionary
            self._clusterid2clusterlabel,

            The following are regenerated on load:
                self._name2node = {}

                self.cc = nk.components.ConnectedComponents(self.g) # this is an object, and is not serialisable
               # the below, which are set by cluster()
                self.clustered_at_change_id = None
                self.cluster2names= {}
                self.name2cluster = {}
        """

        # test the dictionary passed has the expected keys
        if not set(sdict.keys()) == set(
            [
                "snv_threshold",
                "parameters",
                "name",
                "_name2meta",
                "mixed_sample_management",
                "is_simplified",
                "mixed_sample_management",
                "_edges",
                "_node2name",
                "_name2meta",
                "clusterid2clusterlabel",
            ]
        ):
            raise KeyError("Dictionary passed does not have the right keys: got {0}".format(sdict.keys()))

        # create a new graph from serialisation
        self.snv_threshold = sdict["snv_threshold"]
        self.is_simplified = sdict["is_simplified"]
        self._name2meta = sdict["_name2meta"]
        self.mixed_sample_management = sdict["mixed_sample_management"]
        self.parameters = sdict["parameters"]
        self.name = sdict["name"]
        self._clusterid2clusterlabel = self._dictkey2int(sdict["clusterid2clusterlabel"])

        # make new graph
        old_node2name = dict(
            zip([int(x) for x in sdict["_node2name"].keys()], sdict["_node2name"].values())
        )  # keys are integers
        node_order = sorted([int(x) for x in old_node2name.keys()])  # order added

        self.g = nk.graph.Graph(
            weighted=False, directed=False, n=len(node_order)
        )  # empty graph with the correct number of nodes

        # map new node numbers (if different) to old node numbers
        old2new = dict(zip(node_order, [x for x in self.g.iterNodes()]))

        # add edges
        for [u, v] in sdict["_edges"]:
            nu = old2new[u]
            nv = old2new[v]
            self.g.addEdge(nu, nv)

        # add new node names
        self._name2node = {}
        self._node2name = {}
        for u in [int(x) for x in old_node2name.keys()]:  # old
            nu = old2new[u]
            guid = old_node2name[u]
            self._name2node[guid] = nu
            self._node2name[nu] = guid
        self.cc = nk.components.ConnectedComponents(self.g)  # object for finding components

        self.cluster()
        return

    def guids(self) -> set:
        """ returns all guids (sequence identifiers) in the MixtureAwareLinkage object """
        return set(self._name2node.keys())

    def update(self) -> int:
        """ adds any valid guids which are in PERSIST but not in the MixtureAwareLinkage object to the graph"""
        valid_guids = self.PERSIST.guids_valid()
        guids_to_add = valid_guids - self.guids()
        return self.add(guids_to_add)

    def add_sample(self, guid) -> int:
        """adds a single guid, and links between them and existing guids, to the graph.

        Mixture checking is performed after addition of the samples and the edges of mixed samples are removed.
        This algorithm does not itself perform clustering; call .cluster() to do so.

        parameters:
        guid: a guid (sample identifier) which needs adding


        returns:
        change_id: an integer which increases as more samples are added.  useful for identifying samples added after a particular point"""

        return self.add(set([guid]))

    def is_mixed(self, guid) -> bool:
        """returns whether a sample is mixed or not.

        parameters:
        guid: a guid (sample identifier) to check

        returns:
        bool
        if unknown or untested, returns None"""

        try:
            return self._name2meta[guid]["is_mixed"]
        except KeyError:
            return None

    def add(self, guids_to_add: set) -> int:
        """adds guids_to_add, and links between them and existing guids, to the graph.

        Mixture checking is performed after addition of the samples and the edges of mixed samples are removed.
        This algorithm does not itself perform clustering; call .cluster() to do so.

        parameters:
        guids_to_add: a list of guids (sample identifiers) which need adding


        returns:
        change_id: an integer which increases as more samples are added.  returns None if nothing changed.
                    useful for identifying samples added after a particular point"""

        change_id = self.add_without_mixture_checking(guids_to_add)

        if change_id is None:
            # nothing added
            return None
        # now prepare to examine for mixtures, recording degree for each node.
        existing_guids = self.guids()
        guids_potentially_requiring_evaluation = set()
        # print("listing guids which need re-evaluating")
        bar = progressbar.ProgressBar(max_value=len(guids_to_add))
        for i, guid in enumerate(guids_to_add):
            bar.update(i + 1)
            neighbours = self.PERSIST.guid2neighbours(guid, returned_format=1)["neighbours"]

            neighbours = sorted(neighbours, key=lambda x: int(x[1]))
            neighbours = [
                x for x in neighbours if x[0] in existing_guids and x[1] <= self.snv_threshold
            ]  # only consider links to existing guids

            self._name2meta[guid]["nneighbours"] = len(neighbours)
            if len(neighbours) > 0:  # we don't assess samples with zero edges
                guids_potentially_requiring_evaluation.add(guid)
            for neighbour in neighbours:
                guids_potentially_requiring_evaluation.add(neighbour[0])
        bar.finish()

        # print("Counting neighbours")
        to_evaluate = guids_potentially_requiring_evaluation - guids_to_add
        bar = progressbar.ProgressBar(max_value=len(to_evaluate))

        # record the number of neighbours for each of these updated guids.  This is useful for analytics later, although it is not actually required for this process.
        # print("SNV threshold",self.snv_threshold)
        # print("Existing guids n=",len(existing_guids))
        # print("Guids potentially requiring evaluation n=",len(guids_potentially_requiring_evaluation))
        # print("To add guids n=",len(guids_to_add))
        # print("To evaluate, n=",len(to_evaluate))
        for i, guid in enumerate(to_evaluate):
            bar.update(i + 1)
            neighbours = self.PERSIST.guid2neighbours(guid, returned_format=1)["neighbours"]
            neighbours = sorted(neighbours, key=lambda x: int(x[1]))
            neighbours = [
                x for x in neighbours if x[0] in existing_guids and x[1] <= self.snv_threshold
            ]  # only consider links to existing guids
            self._name2meta[guid]["nneighbours"] = len(neighbours)

        bar.finish()

        # remove anything known to be mixed from what needs to be mixture checked.
        already_mixed = set()
        for guid in self._name2meta.keys():
            if self._name2meta[guid]["is_mixed"] is True:  # already mixed
                already_mixed.add(guid)

        guids_potentially_requiring_evaluation = guids_potentially_requiring_evaluation - already_mixed

        # assess whether these are mixed
        logging.info("Checking mixture status")
        bar = progressbar.ProgressBar(max_value=len(guids_potentially_requiring_evaluation))
        # we unlink everything - when new samples are added, they can re-link known mixed samples
        for i, guid in enumerate(guids_potentially_requiring_evaluation):
            bar.update(i + 1)
            mix_result = self.MIXCHECK.is_mixed(guid)
            self._name2meta[guid].update(mix_result)
            self._name2meta[guid]["mix_checked_at_change_id"] = change_id
            if self._name2meta[guid]["is_mixed"]:  # check whether it is mixed
                self._name2meta[guid]["mix_detected_at_change_id"] = change_id
                already_mixed.add(guid)
        bar.finish()
        # if we are told to ignore mixtures, then we are finished.
        if self.mixed_sample_management == "ignore":
            return change_id
        else:
            # now remove all the edges from the mixed samples
            self.remove_edges(already_mixed)
            return change_id

    def add_without_mixture_checking(self, guids_to_add: list) -> int:
        """adds guids_to_add, and links between them and existing guids, to the graph.
        Mixture checking is not performed
        parameters:
        guids_to_add: a list of guids (sample identifiers) which need adding

        returns:
        change_id: an integer which increases as more samples are added.  returns None if nothing added.
        useful for identifying samples added after a particular point"""
        change_id = self._add_nodes_and_links(guids_to_add, check_edges=False)
        return change_id

    def ensure_edges(self, guids_to_add: list) -> int:
        """ensures guids_to_add are nodes in the graph,
        and adds all links between them and existing guids, to the graph.

        parameters:
        guids_to_add: a list of guids (sample identifiers) which need adding

        returns:
        change_id: an integer which increases as more samples are added.  useful for identifying samples added after a particular point"""
        retVal = self._add_nodes_and_links(guids_to_add, check_edges=True)
        if set(guids_to_add) == self.guids():  # it's not longer simplified: all edges are present
            self.is_simplified = False
        return retVal

    def _add_nodes_and_links(self, guids_to_add: list, check_edges=False) -> int:
        """adds guids_to_add, and links between them and existing guids, to the graph
        Mixture checking is not performed.

        parameters:
        guids_to_add: a list of guids (sample identifiers) which need adding

        check_edges:  should normally be false.  if it is necessary to restore all edges of guids_to_add (for example, in a graph which has been simplified)
                        set check_edges=True.  All edges of guids_to_add will be restored.

        returns:
        change_id: an integer which increases as more samples are added.  useful for identifying samples added after a particular point  returns None if nothing added"""

        remaining_guids_to_add = set(guids_to_add) - self.guids()  # don't try to readd anything which is already in the graph
        valid_guids = guids_to_add.union(
            self.guids()
        )  # we only add links to these : the new ones, and anything already there.
        node_ids_added = set()

        # add any guids which are not present as nodes
        for guid_to_add in remaining_guids_to_add:
            node_id = self.g.addNode()
            self._node2name[node_id] = guid_to_add
            self._name2node[guid_to_add] = node_id
            self._name2meta[guid_to_add] = {}

        # add any edges necessary.
        # if check_edges is False, we just add the edges of the new nodes in remaining_guids_to_add
        # if check_edges is True, we  add any edges which are not in the graph for guids_to_add
        if check_edges:
            guids_whose_edges_to_add = guids_to_add
        else:
            guids_whose_edges_to_add = remaining_guids_to_add
        logging.info("build graph with all edges")
        bar = progressbar.ProgressBar(max_value=len(guids_whose_edges_to_add))

        for i, guid_to_add in enumerate(guids_whose_edges_to_add):
            bar.update(i + 1)

            node_id_1 = self._name2node[guid_to_add]
            node_ids_added.add(node_id_1)

            guid2neighbour_dict = self.PERSIST.guid2neighbours(guid_to_add, returned_format=1)
            guid2neighbours = guid2neighbour_dict["neighbours"]

            # only add links which are less or equal to than the cutoff value and which involve vertices (sequences) in the graph
            guid2neighbours = self._filter_guid2neighbours_by_snpcutoff(guid2neighbours, self.snv_threshold)
            guid2neighbours = self._filter_guid2neighbours_by_targets(guid2neighbours, valid_guids)
            for item in guid2neighbours:
                node_id_2 = self._name2node[item[0]]  # lookup the integer node id for the target of the edge
                if not (self.g.hasEdge(node_id_1, node_id_2) or self.g.hasEdge(node_id_2, node_id_1)):
                    self.g.addEdge(node_id_1, node_id_2)
        bar.finish()
        if len(node_ids_added) == 0:
            return -1  # nothing added, no changes

        change_id = max(
            node_ids_added
        )  # this is the changeid: an integer number which is guaranteed to increase as more samples are added
        for guid_to_add in remaining_guids_to_add:
            # record that is has not been mix checked
            self._name2meta[guid_to_add].update({"mix_check_method": "No check", "is_mixed": None, "add_change_id": change_id})
        return change_id

    def cluster(self) -> int:
        """Performs clustering.
        Parameters:
        None

        Returns:
        clustered_at_change_id: see below

        Side effects:
        Sets the following properties of the MixtureAwareLinkage object:
        - clustered_at_change_id:       the change_id up to which clustering has been performed.
        - cluster2names:                a dictionary comprising:  {cluster_id: {'latest_change_id':latest_change_id, 'guids':members, 'member_hash': member_hash} ,... where:
                                                cluster_id is an integer cluster_id
                                                member_hash is a hash representing cluster membership;
                                                latest_change_id is an integer reflecting when this cluster last changed;
                                                members are a list of guids in the cluster
        - name2cluster:                 a dictionary comprising {guid: {'latest_change_id':latest_change_id,'cluster_id'}}

        """

        # cluster by node_id.  This deals with the components with edges.
        # if mixed_sample_management = 'ignore', this is all the samples.
        # otherwise, it's only the unmixed samples.
        existing_guids = self.guids()

        # if there are no samples, we return
        if len(existing_guids) == 0:
            return 0

        clusters = self._connectedComponents(what="node_id")

        node_ids = set()

        # iterate over clusters defined from unmixed samples.
        for key in clusters.keys():
            cluster = clusters[key]
            cluster_id = min(cluster)  # smallest id

            change_id_this_cluster = max(cluster)  # latest time of change
            guids = [self._node2name[x] for x in cluster]  # guids
            node_ids.add(change_id_this_cluster)  # largest change_id in this cluster

            for guid in guids:
                self.name2cluster[guid] = {"latest_change_id": change_id_this_cluster, "cluster_id": [cluster_id]}

        # note the maximum node_id
        if len(node_ids) == 0:
            self.clustered_at_change_id = -1
        else:
            self.clustered_at_change_id = max(node_ids)

        # ignore if we are ignoring mixing samples, or if we are in 'exclude' mode (for which the clustering has already generated the right result - as the edges have all been removed from mixed cases), as in this mode they go into their own clusters as singletons.
        if self.mixed_sample_management == "include":

            # now we iterate over the mixed samples.
            mixed_samples = set()
            for guid in self._name2meta.keys():
                if self._name2meta[guid]["is_mixed"]:
                    mixed_samples.add(guid)

            for guid in mixed_samples:
                # load the neighbours of this sample.
                # find the clusters they belong to
                # add these to the cluster_id list in name2cluster, and cluster2names.

                # At the moment, each sample is in its own cluster
                initial_cluster_membership = set(self.name2cluster[guid]["cluster_id"])
                cluster_membership = set()
                neighbours = self.PERSIST.guid2neighbours(guid, returned_format=1)["neighbours"]
                neighbours = sorted(neighbours, key=lambda x: int(x[1]))

                neighbours = [
                    x
                    for x in neighbours
                    if x[0] in existing_guids and x[1] <= self.snv_threshold and not x[0] in mixed_samples
                ]  # only consider links to existing guids which are not known to be mixed
                for target, snpdist in neighbours:
                    # if we are including the samples, we put them in any matching clusters
                    cl_membership_this_target = self.name2cluster[target]["cluster_id"][0]
                    cluster_membership.add(cl_membership_this_target)

                # if we haven't assigned our mixed sample to any clusters, we assign it to its own
                if len(cluster_membership) == 0:
                    cluster_membership = initial_cluster_membership

                self.name2cluster[guid]["cluster_id"] = list(cluster_membership)

        # build self.cluster2names from self.name2cluster;
        self.cluster2names = {}
        for name in self.name2cluster.keys():
            cluster_ids = self.name2cluster[name]["cluster_id"]
            for cluster_id in cluster_ids:
                if cluster_id not in self.cluster2names.keys():
                    self.cluster2names[cluster_id] = {"guids": []}  # starting
                self.cluster2names[cluster_id]["guids"].append(name)

        # compute change_id and hash
        for cluster_id in self.cluster2names.keys():
            members = self.cluster2names[cluster_id]["guids"]
            hashed_guids = self.iss.make_identifier("cluster", "-", False, members)
            try:
                latest_change_id = max([self._name2node[x] for x in members])
            except ValueError:  # no data
                latest_change_id = -1

            self.cluster2names[cluster_id].update(
                {"latest_change_id": latest_change_id, "member_hash": hashed_guids}
            )  # starting

        return self.clustered_at_change_id

    def centrality(self, what: str = "degree") -> list:
        """returns degree centrality (number of links)
            Note that this reflects what is recorded in the graph; if the graph is simplified, it will not reflect the true number of links.
        parameters:
            what : the kind of centrality: 'degree' is implemented (only) at present
        returns:
            a dictionary, mapping node_id to centrality and Z-normalized centrality
        """
        centrality = self.dc.run().scores()
        # determine non-zero centralitiy
        nz_centrality = [x for x in centrality if x > 0]
        if len(nz_centrality) == 0:  # no data
            return {}
        mad = scipy.stats.median_abs_deviation(nz_centrality)
        if mad < 1:
            mad = 1  # lower bound on mad is 1
        med = statistics.median(nz_centrality)

        # return a dictionary
        retVal = {}
        for (node_id, node_centrality) in zip(self.g.iterNodes(), centrality):

            Z = (node_centrality - med) / mad
            retVal[self._node2name[node_id]] = {
                "node_id": node_id,
                "degree_centrality": node_centrality,
                "degree_Z": Z,
                "degree_valid": self.is_simplified is False,
            }  # warning if the graph is simplified
        return retVal

    def _connectedComponents(self, what: str = "name") -> list:
        """returns connected components (clusters).  A connected component has two or more members.

        parameters:
            what : either 'node_id' (internal node ids) or 'name' (guids)
        returns:
            a dictionary, keyed by a hash of the connected components, with a value of a list of their values
        """
        self.cc.run()
        retVal = {}
        for item in self.cc.getComponents():
            if what == "name":
                # then we return the guids
                returned_list = [self._node2name[x] for x in item]
            else:
                returned_list = item  # we return the node_ids
            sha1_item = self.iss.make_identifier("connections", "-", False, returned_list)

            retVal[sha1_item] = returned_list
        return retVal

    def _filter_guid2neighbours_by_snpcutoff(self, guid2neighbours: list, snpcutoff: int) -> list:
        """removes any elements from the list guid2neighbours if they have distance>snpcutoff

        parameters:
            guid2neighbours : a list of format [[guid1,distance1], [guid2,distance2]..] as returned by PERSIST.guid2neighbours()
            snpcutoff: an integer snp distance
        returns:
            a list of lists"""
        retVal = []

        for guid, snpdist in guid2neighbours:

            if snpdist <= snpcutoff:  # item[1] is the snp distance
                retVal.append([guid, snpdist])
        return retVal

    def _filter_guid2neighbours_by_targets(self, guid2neighbours: list, valid_targets: set) -> list:
        """removes any elements from the list guid2neighbours, which has format [[guid1,distance1], [guid2,distance2]..] as returned by PERSIST.guid2neighbours(), if guid is not in valid_targets

        parameters:
            guid2neighbours : a list of format [[guid1,distance1], [guid2,distance2]..] as returned by PERSIST.guid2neighbours()
            valid_targets: a set of guids (sample identifiers) for which edges are to be returned
        returns:
            a list of lists"""
        retVal = []

        for guid, snpdist in guid2neighbours:
            if guid in valid_targets:  # item[1] is the snp distance
                retVal.append([guid, snpdist])
        return retVal

    def remove_edges(self, to_remove) -> None:
        """remove edges from a list of guids.

        if the graph is simplified (i.e. not all edges are present) then all the edges of the connected components of guids will be restored.

        Parameters:
        to_remove: a string, or iterable contained guids whose edges should be removed.

        Returns:
        None

        """

        if isinstance(to_remove, str):
            to_remove = set([to_remove])  # if a single guid is specified, then we put that in an iterable.

        # if the graph is simplified, then all the edges of all the components are restored.
        if self.is_simplified:

            # identify the connected components of the nodes in to_remove
            to_ensure = set()
            ccs = self._connectedComponents(what="node_id")
            for cc_hash in ccs.keys():
                cc_members = sorted(ccs[cc_hash])
                if len(cc_members) > 1:
                    for element in cc_members:
                        to_ensure.add(self._node2name[element])
            # make sure all the edges of these nodes are present.
            self.ensure_edges(to_ensure)

        # remove the edges
        edges_to_remove = set()
        for guid in to_remove:
            node_id = self._name2node[guid]
            for neighbour in self.g.iterNeighbors(node_id):
                edges_to_remove.add(frozenset([node_id, neighbour]))

        for (u, v) in edges_to_remove:
            self.g.removeEdge(u, v)
        self.g.compactEdges()

        return None

    def simplify(self, after_change_id: int = None) -> None:
        """rewires the graph minimising edges by preserving the connected components

        Parameters:
        after_change_id: only update nodes inserted after_change_id, an integer provided by the .add() method to version additions.  If None, will simplify all parts of the graph.

        Returns:
        None"""

        # if no change_id is specified, we want to include all nodes.  Nodes start from zero.
        if after_change_id is None:
            after_change_id = -1
        ccs = self._connectedComponents(what="node_id")
        node2min_element_in_cluster = {}
        for cc_hash in ccs.keys():
            cc_members = sorted(ccs[cc_hash])
            min_element = cc_members[0]
            if len(cc_members) > 1:
                for element in cc_members:
                    node2min_element_in_cluster[element] = min_element

        dd = nk.centrality.DegreeCentrality(self.g).run().scores()  # compute number of edges for each node

        # iterate over all nodes; identify edges which need updating
        to_remove = set()
        to_add = set()
        for node_id, degree in enumerate(dd):
            if node_id > after_change_id:  # if it's in scope for update
                # test whether the node is part of a cluster
                if node_id in node2min_element_in_cluster:
                    # it is part of a cluster
                    if not node_id == node2min_element_in_cluster[node_id]:  # it's not the first element
                        # if the desired edge is the single edge present
                        if degree == 1 and self.g.hasEdge(node_id, node2min_element_in_cluster[node_id]):
                            pass
                        else:
                            to_add.add(frozenset([node2min_element_in_cluster[node_id], node_id]))
                            for neighbour in self.g.iterNeighbors(node_id):
                                to_remove.add(frozenset([node_id, neighbour]))

        for (u, v) in to_remove:
            self.g.removeEdge(u, v)
        for (u, v) in to_add:
            self.g.addEdge(u, v)
        self.g.compactEdges()
        self.is_simplified = True
        return None

    def serialise_output(self):
        """serialises the results of the clustering.  This can be read into a MixtureAwareClusteringResults object, as done by findNeighbour4.

        The following are serialised:
                guid2clustermeta
                clusterid2clusterlabel
                parameters : a dictionary, which may include
                        snv_threshold
                        uncertain_base_type
                        mixed_sample_management
                        mixture_criterion
                        cutoff

                clustering_time

        """

        return {
            "name": self.name,
            "parameters": self.parameters,
            "clustering_time": datetime.datetime.now().isoformat(),
            "clusterid2clusterlabel": self._dictkey2string(self._clusterid2clusterlabel),
            "guid2clustermeta": self.guid2clustermeta().to_dict(orient="index"),
        }

    def apply_cluster_labels(self, cl2label):
        """sets to lookup between the internal integer cluster_id and an externally supplied, stable cluster_label (such as AA041).
        Internally, the lookup is held in _clusterid2clusterlabel.

        Parameters:
            cl2label, a dictionary of the form {x}:{'cluster_label':'AA041'}
                where x is a cluster_id, and AA041 is the name of the cluster.

        The code runs the following checks:
        - the cluster_labels supplied are unique, and  map 1:1 to a cluster_id
        - the cluster_ids supplied all exist.
        - the values supplied are dictionaries with a 'cluster_label' key

        Application of the supplied cluster labels replaces previous entries.
        """

        # check: we have been supplied a dictionary of the right format.
        if not isinstance(cl2label, dict):
            raise ValueError("Must supply a dictionary, not a {0}".format(type(cl2label)))

        for cluster_id in cl2label.keys():
            label_dict = cl2label[cluster_id]
            if not isinstance(label_dict, dict):
                raise ValueError("Must supply a dictionary, not a {0} as the value of cl2label".format(type(cl2label)))
            if "cluster_label" not in label_dict.keys():
                raise KeyError("cluster_label must be a key; got {0}".format(label_dict))

        # checks: 1:1 mapping of cluster_id to cluster_labels;
        cluster_labels = set([x["cluster_label"] for x in cl2label.values()])
        cluster_ids = set(cl2label.keys())
        if not len(cluster_ids) == len(cluster_labels):
            # there is not a 1:1 mapping between cluster_id and cluster_label
            raise KeyError(
                "There is not a 1:1 mapping between cluster_ids and cluster_labels. In _labels but not in _ids:  {0}; in _ids but not in _labels {1}".format(
                    cluster_labels - cluster_ids, cluster_ids - cluster_labels
                )
            )

        # check: all clusters exist
        non_referenced_cluster_ids = cluster_ids - set(self.cluster2names.keys())
        if not non_referenced_cluster_ids == set([]):  # all clusterids exist
            raise KeyError("Not all clusterids exist; missing ones are :{0}".format(non_referenced_cluster_ids))

        # apply the labels
        # print("applying label",cl2label)
        self._clusterid2clusterlabel = cl2label

    def raise_error(self, token):
        """raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging"""
        raise ZeroDivisionError(token)
