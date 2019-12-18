#!/usr/bin/env python
""" includes classes performing fast single-linkage and mixture-aware clustering """
import networkit as nk
import unittest
import uuid
import random
import os
import json
import hashlib
import datetime
import scipy.stats
import statistics
import warnings
import pandas as pd
import progressbar
import logging
from identify_sequence_set import IdentifySequenceSet

# for unittesting only
from hybridComparer import hybridComparer
from read_config import ReadConfig
from mongoStore import fn3persistence

class MockPersistence():
    """ simulates the fnPersistence class which provides access to findNeighbour stored data;
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

    def __init__(self, n_guids:int):
        """ starts up a MockPersistence object; generates a set of pairwise links compatible with n_guids being part of a set of similar sequences. """
        self.latest_version = 0
        self.latest_version_behaviour = 'increment'

        self.node2name = {}
        self.name2node = {}
        self.name2clusterid = {}
        self.node2clusterid = {}
        self._guid2neighbours={}
        self.store = {}
        self.g = nk.generators.ClusteredRandomGraphGenerator(n_guids, int(n_guids/2), 1, 0).generate()
        for x in self.g.nodes():
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
        for (x,y) in self.g.edges():
            guid1 = self.node2name[x]
            guid2 = self.node2name[y]
            snv = random.sample(range(7),1)[0]     # distances drawn randomly from 0-6
            self._guid2neighbours[guid1].append([guid2,snv])
            self._guid2neighbours[guid2].append([guid1,snv])
    def cluster_latest_version(self, clustering_version):
        """ returns fake version information; increments if latest_version_behaviour is 'increment' """
        if self.latest_version_behaviour=='increment':
            self.latest_version +=1
        return self.latest_version
        
    def guids(self):   
        """ returns all guids (sequence identifiers) in the network """
        return set(self.name2node.keys())

    def guid2neighbours(self, guid:str, returned_format:int=1)->dict:
        """ returns neighbours of a guid in type 1 format [[guid1, snvdist1], [guid2, snvdist2], ...]
        
        note: attempts to change the returned_format result in a NotImplementedError

        Parameters:
        guid: the identifier of the sequence whose neighbours are sought
        returned_type: a placeholder, only a value of 1 is accepted
        """
        if not returned_format==1:
            raise NotImplementedError("the MockPersistence.guid2neighbours() method always returns type 1 output, but type {0} was requested".format(returned_format))
        return {'guid':guid, 'neighbours':self._guid2neighbours[guid]}

    def cluster_store(self, key, serialisation):
        """ stores serialisation in a dictionary using key """
        self.store[key]=serialisation
        return None
    def cluster_read(self, key):
        try:
            return self.store[key]
        except KeyError:
            return None

class MixtureChecker():   
    """ abstract base class for mixture checkers.  do not use """
    def __init__(self, **kwargs):
        """ does nothing """
        self.info="This is the arbitrary base class; do not use "

    def is_mixed(self, guid):
        """ does not check for mixtures """
        return {'mix_check_method':self.info, 'is_mixed':False}
    
    
class NoMixtureChecker(MixtureChecker):   
    """ a class which implements a MixtureChecker which does not check for mixtures """
    def __init__(self):
        """ does nothing """
        self.info="No check"

class TestMixtureChecker(MixtureChecker):   
    """ a class which implements a MixtureChecker which sets guids as mixed if they begin with a number """
    def __init__(self):
        """ does nothing """
        self.info="Test check"
    def is_mixed(self, guid):
        """ does not check for mixtures """
        return {'mix_check_method':self.info, 'is_mixed':guid[0] in ['0','1','2']}

class MixPOREMixtureChecker(MixtureChecker):   
    """ a class which implements a MixtureChecker based on the mixPORE approach """
    def __init__(self, 
                    HYBRIDCOMPARER, 
                    snv_threshold,      
                    mixture_criterion, 
                    cutoff, 
                    uncertain_base_type,
                    max_seqs_mixpore,
                    **kwargs
            ):
        """ sets up the MixtureChecker.

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

                max_seqs_mixpore:  how many sequences the test sequence will be compared with to perform mixpore.  2 minimum.  recommended: 10

                **kwargs: other args are ignored
                Note that in all the below documentation, 'guid' refers to a node identified by a guid,
                    and 'guids' to multiple such nodes.

"""
        self.info="MixPORE"
        self.hc = HYBRIDCOMPARER
        self.snv_threshold= snv_threshold
        self.mixture_criterion = mixture_criterion
        self.cutoff = cutoff 
        self.uncertain_base_type = uncertain_base_type
        self.max_seqs_mixpore =max_seqs_mixpore

        # we need to load enough samples into the hybridComparer to do reliable composition based computations.  The hybridcomparer doesn't contain all the sequences, and doesn't need to, but it does need to contain enough to compute composition data used by msa().  Important: the hybridComparer must have disable_insertion = True, to avoid possible conflicts with different catwalks.

        #if not self.hc.disable_insertion:
        #    raise NotImplementedError("Cannot use a hybridComparer with insertion enabled ")

        self.ensure_composition()

    def ensure_composition(self):
        """ loads enough enough samples into the hybridComparer to do reliable composition based computations.  The hybridcomparer attached to this object doesn't contain all the sequences, and doesn't need to, but it does need to contain enough to compute composition data used by msa(). """
        n_current_sequences = len(self.hc.pc.guids())

        if n_current_sequences < 100:       # the target
            self.hc.repopulate_sample(100)  # will do up to 100

    def is_mixed(self, guid):
        """ check for mixtures using mixpore method."""
        self.ensure_composition()
        neighbours = self.hc.PERSIST.guid2neighbours(guid, returned_format=1)['neighbours']
        neighbours = sorted(neighbours, key = lambda x: int(x[1]))
        neighbours = [x for x in neighbours if x[1]<=self.snv_threshold]

        if len(neighbours)>1:		# need 2 or more to apply mixpore
            neighbours_analysed = neighbours[:self.max_seqs_mixpore]	# take up to specific number
            to_msa = neighbours_analysed.copy()
            to_msa.append([guid,0])			
            guids_to_msa = [x[0] for x in to_msa]
            msa_result = self.hc.multi_sequence_alignment(guids_to_msa,  uncertain_base_type=self.uncertain_base_type)
            res = msa_result.df.loc[guid].to_dict()
            del(res['aligned_seq'])		# not useful, wastes ram
            res['is_mixed'] = False
            if res[self.mixture_criterion] is not None:
                if res[self.mixture_criterion]<self.cutoff:
	                res['is_mixed']= True

            ## return result
            res.update({'mix_check_method':self.info})
            return res
        else:
            # not assessed
            return ({'mix_check_method':self.info, 'is_mixed':None})

class MixtureAwareLinkageResult():
    """ provides access to stored clustering from MixtureAwareLinkage

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
    The class does not itself do any clustering; it uses data written into the database by findNeighbour4-clustering in order to provide clustering information to front end services. 
    """

    def __init__(self, serialisation=None, PERSIST=None, name='Test'):
        """ loads data from serialisation, which is a dictionary generate by the
            MixtureAwareLinkage.serialise_output() method.
            parameters:
                serialisation: a dictionary generated by the MixtureAwareLinkage.serialise_output() method.  May be None; if none, will attempt to load a serialisation from the PERSIST object (see below)
                PERSIST:  a findNeighbour persistence object, which interfaces with Mongo.  If both this an serialisation are None, raises an error.
                name: the name of the clustering process saving the data
            returns:
                None if there is not clustering data; or a changeId reflecting the clustering version.
        """
        self.loaded_version = None
        self.PERSIST= PERSIST
        self.name = name
        self.storage_key = "{0}-{1}".format(self.name, 'output')
        
        self.refresh(serialisation=serialisation)

    def refresh(self, serialisation=None):
        """ reloads data from persistence store, if necessary.
        > If serialisation is a dictionary, loads from that
        > If serialisation is None and self.PERSIST is not None:
            A data source exists; load a new version from it if it exists.
        > If serialisation is None and self.PERSIST is None:
            This is a first run situation.
        """
        if isinstance(serialisation, dict):
            self._deserialise_from_dict(serialisation)
            self.loaded_version_load_time = datetime.datetime.now()
            self.loaded_version = None                                 # not linked to a stored version
            
            return

        latest_version = self.PERSIST.cluster_latest_version(self.storage_key) 
        if serialisation is None and self.PERSIST is not None:   # try to recover data
            if self.loaded_version is not None and self.loaded_version==latest_version:   # no update needed
                return
    
            elif self.loaded_version is None or (not self.loaded_version==latest_version):   # nothing stored in ram,  or latest version not loaded;

                serialisation = self._recover_serialisation()
                if serialisation is not None:
                    self._deserialise_from_dict(serialisation)
                    self.loaded_version= latest_version 
                    
                    self.loaded_version_load_time = datetime.datetime.now()
                    return

        # otherwise
        # First run 
        self._first_run()
  
    def _recover_serialisation(self):
        """ reads a serialisation of  the clustering output from the Mongo db.  

        """

        if self.PERSIST is None:
            raise ValueError("Must define a PERSIST Object to persist to disc")
        return self.PERSIST.cluster_read(self.storage_key)

    def _deserialise_from_dict(self, serialisation):
        """ recovers the clustering data from a serialisation dictionary """
        expected_keys = set(['parameters',
                             'clustering_time',
                             'guid2clustermeta',
                             'clusterid2clusterlabel',
                             'name'])
        if not set(serialisation.keys()) == expected_keys:
            raise KeyError("Expected keys are not present in serialisation output: {0} vs {1}".format(set(serialisation.keys()), expected_keys))

        self.parameters = serialisation['parameters']
        self.refresh_time = datetime.datetime.now().isoformat()
        self.guid2cluster = serialisation['guid2clustermeta']
        self._clusterid2clusterlabel = self._dictkey2int(serialisation['clusterid2clusterlabel'])
    
        self.uncertain_base_type = self.parameters['uncertain_base_type']
        self.snv_threshold = self.parameters['snv_threshold']
        # compute change_id

        try:
            self.change_id = max([x['add_change_id'] for x in self.guid2cluster.values()])
        except ValueError:      # no data
            self.change_id = -1
    
        # build cluster2guid lookup
        self.cluster2guid = {}
        for guid in self.guid2cluster.keys():
            cluster_ids = self.guid2cluster[guid]['cluster_id']
            for cluster_id in cluster_ids:
                if not cluster_id in self.cluster2guid.keys():
                    self.cluster2guid[cluster_id] = []
                self.cluster2guid[cluster_id].append(guid)
        self.current_version_load_time = datetime.datetime.now()
 
    def _first_run(self):
        """ sets up an empty clustering entry """
        self.parameters = {'note':'No data found'}
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
    def guid2clusters(self,guid):
        """ returns which clusters the guid belongs to """
        try:
            cluster_ids = self.guid2cluster[guid]['cluster_id']
        except KeyError:
            return None
        retVal=[]
        for cluster_id in cluster_ids:
            try:
                clusterlabel = self._clusterid2clusterlabel[cluster_id]['cluster_label']
            except KeyError:
                clusterlabel = '-'
            retVal.append({'cluster_id':cluster_id, 'cluster_label':clusterlabel})
        return retVal

    def is_mixed(self, guid, reportUnknownAsFalse=True):
        """ returns mixture status.
            internally, findNeighbour4 scores this as three values: True, False, and None (=not assessable)
            if reportUnknownAsFalse is True, reports none (unknown) values as False (this is FindNeighbour3's behaviour) 
            reports None if guids does not exist"""
        try:
            mix = self.guid2cluster[guid]['is_mixed']
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
            change_id = self.guid2cluster[guid]['latest_change_id']     # refers to cluste
            is_mixed = self.is_mixed(guid)            
            if is_mixed:
                if self.guid2cluster[guid]['mix_detected_at_change_id']> change_id:
                    change_id =  self.guid2cluster[guid]['mix_detected_at_change_id']
                
            for cluster_id in self.guid2cluster[guid]['cluster_id']:
                try:
                    clusterlabel = self._clusterid2clusterlabel[cluster_id]['cluster_label']
                except KeyError:
                    clusterlabel = '-'
                if (after_change_id is None) or (change_id > after_change_id):
                    retVal.append({'guid':guid, 'cluster_id':cluster_id, 'cluster_label':clusterlabel, 'change_id':change_id, 'is_mixed':is_mixed})
        return retVal
 
class MixtureAwareLinkage():
    """ joins samples (identified by guids) which are less than some SNP cutoff apart, but handles mixtures.
        The clusters generated are identified by integer cluster_ids; these cluster_ids are not guaranteed to be stable.
        This class does not generate stable identifiers (which we call cluster_labels); however, it allows these to be set using the set_cluster_labels
        method, and does persist them with the cluster in an atomic manner.

        Mixtures are not computed by this class; they are computed by an external class which is provided to this class to use.

        The implementation initially (optionally) ignores mixed samples and uses an in-memory graph.  An undirected graph is constructed where an edge represents a pairwise distance less than cutoff.  A connected component represents a cluster.  Optionally, mixed samples are subsequently added to these clusters.

        This implementation has very high capacity- it has been tested with up to 100 million samples and about 3 billion edges.  
        However, to reduce the number of edges, the class can reduce edges if edge numbers become large.   The simplify() method does this.  This isn't called automatically at present.

        It preserves the connectedness of samples in a cluster, but reduces the number of edges markedly (for tb data, by about 95%).  After simplification, statistics such as betweenness, centrality and degree cannot be meaningfully applied to the resulting graph, but we don't need to do this.
    """

    def __init__(self, 
                    PERSIST, 
                    MIXCHECK=None, 
                    snv_threshold:int=20, 
                    serialisation=None, 
                    mixed_sample_management = 'ignore',
                    parameters={},
                    name="NoName"):
        """ creates a MixtureAwareLinkage object

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

        self.PERSIST = PERSIST      # store the persistence in the MixtureAwareLinkage object
        self.MIXCHECK = MIXCHECK    # a mixture checker object
        
        if self.MIXCHECK is None:
            self.MIXCHECK = NoMixtureChecker()

        # check it is the right class
        if not isinstance(self.MIXCHECK, MixtureChecker):
            raise TypeError("MIXCHECK must be a Mixture Checker")
        self.iss = IdentifySequenceSet()

        self._node2name = {}
        self._name2node = {}
        self.clustered_at_change_id = None
        self._name2meta={}
        self.cluster2names= {}
        self.name2cluster = {}
        self.is_simplified = False
        self.mixed_sample_management = mixed_sample_management
        self.parameters = parameters        # 
        self.name=name
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
                self.g = nk.graph.Graph(weighted=False, directed=False, n=0)    # empty graph
                self.cc = nk.components.ConnectedComponents(self.g)

            else: 
                self._deserialise_from_dict(serialisation)

        else:
            self._deserialise_from_dict(serialisation)
        self.dc = nk.centrality.DegreeCentrality(self.g, normalized=False, ignoreSelfLoops=True)
    def name2meta(self):
        """ returns the guid to metadata (including mixture and if appropriate clustering data) information as a pandas dataframe.
            there is one row per guid.
            clusters are presented as a list; if guid aa1234 is in cluster 1, clusters are [1] """
        res = pd.DataFrame.from_dict(self._name2meta, orient='index')

        # add in clustering if available
        if len(self.name2cluster) > 0:
            # there is clustering data
            res = res.merge(pd.DataFrame.from_dict(self.name2cluster, orient='index'), left_index=True, right_index=True, how='left')
        
        return res

    def guid2cluster_labels(self):
        """ List the labelled clusters in which guids belong.  If no labels for the cluster have
            been assigned, none will be listed.  Clusters will only be included if at least two 
            cluster members exist.

            Parameters:
            None

            Returns:
            For guids in clusters of size at least 2, returns a dictionary of the type
            {'guid1':['AA001','AA002'], 'guid2':['AA003','AA004']} """
        retVal = {}
        for cl in self.cluster2names.keys():
            if cl in self._clusterid2clusterlabel.keys():
                members = self.cluster2names[cl]['guids']
                if len(members)>=2:
                    for guid in members:
                        try:
                            retVal[guid]=set()
                        except KeyError:
                            pass
                        retVal[guid].add(self._clusterid2clusterlabel[cl]['cluster_label'])
            for guid in retVal:
                retVal[guid] = list(retVal[guid])
        return retVal
                    
    def existing_labels(self):
        """ returns a list of the existing cluster labels """
        existing_labels = set()
        for labels in self._clusterid2clusterlabel.values():
            existing_labels.add(labels['cluster_label'])
        return list(sorted(existing_labels))

    def guid2clustermeta(self):
        """ Parameters:
            None 

            Returns
a guid to metadata lookup (including mixture and if appropriate clustering data) information as a pandas dataframe.
            one guid can exist in more than one cluster.
            if guid aa1234 is in clusters 1 and 2, there are two rows added"""
        res = pd.DataFrame.from_dict(self._name2meta, orient='index')

        # add in clustering if available
        if len(self.name2cluster) > 0:
            # there is clustering data
            res = res.merge(pd.DataFrame.from_dict(self.name2cluster, orient='index'), left_index=True, right_index=True, how='left')
        
        return res

    def serialise(self):
        """ serialises the graph to a dictionary.
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
        for (u,v) in self.g.edges():
            edges.append([u,v])

        # write it back to file
        retVal = {'snv_threshold':self.snv_threshold,
                  '_edges':edges,
                  'name':self.name,
                  'is_simplified':self.is_simplified, 
                  '_node2name':self._node2name,
                  '_name2meta':self._name2meta,
                  'mixed_sample_management':self.mixed_sample_management,
                  'parameters':self.parameters,
                  'clusterid2clusterlabel':self._dictkey2string(self._clusterid2clusterlabel)}

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
        """ stores a serialisation of either the clustering graph itself, or the output in the Mongo db.  

            what is one of 'graph' or 'output'
            Will store the serialisation in the clusters collection, using the key
            {name}-{what}

            Note: will not remove old versions
            call .remove_legacy() to do this.

        """
        if self.PERSIST is None:
            raise ValueError("Must define a PERSIST Object to persist to disc")
        storage_key = "{0}-{1}".format(self.name, what)
        if what == 'graph':
            serialisation = self.serialise()
        elif what == 'output':
            serialisation = self.serialise_output()
        else:
            raise ValueError("what must be one of graph, output")
        self.PERSIST.cluster_store(storage_key, serialisation)

        return

    def remove_legacy(self):
        """ removes any legacy versions of clustering.  
            This is a good idea because otherwise every time a clustering is performed, a new version will be stored (which happens routinely) and old ones retained (which is not normally relevant): while an audit trail may be helpful during development, this in general undesirable as very large amounts of disc space (hundreds of gigabytes) can readily be consumed """
        self.PERSIST.cluster_delete_legacy(self.name)

    def _recover_serialisation(self):
        """ reads a serialisation of  the clustering graph from the Mongo db.  


            Will recover the serialisation from the clusters collection, using the key
            {name}-graph

        """
        what='graph'
        if self.PERSIST is None:
            raise ValueError("Must define a PERSIST Object to persist to disc")
        storage_key = "{0}-{1}".format(self.name, what)
        return self.PERSIST.cluster_read(storage_key)

    def to_dict(self):
        """ serialises the object to a dictionary.  synonym of serialise() """
        return self.serialise()

    def _deserialise_from_dict(self, sdict):
        """ deserialises the graph from dictionary sdict.
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
        if not set(sdict.keys()) == set(['snv_threshold','parameters','name','_name2meta','mixed_sample_management','is_simplified','mixed_sample_management','_edges','_node2name', '_name2meta','clusterid2clusterlabel']):
            raise KeyError("Dictionary passed does not have the right keys: got {0}".format(sdict.keys()))

        # create a new graph from serialisation
        self.snv_threshold = sdict['snv_threshold']
        self.is_simplified = sdict['is_simplified']
        self._name2meta = sdict['_name2meta']
        self.mixed_sample_management = sdict['mixed_sample_management']
        self.parameters= sdict['parameters']
        self.name= sdict['name']
        self._clusterid2clusterlabel= self._dictkey2int(sdict['clusterid2clusterlabel'])
 
        # make new graph
        old_node2name = dict(zip([int(x) for x in sdict['_node2name'].keys()],sdict['_node2name'].values()))        # keys are integers
        node_order = sorted([int(x) for x in  old_node2name.keys()])        # order added

        self.g = nk.graph.Graph(weighted=False, directed=False, n=len(node_order))    # empty graph with the correct number of nodes

        # map new node numbers (if different) to old node numbers
        old2new = dict(zip(node_order, [x for x in self.g.nodes()]))

        # add edges
        for [u,v] in sdict['_edges']:
            nu = old2new[u]
            nv = old2new[v]
            self.g.addEdge(nu,nv)

        # add new node names
        self._name2node = {}
        self._node2name = {}
        for u in [int(x) for x in  old_node2name.keys()]:       # old
            nu = old2new[u]
            guid = old_node2name[u]
            self._name2node[guid]=nu
            self._node2name[nu]=guid
        self.cc = nk.components.ConnectedComponents(self.g)  # object for finding components

        self.cluster()
        return

    def guids(self)->set:   
        """ returns all guids (sequence identifiers) in the MixtureAwareLinkage object """
        return set(self._name2node.keys())

    def update(self)->int:
        """ adds any guids which are in PERSIST but not in the MixtureAwareLinkage object to the graph"""
        all_guids = self.PERSIST.guids()
        guids_to_add = all_guids-self.guids()
        return self.add(guids_to_add)

    def add_sample(self, guid)->int:
        """ adds a single guid, and links between them and existing guids, to the graph.
            
            Mixture checking is performed after addition of the samples and the edges of mixed samples are removed.
            This algorithm does not itself perform clustering; call .cluster() to do so.

            parameters:
            guid: a guid (sample identifier) which needs adding

            
            returns:
            change_id: an integer which increases as more samples are added.  useful for identifying samples added after a particular point"""

        return self.add(set([guid]))

    def is_mixed(self,guid)->bool:
        """ returns whether a sample is mixed or not.  

            parameters:
            guid: a guid (sample identifier) to check

            returns:
            bool
            if unknown or untested, returns None """

        try:
            return self._name2meta[guid]['is_mixed']
        except KeyError:
            return None

    def add(self, guids_to_add:set)->int:
        """ adds guids_to_add, and links between them and existing guids, to the graph.
            
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
        #print("listing guids which need re-evaluating") 
        bar = progressbar.ProgressBar(max_value=len(guids_to_add))
        for i,guid in enumerate(guids_to_add):
            bar.update(i+1)
            neighbours = self.PERSIST.guid2neighbours(guid, returned_format=1)['neighbours']

            neighbours = sorted(neighbours, key = lambda x: int(x[1]))
            neighbours = [x for x in neighbours if x[0] in existing_guids and x[1]<=self.snv_threshold]        # only consider links to existing guids
            
            self._name2meta[guid]['nneighbours'] = len(neighbours)
            if len(neighbours)>0:		# we don't assess samples with zero edges
    	        guids_potentially_requiring_evaluation.add(guid)
            for neighbour in neighbours:  
    	        guids_potentially_requiring_evaluation.add(neighbour[0])
        bar.finish()

        #print("Counting neighbours")
        to_evaluate = guids_potentially_requiring_evaluation-guids_to_add
        bar = progressbar.ProgressBar(max_value=len(to_evaluate))
        
        # record the number of neighbours for each of these updated guids.  This is useful for analytics later, although it is not actually required for this process.
        #print("SNV threshold",self.snv_threshold)
        #print("Existing guids n=",len(existing_guids))
        #print("Guids potentially requiring evaluation n=",len(guids_potentially_requiring_evaluation))
        #print("To add guids n=",len(guids_to_add))
        #print("To evaluate, n=",len(to_evaluate))
        for i,guid in enumerate(to_evaluate):
            bar.update(i+1)
            neighbours = self.PERSIST.guid2neighbours(guid, returned_format=1)['neighbours']
            neighbours = sorted(neighbours, key = lambda x: int(x[1]))
            neighbours = [x for x in neighbours if x[0] in existing_guids and x[1]<=self.snv_threshold]        # only consider links to existing guids
            self._name2meta[guid]['nneighbours'] = len(neighbours)
            
        bar.finish()

       # remove anything known to be mixed from what needs to be mixture checked.
        already_mixed = set()
        for guid in self._name2meta.keys():
	        if self._name2meta[guid]['is_mixed'] is True:		# already mixed
		        already_mixed.add(guid)

        guids_potentially_requiring_evaluation = guids_potentially_requiring_evaluation - already_mixed
        
        # assess whether these are mixed
        logging.info("Checking mixture status")
        bar = progressbar.ProgressBar(max_value=len(guids_potentially_requiring_evaluation))
        # we unlink everything - when new samples are added, they can re-link known mixed samples
        for i,guid in enumerate(guids_potentially_requiring_evaluation):
            bar.update(i+1) 
            mix_result = self.MIXCHECK.is_mixed(guid) 
            self._name2meta[guid].update(mix_result)
            self._name2meta[guid]['mix_checked_at_change_id']= change_id
            if self._name2meta[guid]['is_mixed']:        # check whether it is mixed
                self._name2meta[guid]['mix_detected_at_change_id']= change_id
                already_mixed.add(guid)
        bar.finish()
        # if we are told to ignore mixtures, then we are finished.
        if self.mixed_sample_management == 'ignore':
           return change_id
        else:
            # now remove all the edges from the mixed samples
            self.remove_edges(already_mixed)
            return change_id

    def add_without_mixture_checking(self, guids_to_add:list)->int:
        """ adds guids_to_add, and links between them and existing guids, to the graph.
            Mixture checking is not performed
            parameters:
            guids_to_add: a list of guids (sample identifiers) which need adding
         
            returns:
            change_id: an integer which increases as more samples are added.  returns None if nothing added.
            useful for identifying samples added after a particular point"""
        change_id = self._add_nodes_and_links(guids_to_add, check_edges =False)
        return change_id
    def ensure_edges(self, guids_to_add:list)->int:
        """ ensures guids_to_add are nodes in the graph,
            and adds all links between them and existing guids, to the graph.

            parameters:
            guids_to_add: a list of guids (sample identifiers) which need adding
            
            returns:
            change_id: an integer which increases as more samples are added.  useful for identifying samples added after a particular point"""
        retVal = self._add_nodes_and_links(guids_to_add, check_edges =True)
        if set(guids_to_add) == self.guids():         # it's not longer simplified: all edges are present
            self.is_simplified = False
        return retVal

    def _add_nodes_and_links(self, guids_to_add:list, check_edges=False)->int:
        """ adds guids_to_add, and links between them and existing guids, to the graph 
            Mixture checking is not performed.

            parameters:
            guids_to_add: a list of guids (sample identifiers) which need adding

            check_edges:  should normally be false.  if it is necessary to restore all edges of guids_to_add (for example, in a graph which has been simplified)
                            set check_edges=True.  All edges of guids_to_add will be restored.

            returns:
            change_id: an integer which increases as more samples are added.  useful for identifying samples added after a particular point  returns None if nothing added"""
        

        remaining_guids_to_add = set(guids_to_add)-self.guids()     # don't try to readd anything which is already in the graph
        valid_guids= guids_to_add.union(self.guids())               # we only add links to these : the new ones, and anything already there.
        node_ids_added = set()

        # add any guids which are not present as nodes
        nLoaded =0
        for guid_to_add in remaining_guids_to_add:
            node_id = self.g.addNode()
            self._node2name[node_id]=guid_to_add
            self._name2node[guid_to_add] =node_id
            self._name2meta[guid_to_add]={}
        
        # add any edges necessary.
        # if check_edges is False, we just add the edges of the new nodes in remaining_guids_to_add
        # if check_edges is True, we  add any edges which are not in the graph for guids_to_add
        if check_edges:
            guids_whose_edges_to_add = guids_to_add
        else:
            guids_whose_edges_to_add = remaining_guids_to_add
        logging.info("build graph with all edges")
        bar = progressbar.ProgressBar(max_value = len(guids_whose_edges_to_add))
   
        for i,guid_to_add in enumerate(guids_whose_edges_to_add):
            bar.update(i+1)
     
            node_id_1 = self._name2node[guid_to_add]
            node_ids_added.add(node_id_1)

            guid2neighbour_dict = self.PERSIST.guid2neighbours(guid_to_add, returned_format =1)
            guid2neighbours = guid2neighbour_dict['neighbours']

            # only add links which are less or equal to than the cutoff value and which involve vertices (sequences) in the graph
            guid2neighbours = self._filter_guid2neighbours_by_snpcutoff(guid2neighbours, self.snv_threshold)
            guid2neighbours = self._filter_guid2neighbours_by_targets(guid2neighbours, valid_guids)
            for item in guid2neighbours:
                node_id_2 = self._name2node[item[0]]        # lookup the integer node id for the target of the edge
                if not (self.g.hasEdge(node_id_1, node_id_2) or self.g.hasEdge(node_id_2, node_id_1)):
                    self.g.addEdge(node_id_1, node_id_2)
        bar.finish()
        if len(node_ids_added ) ==0:
            return -1                      # nothing added, no changes

        change_id = max(node_ids_added)      # this is the changeid: an integer number which is guaranteed to increase as more samples are added
        for guid_to_add in remaining_guids_to_add:
            # record that is has not been mix checked
            self._name2meta[guid_to_add].update({'mix_check_method':'No check', 'is_mixed':None, 'add_change_id':change_id})
        return change_id


    def cluster(self)->int:
        """ Performs clustering. 
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
        clusters = self._connectedComponents(what = 'node_id')

        node_ids = set()
        
        # iterate over clusters defined from unmixed samples.
        for key in clusters.keys():
            cluster = clusters[key]
            cluster_id = min(cluster)                       # smallest id
            
            change_id_this_cluster = max(cluster)           # latest time of change
            guids = [self._node2name[x] for x in cluster]   # guids 
            node_ids.add(change_id_this_cluster)            # largest change_id in this cluster

            for guid in guids:
                self.name2cluster[guid] = {'latest_change_id':change_id_this_cluster, 'cluster_id':[cluster_id]}

        # note the maximum node_id
        if len(node_ids) == 0:
            self.clustered_at_change_id = -1
        else:        
            self.clustered_at_change_id = max(node_ids)

        # ignore if we are ignoring mixing samples, or if we are in 'exclude' mode (for which the clustering has already generated the right result - as the edges have all been removed from mixed cases), as in this mode they go into their own clusters as singletons.
        if self.mixed_sample_management == 'include':
        
            # now we iterate over the mixed samples.
            mixed_samples = set()
            for guid in self._name2meta.keys():
                if self._name2meta[guid]['is_mixed']:
                    mixed_samples.add(guid)

            for guid in mixed_samples:              
                    # load the neighbours of this sample.
                    # find the clusters they belong to
                    # add these to the cluster_id list in name2cluster, and cluster2names.

                    # At the moment, each sample is in its own cluster
                    initial_cluster_membership = set(self.name2cluster[guid]['cluster_id'])
                    cluster_membership = set()
                    neighbours = self.PERSIST.guid2neighbours(guid, returned_format=1)['neighbours']
                    neighbours = sorted(neighbours, key = lambda x: int(x[1]))
                    
                    neighbours = [x for x in neighbours if x[0] in existing_guids and x[1]<=self.snv_threshold and not x[0]in mixed_samples]        # only consider links to existing guids which are not known to be mixed
                    for target, snpdist in neighbours:
                            # if we are including the samples, we put them in any matching clusters
                            cl_membership_this_target = self.name2cluster[target]['cluster_id'][0]
                            cluster_membership.add(cl_membership_this_target)

                    # if we haven't assigned our mixed sample to any clusters, we assign it to its own
                    if len(cluster_membership) == 0:
                        cluster_membership = initial_cluster_membership

                    self.name2cluster[guid]['cluster_id'] = list(cluster_membership) 

        # build self.cluster2names from self.name2cluster;
        self.cluster2names={}
        for name in self.name2cluster.keys():
            cluster_ids = self.name2cluster[name]['cluster_id']
            for cluster_id in cluster_ids:
                if not cluster_id in self.cluster2names.keys():
                    self.cluster2names[cluster_id] = {'guids':[]}        # starting
                self.cluster2names[cluster_id]['guids'].append(name)

        # compute change_id and hash
        for cluster_id in self.cluster2names.keys():
            members = self.cluster2names[cluster_id]['guids']
            hashed_guids = self.iss.make_identifier('cluster','-',False,members)
            try:
                latest_change_id = max([self._name2node[x] for x in members])
            except ValueError:      # no data
                latest_change_id = -1

            self.cluster2names[cluster_id].update({'latest_change_id':latest_change_id, 'member_hash': hashed_guids})        # starting

        return self.clustered_at_change_id

    def centrality(self, what:str='degree')->list:
        """ returns degree centrality (number of links)
            Note that this reflects what is recorded in the graph; if the graph is simplified, it will not reflect the true number of links.
        parameters:
            what : the kind of centrality: 'degree' is implemented (only) at present
        returns:
            a dictionary, mapping node_id to centrality and Z-normalized centrality
        """
        centrality = self.dc.run().scores()
        # determine non-zero centralitiy
        nz_centrality = [x for x in centrality if x>0]
        if len(nz_centrality)==0:       # no data
            return {}
        mad = scipy.stats.median_absolute_deviation(nz_centrality)
        if mad < 1:
            mad = 1     # lower bound on mad is 1 
        med = statistics.median(nz_centrality)
        
        # return a dictionary
        retVal = {}
        for (node_id, node_centrality) in zip(self.g.nodes(), centrality):

            Z= (node_centrality-med)/mad
            retVal[self._node2name[node_id]] = {'node_id':node_id, 
                                                'degree_centrality':node_centrality, 
                                                'degree_Z':Z,
                                                'degree_valid': self.is_simplified == False}            # warning if the graph is simplified
        return retVal

    def _connectedComponents(self, what:str='name')->list:
        """ returns connected components (clusters).  A connected component has two or more members.

        parameters:
            what : either 'node_id' (internal node ids) or 'name' (guids)
        returns:
            a dictionary, keyed by a hash of the connected components, with a value of a list of their values
        """
        self.cc.run()
        retVal = {}
        for item in self.cc.getComponents():
            min_node_id = min(item)

            if what == 'name':
                # then we return the guids
                returned_list = [self._node2name[x] for x in item]
            else:
                returned_list = item    # we return the node_ids
            sha1_item = self.iss.make_identifier('connections','-',False, returned_list)
    
            retVal[sha1_item] = returned_list
        return retVal   
    def _filter_guid2neighbours_by_snpcutoff(self, guid2neighbours:list, snpcutoff:int)->list:
        """ removes any elements from the list guid2neighbours if they have distance>snpcutoff

        parameters:
            guid2neighbours : a list of format [[guid1,distance1], [guid2,distance2]..] as returned by PERSIST.guid2neighbours()
            snpcutoff: an integer snp distance
        returns:
            a list of lists """
        retVal = []
 
        for guid, snpdist in guid2neighbours:

           if snpdist<=snpcutoff:      # item[1] is the snp distance
                retVal.append([guid,snpdist])
        return retVal

    def _filter_guid2neighbours_by_targets(self, guid2neighbours:list, valid_targets:set)->list:
        """ removes any elements from the list guid2neighbours, which has format [[guid1,distance1], [guid2,distance2]..] as returned by PERSIST.guid2neighbours(), if guid is not in valid_targets 

        parameters:
            guid2neighbours : a list of format [[guid1,distance1], [guid2,distance2]..] as returned by PERSIST.guid2neighbours()
            valid_targets: a set of guids (sample identifiers) for which edges are to be returned
        returns:
            a list of lists """
        retVal = []

        for guid, snpdist in guid2neighbours:
            if guid in valid_targets:      # item[1] is the snp distance
                retVal.append([guid, snpdist])
        return retVal

    def remove_edges(self, to_remove)->None:
        """ remove edges from a list of guids.
    
        if the graph is simplified (i.e. not all edges are present) then all the edges of the connected components of guids will be restored.

        Parameters:  
        to_remove: a string, or iterable contained guids whose edges should be removed. 

        Returns:
        None

        """

        if isinstance(to_remove, str):
            to_remove = set([to_remove])        # if a single guid is specified, then we put that in an iterable.

        # if the graph is simplified, then all the edges of all the components are restored.
        if self.is_simplified:

            #identify the connected components of the nodes in to_remove
            to_ensure = set()
            ccs = self._connectedComponents(what='node_id')
            for cc_hash in ccs.keys():
                cc_members = sorted(ccs[cc_hash])
                if len(cc_members)>1:
                    for element in cc_members:
                        to_ensure.add(self._node2name[element])
            # make sure all the edges of these nodes are present.
            self.ensure_edges(to_ensure)

        # remove the edges
        edges_to_remove = set()
        for guid in to_remove:
            node_id = self._name2node[guid]
            for neighbour in self.g.neighbors(node_id):
                edges_to_remove.add(frozenset([node_id,neighbour]))
       
        for (u,v) in edges_to_remove:
            self.g.removeEdge(u,v)
        self.g.compactEdges()
            
        return None
    def simplify(self, after_change_id:int=None)->None:
        """ rewires the graph minimising edges by preserving the connected components 
    
        Parameters:  
        after_change_id: only update nodes inserted after_change_id, an integer provided by the .add() method to version additions.  If None, will simplify all parts of the graph.

        Returns:
        None"""


        # if no change_id is specified, we want to include all nodes.  Nodes start from zero.
        if after_change_id is None:
            after_change_id = -1
        ccs = self._connectedComponents(what='node_id')
        node2min_element_in_cluster = {}
        for cc_hash in ccs.keys():
            cc_members = sorted(ccs[cc_hash])
            min_element = cc_members[0]
            if len(cc_members)>1:
                for element in cc_members:
                    node2min_element_in_cluster[element]=min_element

        dd= nk.centrality.DegreeCentrality(self.g).run().scores()       # compute number of edges for each node
        
        # iterate over all nodes; identify edges which need updating
        to_remove = set()
        to_add = set()
        for node_id,degree in enumerate(dd):
            if node_id > after_change_id:     # if it's in scope for update
                # test whether the node is part of a cluster
                if node_id in node2min_element_in_cluster:
                    # it is part of a cluster
                    if not node_id == node2min_element_in_cluster[node_id]:     # it's not the first element
                        rewire_edges = True
                        # if the desired edge is the single edge present
                        if degree == 1 and self.g.hasEdge(node_id,node2min_element_in_cluster[node_id]):
                            rewire_edges=  False
                        else:
                            rewire_edges = True
                            to_add.add(frozenset([node2min_element_in_cluster[node_id],node_id]))
                            for neighbour in self.g.neighbors(node_id):
                                to_remove.add(frozenset([node_id,neighbour]))
       
        for (u,v) in to_remove:
            self.g.removeEdge(u,v)
        for (u,v) in to_add:
            self.g.addEdge(u,v)
        self.g.compactEdges()
        self.is_simplified =True   
        return None

    def serialise_output(self):
        """ serialises the results of the clustering.  This can be read into a MixtureAwareClusteringResults object, as done by findNeighbour4.

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
            'name':self.name, 
            'parameters':self.parameters, 
            'clustering_time':datetime.datetime.now().isoformat(),
            'clusterid2clusterlabel':self._dictkey2string(self._clusterid2clusterlabel), 
            'guid2clustermeta':self.guid2clustermeta().to_dict(orient='index')
        }
    
    def apply_cluster_labels(self, cl2label):
        """ sets to lookup between the internal integer cluster_id and an externally supplied, stable cluster_label (such as AA041).
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
            if not 'cluster_label' in label_dict.keys():
                raise KeyError("cluster_label must be a key; got {0}".format(label_dict))

        # checks: 1:1 mapping of cluster_id to cluster_labels;
        cluster_labels = set([x['cluster_label'] for x in cl2label.values()])
        cluster_ids = set(cl2label.keys())
        if not len(cluster_ids) == len(cluster_labels):
            # there is not a 1:1 mapping between cluster_id and cluster_label
            raise KeyError("There is not a 1:1 mapping between cluster_ids and cluster_labels. In _labels but not in _ids:  {0} ; in _ids but not in _labels {1)".format(cluster_labels-cluster_ids, cluster_ids-cluster_labels))

        # check: all clusters exist
        non_referenced_cluster_ids = cluster_ids - set(self.cluster2names.keys())
        if not non_referenced_cluster_ids == set([]):      # all clusterids exist 
            raise KeyError("Not all clusterids exist; missing ones are :{0}".format(non_referenced_cluster_ids))          

        # apply the labels
        self._clusterid2clusterlabel = cl2label

     
    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging """
        raise ZeroDivisionError(token)

class Test_MP(unittest.TestCase):
    """ tests the MockPersistance object """
    def runTest(self):

        p = MockPersistence(n_guids =20)
        self.assertEqual(20, len(p.node2name))
        self.assertEqual(20, len(p.name2node))      

        # get guids
        guids = p.guids()
        self.assertEqual(len(guids),20)
        self.assertTrue(isinstance(guids, set))
        
        guid = min(guids)
        with self.assertRaises(NotImplementedError):
            res  = p.guid2neighbours(guid, returned_format=2)

        res  = p.guid2neighbours(guid, returned_format=1)
        self.assertTrue(isinstance(res,dict))

        to_store = {'one':1, 'two':2}   
        res = p.cluster_store('myKey',to_store)
        self.assertIsNone(res)
        res = p.cluster_read('NoKey')
        self.assertIsNone(res)
        res = p.cluster_read('myKey')
        
        self.assertEqual(to_store, res)
class Test_MAL_1(unittest.TestCase):
    """ tests the MixtureLinkage startup"""
    def runTest(self):
        p = MockPersistence(n_guids =20)
        m = MixtureAwareLinkage(PERSIST=p)

        # test there are no guids in the object on instatiation
        guids = m.guids()
        self.assertEqual(len(guids),0)
        self.assertTrue(isinstance(guids, set))
 
        with self.assertRaises(TypeError):       
            m = MixtureAwareLinkage(PERSIST=p, MIXCHECK=p)      # MIXCHECK is not a mixture checker


class Test_MAL_2(unittest.TestCase):
    """ tests the MixtureLinkage _filter_guid2neighbours_by_cutoff method"""
    def runTest(self):
        p = MockPersistence(n_guids =20)
        m = MixtureAwareLinkage(PERSIST=p)

        # test there are no guids in the object on instatiation
        input = [['a',6],['b',5],['c',4]]
        output = m._filter_guid2neighbours_by_snpcutoff(input, 5)
        self.assertEqual(output, [['b',5],['c',4]])

class Test_MAL_3(unittest.TestCase):
    """ tests the MixtureLinkage _filter_guid2neighbours_by_targets method"""
    def runTest(self):
        p = MockPersistence(n_guids =20)
        m = MixtureAwareLinkage(PERSIST=p)

        # test there are no guids in the object on instatiation
        input = [['a',6],['b',5],['c',4]]
        output = m._filter_guid2neighbours_by_targets(input, ['c'])
        self.assertEqual(output, [['c',4]])

class Test_MAL_4(unittest.TestCase):
    """ tests the MixtureLinkage add method"""
    def runTest(self):
        p = MockPersistence(n_guids =20)
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        change_id = m.add(guids_to_add)     
        self.assertEqual(m.guids(),p.guids())       # check all were added
        self.assertEqual(m.g.numberOfEdges(), p.g.numberOfEdges())       # same number of edges
        self.assertEqual(m.g.numberOfNodes(), p.g.numberOfNodes())       # same number of nodes
        self.assertEqual(m.g.numberOfNodes()-1, change_id)       # change_id is the nodeid of the last item inserted, which is zero indexed
        self.assertEqual(len(m._name2meta), p.g.numberOfNodes())         # one meta data entry per guid

class Test_MAL_4b(unittest.TestCase):
    """ tests the MixtureLinkage add_sample method"""
    def runTest(self):
        p = MockPersistence(n_guids =20)
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        for guid in guids_to_add:
            change_id = m.add_sample(guid)     
        self.assertEqual(m.guids(),p.guids())       # check all were added
        self.assertEqual(m.g.numberOfEdges(), p.g.numberOfEdges())       # same number of edges
        self.assertEqual(m.g.numberOfNodes(), p.g.numberOfNodes())       # same number of nodes
        self.assertEqual(m.g.numberOfNodes()-1, change_id)       # change_id is the nodeid of the last item inserted, which is zero indexed
        self.assertEqual(len(m._name2meta), p.g.numberOfNodes())         # one meta data entry per guid

class Test_MAL_5(unittest.TestCase):
    """ tests the _connectedComponents() method"""
    def runTest(self):
        p = MockPersistence(n_guids =20)
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        m.add(guids_to_add)     
        
        result = m._connectedComponents(what = 'node_id')
        for key in result.keys():       # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], int))

        result = m._connectedComponents(what = 'name')
        for key in result.keys():       # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], str))

class Test_MAL_7(unittest.TestCase):
    """ tests the _connectedComponents() method, at scale, with check of whether it got the right answer"""
    def runTest(self):

        n_guids = 20
        p = MockPersistence(n_guids = n_guids)
 
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)

        # add guids
        guids_to_add = p.guids()
        m.add(guids_to_add)     
  
        # we compute the number of clusters in the original (as defined for node by p.node2clusterid[node])

        # for each cluster in the test (defined by the keys of result)
        resultkey2original = dict()

        result = m._connectedComponents(what = 'name')
        for key in result.keys():       # keys are hashes of contents
            resultkey2original[key]=set()
        
        for key in result.keys():       # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], str))        # contents are lists of integers
            
            for i,name in enumerate(result[key]):                   # for each of our clusters
                resultkey2original[key].add(p.name2clusterid[name]) # we expect exactly one original cluster: i.e. the clustering works
        
        for key in resultkey2original.keys():
            self.assertEqual(len(resultkey2original[key]),1)
       

        # test guid conversion works
        result = m._connectedComponents(what = 'node_id')
        for key in result.keys():       # keys are hashes of contents
            self.assertTrue(isinstance(result[key][0], int))

class Test_MAL_8(unittest.TestCase):
    """ tests the update() method, which adds any guids not in the MixtureAwareLinkage object to it."""
    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())
        first_five = set(guids_to_add[0:5])
        m.add(first_five)
    
        self.assertEqual(len(m.guids()), 5)     
        self.assertEqual(m.guids(), first_five)     
      
        m.update()
        self.assertEqual(len(m.guids()), n_guids)     
      

        # check update adds all guids if none present
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        self.assertEqual(len(m.guids()), 0)     
        m.update()
        self.assertEqual(len(m.guids()), n_guids)     

class Test_MAL_9(unittest.TestCase):
    """ tests the simplify() and ensure_edges() methods, which respectively
            rewires a graph to minimise edges while preserving connected components.
            restore all edges"""
    def runTest(self):

        n_guids = 2000
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())    
        m.update()
        self.assertEqual(len(m.guids()), n_guids)     
      
        # ensure the connected components remain the same after simplication
        pre = m._connectedComponents(what = 'node_id')
        edges_pre = m.g.numberOfEdges()

        self.assertFalse(m.is_simplified)
        m.simplify()  
        self.assertTrue(m.is_simplified)
 
        post = m._connectedComponents(what = 'node_id')
        edges_post = m.g.numberOfEdges()
                
        self.assertEqual(pre, post)
        self.assertTrue(edges_pre > edges_post)

        m.ensure_edges(set(guids_to_add[0:10]))     # just a small number
        self.assertTrue(m.is_simplified)    # still simplified
        
        m.ensure_edges(set(guids_to_add))
        edges_restored = m.g.numberOfEdges()
        self.assertTrue(edges_restored > edges_post)
        self.assertEqual(edges_restored, edges_pre)
        self.assertFalse(m.is_simplified)
 
class Test_MAL_10(unittest.TestCase):
    """ tests the cluster() method, which persists transformation of the data generated by _connectedComponents()"""
    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())    
        m.update()
        self.assertEqual(len(m.guids()), n_guids)     


        self.assertTrue(m.clustered_at_change_id is None)
        self.assertEqual(0, len(m.name2cluster.keys()))
        self.assertEqual(0, len(m.cluster2names.keys()))
      
        # ensure the connected components remain the same after simplication
        m.cluster()   

        self.assertEqual(n_guids-1, m.clustered_at_change_id)
        self.assertEqual(n_guids, len(m.name2cluster.keys()))
        self.assertEqual(len(set([x for x in p.node2clusterid.values()])), len(m.cluster2names.keys()))

class Test_MAL_11(unittest.TestCase):
    """ tests the serialise() method, which returns a representation of the MixtureAwareLinkage object as a json string"""
    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())    
        m.update()
        
        retVal = m.serialise()
        self.assertEqual(set(retVal.keys()), set(['clusterid2clusterlabel','parameters','is_simplified','name','snv_threshold','mixed_sample_management','_edges','_node2name','_name2meta']))

        retVal2 = m.to_dict()       # synonym
        self.assertEqual(retVal, retVal2)

        # test jsonification round trip
        json_repr = json.dumps(retVal)
        retVal2 = json.loads(json_repr)

        # check we recover the graph
        m._deserialise_from_dict(retVal2)

        self.assertEqual(m.g.numberOfNodes(), p.g.numberOfNodes())
        self.assertEqual(m.g.numberOfEdges(), p.g.numberOfEdges())

        # check edges are all inplace
        for u,v in m.g.edges():
            gu = m._node2name[u]
            gv = m._node2name[v]
        
            ou = p.name2node[gu]
            ov = p.name2node[gv]
            self.assertTrue(p.g.hasEdge(ou,ov))     # check we wound up with what we started with


class Test_MAL_12(unittest.TestCase):
    """ tests recreation of MixtureAwareLinkage from serialisation """
    def runTest(self):

        n_guids = 50
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        #print("Creating",datetime.datetime.now())    
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())    
        m.update()
        #print("Done",datetime.datetime.now())
         
        retVal = m.serialise()
        self.assertEqual(set(retVal.keys()), set(['clusterid2clusterlabel','snv_threshold','name','parameters','is_simplified','_edges','mixed_sample_management','_node2name','_name2meta']))

        # test jsonification round trip
        json_repr = json.dumps(retVal)

        retVal2 = json.loads(json_repr)

        # check we recover the graph
        #print("Recovering",datetime.datetime.now())
        m2 = MixtureAwareLinkage(PERSIST = p, serialisation = retVal2)
        #print("Done",datetime.datetime.now())
 
        self.assertEqual(m2.g.numberOfNodes(), p.g.numberOfNodes())
        self.assertEqual(m2.g.numberOfEdges(), p.g.numberOfEdges())

        # check edges are all inplace
        for u,v in m2.g.edges():
            gu = m2._node2name[u]
            gv = m2._node2name[v]
        
            ou = p.name2node[gu]
            ov = p.name2node[gv]
            self.assertTrue(p.g.hasEdge(ou,ov))     # check we wound up with what we started with

class Test_MAL_13(unittest.TestCase):
    """ tests centrality measurement """
    def runTest(self):

        n_guids = 50
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        #print("Creating",datetime.datetime.now())    
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())    
        m.update()
       
        retVal = m.centrality()
        self.assertEqual(len(retVal.keys()), len(p.guids()))     # one for each guid

class Test_MAL_14a(unittest.TestCase):
    """ tests the remove_edges() method, which 
            removes the edges of one or more guids"""
    def runTest(self):

        n_guids = 2000
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20, mixed_sample_management = 'include')  # with this setting, mixtures are evaluated and edge removals occur
        guids_to_add = list(p.guids())    
        m.update()
        self.assertEqual(len(m.guids()), n_guids)     
      
        # the graph is not simplified.
        self.assertFalse(m.is_simplified)
 
        # remove the edges of  200 guids
        edges_pre = m.g.numberOfEdges() 
        guids_to_remove = guids_to_add[0:200]

        # check it works
        pre = m._connectedComponents(what = 'node_id')
        m.remove_edges(set(guids_to_remove))
        edges_post = m.g.numberOfEdges()
        self.assertTrue(edges_post < edges_pre)
     
        # simplify and check again
        self.assertFalse(m.is_simplified)
        m.simplify()  
        self.assertTrue(m.is_simplified)
 
        # remove the edges of  200 guids
        edges_pre = m.g.numberOfEdges() 
        guids_to_remove = guids_to_add[200:400]

        # check it works post simplification
        pre = m._connectedComponents(what = 'node_id')
        m.remove_edges(set(guids_to_remove))
        edges_post = m.g.numberOfEdges()


class Test_MAL_14b(unittest.TestCase):
    """ tests the add() method, which 
            removes the edges of one or more guids if they are mixed and mixed_sample_management is set"""
    def runTest(self):

        n_guids = 2000
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m1 = MixtureAwareLinkage(PERSIST=p, MIXCHECK=TestMixtureChecker(), snv_threshold=20, mixed_sample_management = 'include')  # with this setting, mixtures are  not evaluated and edge removals don't occur
        guids_to_add = list(p.guids())    
        m1.update()
        self.assertEqual(len(m1.guids()), n_guids)     
      
        # the graph is not simplified.
        self.assertFalse(m1.is_simplified)
 
        # check update adds remaining guids
        m2 = MixtureAwareLinkage(PERSIST=p, MIXCHECK=TestMixtureChecker(), snv_threshold=20, mixed_sample_management = 'ignore')  # with this setting, mixtures are  not evaluated and edge removals don't occur
        guids_to_add = list(p.guids())    
        m2.update()
        self.assertEqual(len(m2.guids()), n_guids)     
      
        # the graph is not simplified.
        self.assertFalse(m2.is_simplified)
 
        self.assertTrue(m1.g.numberOfEdges() < m2.g.numberOfEdges())
class Test_MAL_15(unittest.TestCase):
    """ tests the name2meta method, which exports metadata about the samples"""
    def runTest(self):

        n_guids = 10
        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, snv_threshold=20)
        guids_to_add = list(p.guids())    
        m.update()
        self.assertEqual(len(m.guids()), n_guids)     
        self.assertEqual(len(m.name2meta()), n_guids)
    
        # check output after clustering
        m.cluster()   

        self.assertEqual(len(m.name2meta()), n_guids)

class Test_TMT_1(unittest.TestCase):
    """ tests the setting of mixed samples by the test mixture checker"""
    def runTest(self):

        MIXCHECK = TestMixtureChecker()
        res = MIXCHECK.is_mixed('a12')
        self.assertFalse(res['is_mixed'])       # doesn't begin with a number
        res = MIXCHECK.is_mixed('1a2')
        self.assertTrue(res['is_mixed'])       # does begin with a number


class Test_MAL_16(unittest.TestCase):
    """ tests the setting of mixed samples in 'include' mode"""
    def runTest(self):

        n_guids = 1000              # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
                                    # this is very much a worst case scenario:  
                                    # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include')
        guids_to_add = list(p.guids())    
        m.update()

        self.assertEqual(len(m.guids()), n_guids)     
        self.assertEqual(len(m.name2meta()), n_guids)
    
        # check output after clustering
        m.cluster()   

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta() 

        for guid in res.index:
            begins_with_number = guid[0] in ['0','1','2']
            has_neighbours = res.loc[guid, 'nneighbours'] >0
            if begins_with_number and has_neighbours:
                self.assertTrue(res.at[guid,'is_mixed'])
            elif not begins_with_number and has_neighbours:
                self.assertFalse(res.at[guid,'is_mixed'])
            else:
                self.assertIsNone(res.at[guid,'is_mixed'])

        cl = m.cluster()

        # now: the guids which begin with 0-2 should be mixed.
        # the others are not mixed.
        # the unmixed samples should all be in 1 cluster.
        in_clusters = set()
        for guid in res.index:
            begins_with_number = guid[0] in ['0','1','2']
            has_neighbours = res.loc[guid, 'nneighbours'] >0
            if not begins_with_number:
                self.assertTrue(len(m.name2cluster[guid]['cluster_id']),1)      # not mixed
            else:
                in_clusters.add(len(m.name2cluster[guid]['cluster_id']))

        #    all the samples should be clustered, either on their own or separately.
        self.assertTrue(min(in_clusters)>0)


class Test_MAL_17(unittest.TestCase):
    """ tests the setting of mixed samples in 'ignore' mode"""
    def runTest(self):

        n_guids = 1000              # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
                                    # this is very much a worst case scenario:  
                                    # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'ignore')
        guids_to_add = list(p.guids())    
        m.update()

        self.assertEqual(len(m.guids()), n_guids)     
        self.assertEqual(len(m.name2meta()), n_guids)
    
        # check output after clustering
        m.cluster()   

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta() 

        # it does get scored in ignore mode
        for guid in res.index:
            begins_with_number = guid[0] in ['0','1','2']
            has_neighbours = res.loc[guid, 'nneighbours'] >0
            if begins_with_number and has_neighbours:
                self.assertTrue(res.at[guid,'is_mixed'])
            elif not begins_with_number and has_neighbours:
                self.assertFalse(res.at[guid,'is_mixed'])
            else:
                self.assertIsNone(res.at[guid,'is_mixed'])

        cl = m.cluster()

        # now: the guids which begin with 0-2 should be ignored and treated like all the others.
        # the others are not mixed.
        # the unmixed samples should all be in 1 cluster.
        in_clusters = set()
        for guid in res.index:
            begins_with_number = guid[0] in ['0','1','2']
            has_neighbours = res.loc[guid, 'nneighbours'] >0
            if not begins_with_number:
                self.assertTrue(len(m.name2cluster[guid]['cluster_id']),1)
            else:
                self.assertTrue(len(m.name2cluster[guid]['cluster_id']),1)      # each sample is in its own cluster; no cross cluster samples.
  
class Test_MAL_18(unittest.TestCase):
    """ tests the setting of mixed samples in 'exclude' mode"""
    def runTest(self):

        n_guids = 1000              # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
                                    # this is very much a worst case scenario:  
                                    # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'exclude')
        guids_to_add = list(p.guids())    
        m.update()

        self.assertEqual(len(m.guids()), n_guids)     
        self.assertEqual(len(m.name2meta()), n_guids)
    
        # check output after clustering
        m.cluster()   

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta() 

        for guid in res.index:
            begins_with_number = guid[0] in ['0','1','2']
            has_neighbours = res.loc[guid, 'nneighbours'] >0
            if begins_with_number and has_neighbours:
                self.assertTrue(res.at[guid,'is_mixed'])
                self.assertEqual(res.at[guid,'is_mixed'],m.is_mixed(guid))
            elif not begins_with_number and has_neighbours:
                self.assertFalse(res.at[guid,'is_mixed'])
            else:
                self.assertIsNone(res.at[guid,'is_mixed'])
        cl = m.cluster()

        # now: the guids which begin with 0-2 should be mixed.
        # the others are not mixed.
        # all samples should be in 1 cluster.
        in_clusters = set()
        for guid in res.index:
            begins_with_number = guid[0] in ['0','1','2']
            has_neighbours = res.loc[guid, 'nneighbours'] >0
            if not begins_with_number:
                self.assertTrue(len(m.name2cluster[guid]['cluster_id']),1)
            else:
                self.assertTrue(len(m.name2cluster[guid]['cluster_id']),1)
 

class Test_MAL_19(unittest.TestCase):
    """ tests output functions"""
    def runTest(self):

        n_guids = 10              # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
                                    # this is very much a worst case scenario:  
                                    # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include', parameters={'Param1':1,'Param2':2,'snv_threshold':12,'uncertain_base_type':'M'}, name='MAL_19')
        guids_to_add = list(p.guids())    
        m.update()

        self.assertEqual(len(m.guids()), n_guids)     
        self.assertEqual(len(m.name2meta()), n_guids)
    
        # check output after clustering
        m.cluster()   

        self.assertEqual(len(m.name2meta()), n_guids)

        # if the guid begins with a number and degree_centrality is >0, it should be mixed
        res = m.name2meta() 
        self.assertEqual(len(res.index), len(p.guids()))

        res = m.serialise_output()
        self.assertTrue(res['parameters'] is not None)
        self.assertTrue(res['guid2clustermeta'] is not None)

class Test_MAL_20(unittest.TestCase):
    """ tests persistance and recovery functions"""
    def runTest(self):

        n_guids = 10              # increase this for a performance test of everything except the mixture detection algorithm and the persistence backend
                                    # this is very much a worst case scenario:  
                                    # denovo clustering of 100k samples with ~ 60% mixtures in about 5 mins.  Algorithm is incremental, so this is a one-off cost

        p = MockPersistence(n_guids = n_guids)
 
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include', parameters={'Param1':1,'Param2':2,'snv_threshold':12,'uncertain_base_type':'M'}, name='MAL_20')
        guids_to_add = list(p.guids())    
        m.update()

        # check output after clustering
        m.cluster()   

        m.persist(what='graph')
        m.persist(what='output')
        m.remove_legacy()

        # reload from persistence
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include', parameters={'Param1':1,'Param2':2,'snv_threshold':12,'uncertain_base_type':'M'}, name='MAL_20')
        self.assertEqual(m.guids(),p.guids())

        self.assertEqual(m.guid2cluster_labels(), {})       # no labels assigned

        # apply labels; make some up
        cl2label = {}
        expected_labels = list()
        for cl in m.cluster2names:
            if len(m.cluster2names[cl]['guids'])>1:
                new_label = 'LABEL-{0}'.format(cl)
                cl2label[cl] = {'cluster_label':new_label}
                expected_labels.append(new_label)

        # apply the labels
        m._clusterid2clusterlabel = cl2label

        # check they are used correctly

        for guid in m.guid2cluster_labels():
            self.assertEqual("LABEL-{0}".format(m.name2cluster[guid]['cluster_id'][0]),m.guid2cluster_labels()[guid][0])
        self.assertEqual(set(m.existing_labels()), set(expected_labels))
        # check that the labels get persisted
        m.persist(what='graph')
        m.persist(what='output')
        m.remove_legacy()

        # reload from persistence
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include', parameters={'Param1':1,'Param2':2,'snv_threshold':12,'uncertain_base_type':'M'}, name='MAL_20')

        # test that the guid2cluster maps are the same
        self.assertEqual(m._clusterid2clusterlabel, cl2label)

        # test apply_cluster_labels()
        # check update adds remaining guids
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include', parameters={'Param1':1,'Param2':2,'snv_threshold':12,'uncertain_base_type':'M'}, name='MAL_20_v2')
        guids_to_add = list(p.guids())    
        m.update()

        # check output after clustering
        m.cluster()   
    
        # there are no labels
        self.assertEqual(m.guid2cluster_labels(), {})       # no labels assigned


        # apply labels; make some up
        cl2label = {}
        for cl in m.cluster2names:
            if len(m.cluster2names[cl]['guids'])>1:
                cl2label[cl] = {'cluster_label':'LABEL-{0}'.format(cl)}

        # apply the labels
        m.apply_cluster_labels(cl2label)

        # check they are used correctly
        for guid in m.guid2cluster_labels():
            #print(guid, m.name2cluster[guid]['cluster_id'][0], m.guid2cluster_labels()[guid])
            self.assertEqual("LABEL-{0}".format(m.name2cluster[guid]['cluster_id'][0]),m.guid2cluster_labels()[guid][0])
 
        
class test_Raise_error(unittest.TestCase):
    """ tests raise_error"""
    def runTest(self):
       p = MockPersistence(n_guids = 50)
 
       # check update adds remaining guids
       m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20)
                     
       with self.assertRaises(ZeroDivisionError):
            m.raise_error("token")

class test_MIXCHECK_1(unittest.TestCase):
    """ tests mixpore mixture checker"""
    def runTest(self):

        rc = ReadConfig()
        CONFIG = rc.read_config("../config/default_test_config.json")

        # get a clustering object's settings
        for clustering_name in CONFIG['CLUSTERING'].keys():
            clustering_setting = CONFIG['CLUSTERING'][clustering_name]

            PERSIST=fn3persistence(dbname = CONFIG['SERVERNAME'],
                connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                debug=CONFIG['DEBUGMODE']
                   )
            PERSIST._delete_existing_data()

            hc = hybridComparer(
                reference=CONFIG['reference'],
                maxNs=CONFIG['MAXN_STORAGE'],
                snpCeiling=  CONFIG['SNPCEILING'],
                excludePositions=CONFIG['excluded'],
                preComparer_parameters=CONFIG['PRECOMPARER_PARAMETERS'],
                PERSIST=PERSIST, 
                unittesting=True
                )

            #print(clustering_setting)
            mpmc = MixPOREMixtureChecker(hc, **clustering_setting) 

            # check update adds remaining guids
            m = MixtureAwareLinkage(PERSIST=PERSIST, 
                                    MIXCHECK = mpmc,
                                    mixed_sample_management = clustering_setting['mixed_sample_management'], 
                                    snv_threshold=clustering_setting['snv_threshold'])

            # add fake data
            guids_inserted = list()			
            for i in range(1,10):
                #print("Inserting",i)

                seq = list(str(CONFIG['reference']))

                if i % 3 == 0:      # every third sample is mixed
                    is_mixed = True
                    guid_to_insert = "mixed_{0}".format(i)
                else:
                    is_mixed = False
                    guid_to_insert = "nomix_{0}".format(i)	
        
                # make i mutations at position 500,000

                offset = 700000
                for j in range(i+5):
                    mutbase = offset+j
                    ref = seq[mutbase]
                    if is_mixed == False:
                        if not i % 2 ==0:
                            if not ref == 'T':
                                seq[mutbase] = 'T'
                            if not ref == 'A':
                                seq[mutbase] = 'A'
                        else:
                            if not ref == 'C':
                                seq[mutbase] = 'C'
                            if not ref == 'G':
                                seq[mutbase] = 'G'
                    if is_mixed == True:
                        seq[mutbase] = 'M'					
                seq = ''.join(seq)
                #print(i,guid_to_insert, seq[699995:700020])
                guids_inserted.append(guid_to_insert)
                obj = hc.compress(seq)	
                loginfo = hc.persist(obj, guid_to_insert)

            # test the mixporemixture checker.
            for guid in guids_inserted:
                res = mpmc.is_mixed(guid)
                #print(guid, res)
                #self.assertEqual("mixed" in guid, res['is_mixed'])      # should identify all mixed guids 

            m.update()
            m.cluster()

            res = m.name2meta()
            for guid in res.index:
                self.assertEqual("mixed" in guid, res.at[guid,'is_mixed'])      # should identify all mixed guids              
                
                # check whether it's what we expect
                if clustering_setting['mixed_sample_management'] in ['ignore', 'exclude']:
                    self.assertEqual(len(res.at[guid,'cluster_id']), 1)     # samples are in one cluster
                if clustering_setting['mixed_sample_management']=='include':
                    if 'mixed' in guid:
                        self.assertTrue(len(res.at[guid,'cluster_id'])>0)
                    else:
                        self.assertEqual(len(res.at[guid,'cluster_id']), 1)
            
class Test_MALR_1(unittest.TestCase):
    """ tests MixtureAwareLinkageResults"""
    def runTest(self):

        # generate something to analyse
        n_guids = 10  
        p = MockPersistence(n_guids = n_guids)
        m = MixtureAwareLinkage(PERSIST=p, MIXCHECK = TestMixtureChecker(), snv_threshold=20,
                                    mixed_sample_management = 'include', parameters={'Param1':1,'Param2':2,'snv_threshold':12,'uncertain_base_type':'M'}, name='test')
        guids_to_add = list(p.guids())    
        m.update()
        m.cluster()   
        serialisation = m.serialise_output()

        # create and test all methods
        malr = MixtureAwareLinkageResult(serialisation)
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            for item in res:
                keys = set(item.keys())
                self.assertEqual(keys, set(['cluster_id','cluster_label']))
            self.assertIsNotNone(res)       # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))
        self.assertTrue(isinstance(malr.snv_threshold,int))

        m.persist(what='output')
        

        malr = MixtureAwareLinkageResult(PERSIST=p, name='test')
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            self.assertIsNotNone(res)       # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))


	    # load an empty or missing data set
        malr = MixtureAwareLinkageResult(PERSIST=p, name='nodata')
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), set())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        self.assertTrue(malr.snv_threshold is None)

        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), 0)


	    # load an empty or missing data set
        p.latest_version_behaviour='nochange'       # don't increment version

        malr = MixtureAwareLinkageResult(PERSIST=p, name='test')
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            self.assertIsNotNone(res)       # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))
        t1= malr.current_version_load_time
        malr.refresh()
        t2 = malr.current_version_load_time		# should not reload
        self.assertEqual(t1,t2)


	    # load an empty or missing data set
        p.latest_version_behaviour='increment'       # increment version

        malr = MixtureAwareLinkageResult(PERSIST=p, name='test')
        self.assertTrue(isinstance(malr.change_id, int))
        self.assertTrue(isinstance(malr.guids(), set))
        self.assertEqual(malr.guids(), p.guids())
        self.assertTrue(isinstance(malr.clusters2guid(), dict))
        for guid in p.guids():
            res = malr.guid2clusters(guid)
            self.assertIsNotNone(res)       # should all be clustered
        for guid in p.guids():
            res = malr.is_mixed(guid)
            res2 = malr.is_mixed(guid, reportUnknownAsFalse=False)
            if res2 is None:
                self.assertEqual(res, False)
            else:
                self.assertEqual(res, res2)
        res = malr.clusters2guidmeta()
        self.assertEqual(len(res), len(p.guids()))
        t1= malr.current_version_load_time
        malr.refresh()
        t2 = malr.current_version_load_time
        self.assertNotEqual(t1,t2)      # updated
