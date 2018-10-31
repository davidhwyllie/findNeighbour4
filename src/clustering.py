import networkx as nx
import datetime
import unittest
import logging
import pandas as pd

class snv_clustering():
    """ maintains clusters of samples """
    def __init__(self, saved_result=None, snv_threshold=None,  mixed_sample_management='ignore'):
        """ makes clusters of samples such that all samples can be
        reached from each other with <= snv_threshold SNV
        
        saved_result is either:
        - None
        - a dictionary, generated by .to_dict(), containing a previously stored graph and other parameters used
        to restore a snv_clustering object.
        
        If a saved_result is supplied, all the other parameters are is ignored.
        
        In a first run situation (i.e. creating a new snv_clustering object), snv 
        and mixed_sample_management are set.  They cannot be modified subsequently.
        
        snv - the clustering threshold; samples are joined if <= snv from each other.

        mixed_sample_management: dictates how mixed samples are dealt with
        'ignore': the clustering ignores the 'known_mixed' property.  This is the behaviour of standard 'snp address' and related approaches.
        'include': clusters include the guid with 'known_mixed' property,
                    but no further exploration from this guid is attempted. is reached;
                    no further guids are explored from this guid.  One guid can belong to more than one cluster. 
        'exclude':  the clusters do not include guids with 'known_mixed'=True properties.
        
        Note that in all the below documentation, 'guid' refers to a node identified by a guid,
        and 'guids' to multiple such nodes.
        """
    
        # if provided, reload the graph.        
        if isinstance(saved_result, dict):
            logging.info("Reloading saved clustering result")

            self.G =  nx.json_graph.node_link_graph(saved_result['G'])
            self.change_id = saved_result['change_id']
            self.cluster_id = saved_result['cluster_id']
            self.mixed_sample_management = saved_result['mixed_sample_management']
            self.snv_threshold = saved_result['snv_threshold']
            
        elif saved_result is None:
            logging.info("Creating new snv_clustering object")
            if snv_threshold is None:
                raise ValueError("Asked to new clustering object but snv_threshold is none.  This is not allowed.")
            self.G= nx.Graph()
            self.change_id =0
            self.cluster_id = 0
            self.mixed_sample_management = mixed_sample_management
            self.snv_threshold = snv_threshold
            if not mixed_sample_management in ['ignore','include','exclude']:
                raise ValueError("On startup, mixed_sample_management must be one of all, include, exclude")
        else:
            raise TypeError("Do not know how to reload a clustering result from an object of class {0}; a dict is expected".format(saved_result))
    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging """
        raise ZeroDivisionError(token)
    def to_dict(self):
        """ converts snv_clustering object to a dictionary.  """
        retVal = {}
        retVal['G']= nx.json_graph.node_link_data(self.G)
        retVal['change_id']=self.change_id 
        retVal['cluster_id']=self.cluster_id
        retVal['mixed_sample_management'] = self.mixed_sample_management
        retVal['snv_threshold'] = self.snv_threshold
        return retVal
    def _new_cluster_id(self):
        """ provides unused integer numbers for assignation to new clusters.
        The numbers do not automatically increase; if integer cluster_ids exist
        which are not used, they may be reused.  For example, if
        integers 1,3,4,5 are used, 2 may be assigned.
        
        note that _new_cluster_id() is new relative to the graph.  If you call _new_cluster_id() repeatedly,
        it will return the same number until the _new_cluster_id() is actually inserted.
        """
        
        new_cluster_id = None
        existing_cluster_ids = set()
        for guid, this_cluster_id_list in self.G.nodes.data('cluster_id'):
            if this_cluster_id_list is not None:
                if not isinstance(this_cluster_id_list, list):
                    raise TypeError("Expected a cluster_id list for {0}, but got {1}".format(guid, this_cluster_id_list))
                for this_cluster_id in this_cluster_id_list:
                    existing_cluster_ids.add(this_cluster_id)
        if len(existing_cluster_ids) == 0:
            return 1        # the first cluster_id allocated.
        else:
            for i in range(1, 1+max(existing_cluster_ids)):
                if not i in existing_cluster_ids:
                    return i
            return 1+ max(existing_cluster_ids)

    def set_mixed(self,guid, neighbours):
        """ marks guid as being mixed.  neighbours are the current links of guid.
        
        Since all current links of guid are not stored, only which set of samples it forms part of it,
        in order to correctly split up clusters when a mixed sample is detected, all edges have to be added first.
        
        """

        if self.is_mixed(guid):
            return      # don't have to anything - already mixed
        
        # increase the change counter
        self.change_id += 1
        
        # update cluster_ids in any linked clusters.
        in_cluster_guids = self._traverse_from(guid, how='ignore')        # scan all linked clusters
        self._update_clusterid_to_largest_cluster(in_cluster_guids)
   
        # add all the edges of the mixed guid, having removed the links we've inserted before
        # which connect the members of a cluster, but don't necessarily reflect the real links between
        # the guids.
        new_E = []
        existing_E = []
        for adjacent_node in list(self.G.adj[guid]):
            existing_E.append([guid, adjacent_node])
        self.G.remove_edges_from(existing_E)
        
        for neighbour in neighbours:
            if not self.is_mixed(neighbour):        # don't link to other mixed samples
                new_E.append((guid,neighbour))
        self.G.add_edges_from(new_E)
             
        # set is_mixed attribute
        self._change_guid_attribute(guid, 'is_mixed', True)    ## use custom function tracking history
        
        # now consider what to do next.
        if self.mixed_sample_management == 'ignore':
            # we change nothing
            return
        else:
            # define the cluster(s) which have been generated by this operation
            # to do so we traverse from each guid
            # and identify the unique sets of guids resulting.
            # we have to do this because with some kinds of clustering,
            # guids can be in multiple different clusters.
            cluster_contents = {}
            for guid in in_cluster_guids:   # traverse from each guid unless mixed
                if not self.is_mixed(guid):
                    cluster_content = frozenset(self._traverse_from(guid, how= self.mixed_sample_management))        
                else:
                    # we return a set containing only the mixed sample.
                    cluster_content = frozenset([guid])

                cluster_contents[cluster_content] = 1   # store in a dictionary to deduplicate
                                       
            # now assign each new cluster a new number
            guid2clid = {}
            for guid in in_cluster_guids:                           # assign empty list to each guid
                self._change_guid_attribute(guid, 'cluster_id',[])  # empty list; one sequence can belong to multiple clusters
    
            # assign new cluster identifiers  
            for i,cluster in enumerate(cluster_contents.keys()):
                new_cluster_id = self._new_cluster_id()
            
                for guid in cluster:
                    new_clustering = self.G.node[guid]['cluster_id']+[new_cluster_id]  # append operation
                    self._change_guid_attribute(guid, 'cluster_id',new_clustering)
        
    def _update_clusterid_to_largest_cluster(self, in_cluster_guids):
        """ set the cluster_ids of in_cluster_guids to that of
        the largest existing group of samples within in_cluster_guids
        for example, if we had five guids 1,2,3,4,5
        and guids 1,2, and 3 had a cluster_id of 1, while 4,5 had a cluster_id of 2,
        then all samples would be assigned a cluster_id of 1.
        
        If some guids are identified as possibly mixed (is_mixed == True), and so can legitimately
        belong to > 1 cluster, this function will operate  correctly.  However, in this case, it will only
        analyse guids which are not mixed when determining the cluster_ids to update.  Because of this,
        it's essential to call _update_cluster_id_to_largest_cluster BEFORE
        setting is_mixed==True.  The set_mixed() function does this automatically.
        """
              
        # assess the largest cluster, as currently recorded, within the
        # guids reachable from starting_guid.  These guids are in in_cluster_guids.
        guid2cl = nx.get_node_attributes(self.G,'cluster_id')        # links guid to cluster id for all nodes.
        cl2n = {}           # a dictionary linking cluster to number of guids within the cluster,
                            # as currently recorded in the graph.
        # compute the number of guids per cluster, excluding any mixed samples;
        for guid in in_cluster_guids:
            if not self.is_mixed(guid):
                try:
                    this_cluster_id_list = self.G.node[guid]['cluster_id']
                except KeyError:
                    this_cluster_id_list = []       # nil entered
                    
                for this_cluster_id in this_cluster_id_list:
                    if not this_cluster_id in cl2n.keys():
                        cl2n[this_cluster_id]=0
                    cl2n[this_cluster_id]+=1
        
        # find the cluster_id with the largest number of clusters;
        largest_cluster_size = max(cl2n.values())
        for new_cluster_id in sorted(cl2n.keys()):      # enforces deterministic behaviour
            if cl2n[new_cluster_id]== largest_cluster_size:
                break

        # do the update of the cluster identifiers#
        # find the existing cluster designations which are not in mixed samples, and are
        # not the new cluster_id.  These are the ones we're going to update.
        existing_cluster_designations = set(cl2n.keys()) - set([new_cluster_id])

        # update the cluster designations                    
        for guid in in_cluster_guids:
                try:
                    this_guid2cl = guid2cl[guid]        # get the current cluster

                    made_update = False
                    for i in range(len(this_guid2cl)):
                        if this_guid2cl[i] in existing_cluster_designations:
                            this_guid2cl[i] = new_cluster_id
                            made_update = True
                    if made_update:
                        this_guid2cl = set(this_guid2cl)      # ensure unique elements
                        self._change_guid_attribute(guid, 'cluster_id', sorted(list(this_guid2cl)))        ## keep deterministic; use custom function tracking history
                except KeyError:
                    # the is currently no cluster_id assigned for in_cluster_guids element guid
                    # this is because the guid doesn't exist in the graph
                    pass
                    
    def _change_guid_attribute(self,guid,attribute,value):
        """ updates a guid's attribute with a new value.
        Keeps a track of what changed, and when """
   
        # get a copy of the history of the object
        if not 'history' in self.G.node[guid].keys():
            h = [(self.change_id,"created change history post creation")]
        else:      
            h = self.G.node[guid]['history']        # read existing history
        if not attribute in self.G.node[guid].keys():
            # doesn't exist
            h.append((self.change_id, "created {0} as {1}".format(attribute,value)))
            self.G.node[guid][attribute]=value
            self.G.node[guid]['history']=h
            self.G.node[guid]['change_id']=self.change_id 
            return 
        else:
            existing_value = self.G.node[guid][attribute]
            if not existing_value == value:
                h.append((self.change_id, "changed {0} from {1} to {2}".format(attribute,existing_value, value)))
                self.G.node[guid][attribute]=value
                self.G.node[guid]['history']=h
                self.G.node[guid]['change_id']=self.change_id            
    def is_mixed(self,guid):
        """ returns True if the is_mixed attribute is set True """
        try:
            if self.G.node[guid]['is_mixed']==True:
                return True
        except KeyError:
            # there's no 'is_mixed' attribute
            pass
        return False     

    def _traverse_from(self, starting_guid, how=None):
        """ traverse the graph moving outward from a starting point, starting_guid.
        self.mixed_sample_management describes the behaviour used.
        Normally, we don't set how; it is by default set to self.mixed_sample_management.
        A few operations (such as set_mixed()) set it as part of their operation.
        
        how = 'all': the traverse occurs irrespective of whether the 'known_mixed' property is 'True'
        how = 'include': the traverse occurs until a guid with 'known_mixed' property is reached;
                        no further guids are explored from this guid.  This guid is included in the cluster.
        how = 'exclude':  the traverse does not return guids with 'known_mixed'=True properties.
        
        Note: an alternative, nx.connected_components() is built-in in networkx;
        However, this doesn't allow custom traversal (i.e. only traversing under specific circumstances)
        or starting from a single guid.
        """
        if how is None:
            how = self.mixed_sample_management
        else:
            # validate how
            if not how in ['exclude','include','ignore']:
                raise ValueError("how must be one of exclude,include or ignore.")
          
        in_cluster_guids = set()
        to_visit = set([starting_guid])
        while len(to_visit)>0:
            current_guid = to_visit.pop()
            current_guid_is_mixed = self.is_mixed(current_guid)
            
            # include this guid in the cluster unless it's mixed, and we're told to exclude such
            if not (current_guid_is_mixed and self.mixed_sample_management=='exclude'):
                in_cluster_guids.add(current_guid)
            
            # traverse onward unless mixed, and 'exclude' or 'include'[but don't traverse onward] required
            if  (current_guid_is_mixed and self.mixed_sample_management in ['include','exclude']):
                # we do not traverse further
                pass
            else:
                # traverse onwards
                for adjacent_guid in self.G.adj[current_guid]:
                    if not adjacent_guid in in_cluster_guids:       # if not already visited
                        to_visit.add(adjacent_guid)
        return(in_cluster_guids)
    
    def _minimise_edges(self, in_cluster_guids):
        """ rearrange the graph, ensuring that
         all the guids in in_cluster_guids are still reachable from each other,
         while also ensuring that the number of edges is minimised.
         This strategy keeps down the number of edges which needs to be stored - for n
         similar nodes we only have to store (n-1) edges, rather than n^2.
         
         The only exception is that any edges involving mixed samples are left untouched.
         """
         
        new_E = []
        existing_E = []

        connected_guids = sorted(list(in_cluster_guids))  # sort to make function deterministic 
        # make edges which link the nodes in the component, stringwise
        for i in range(len(connected_guids)-1):
            new_E.append((connected_guids[i], connected_guids[i+1]))
            for adjacent_node in list(self.G.adj[connected_guids[i]]):
                if not self.is_mixed(connected_guids[i]) and not self.is_mixed(adjacent_node):
                    existing_E.append([connected_guids[i], adjacent_node])

        # remove any edges from these nodes
        self.G.remove_edges_from(existing_E)
        # and replace them with a simplified version linking the same samples
        self.G.add_edges_from(new_E)
        
    def add_sample(self, starting_guid, neighbours=[]):
        """ adds a sample, guid, linked to neighbours.
        - guid should be a string
        - neighbours should be a list
        """
        # a counter used to identify the order of changes
        self.change_id +=1

        # create a list of guid - neighbour links,
        # suitable for importing into networkx,
        # from the input data  
        E = []
        for item in neighbours:
            E.append([starting_guid,item])
        self.G.add_node(starting_guid, **{'change_id':self.change_id, 'cluster_id':[self._new_cluster_id()], 'history':[(self.change_id, 'added')]})      # in it's own, new cluster.
        self.G.add_edges_from(E)
      
        # traverse, starting with guid, until we find the edge of the cluster      
        in_cluster_guids = self._traverse_from(starting_guid)
      
        # relabel clusters post addition
        self._update_clusterid_to_largest_cluster(in_cluster_guids)
        
        # minimise the number of edges, to keep the graph small but everything connected.
        self._minimise_edges(in_cluster_guids)
    
    def guids(self):
        """ returns a set of all guids in the graph """
        return set(self.G.nodes)
    
    def guid2clusters(self,guid):
        """ returns the cluster(s) a guid belongs to """
        if not guid in self.G.nodes:
            return None     # no guid
        try:
            return self.G.node[guid]['cluster_id']
        except KeyError:
            raise KeyError("Was asked to return the cluster_id of {0}; the guid exists, but does nto have a 'cluster_id' key (likely software error) {1}".format(guid, self.G.node.data()))

    def clusters2guid(self):
        """ returns a cluster -> guid mapping """
        retVal = {}
        for guid in sorted(self.G.nodes):
            try:
                for cluster_id in self.G.node[guid]['cluster_id']:
                    if not cluster_id in retVal.keys():
                        retVal[cluster_id] = []  
                    retVal[cluster_id].append(guid)
            except KeyError:
                # no cluster_id
                pass
        return retVal        

    def clusters2guidmeta(self, after_change_id=None):
        """ returns a cluster -> guid mapping """
        retVal = []
        for guid in sorted(self.G.nodes):
            for cluster_id in self.G.node[guid]['cluster_id']:
                change_id = self.G.node[guid]['change_id']
                is_mixed = self.is_mixed(guid)
                if (after_change_id is None) or (change_id > after_change_id):
                    retVal.append({'guid':guid, 'cluster_id':cluster_id,'change_id':change_id, 'is_mixed':is_mixed})
        return retVal
    
# unittests
class test_snvc_init(unittest.TestCase):
    """ tests init method of snv_clustering """
    def runTest(self):
        """ tests init """
        
        # no parameters except SNV threshold
        snvc = snv_clustering(snv_threshold=12)
        
        self.assertEqual(type(snvc.G), nx.classes.graph.Graph)
        self.assertEqual(snvc.change_id, 0)

class test_to_dict(unittest.TestCase):
    """ tests serialisation to dict """
    def runTest(self):
        snvc = snv_clustering(snv_threshold=12)
        snvc.G.add_node('a', is_mixed=True)
        snvc.G.add_node('b', is_mixed=False)
        d = snvc.to_dict()

        snvc2 = snv_clustering(saved_result=d)
        self.assertTrue(snvc.G, snvc2.G)
          
class test_set_mixed(unittest.TestCase):
    """ tests _is_mixed function """
    def runTest(self):
        # set up
        snvc = snv_clustering(snv_threshold=12)
        snvc.add_sample('c')
        self.assertFalse(snvc.is_mixed('c'))
        snvc.set_mixed('c', [])
        self.assertTrue(snvc.is_mixed('c'))
         
class test_is_mixed(unittest.TestCase):
    """ tests _is_mixed function """
    def runTest(self):
        # set up
        snvc = snv_clustering(snv_threshold=12)
        snvc.G.add_node('a', is_mixed=True)
        snvc.G.add_node('b', is_mixed=False)
        snvc.G.add_node('c')
        
        self.assertTrue(snvc.is_mixed('a'))
        self.assertFalse(snvc.is_mixed('b'))
        self.assertFalse(snvc.is_mixed('c'))
         
class test_snvc_minimise_edges(unittest.TestCase):
    """ tests edge minimisation """
    def runTest(self):
        # set up
        snvc = snv_clustering(snv_threshold=12)
        nodes = ['n1','n2','n3','n4','n5']
        E=[]
        for i in range(len(nodes)):

            for j in range(i):
                E.append((nodes[i],nodes[j]))
            snvc.G.add_node(nodes[i], cluster_id=[1])
        snvc.G.add_edges_from(E)
        
        neighbours_pre = 0
        for node in snvc.G:
            neighbours_pre = neighbours_pre + len(snvc.G.adj[node])
            
        cc_pre = list(nx.connected_components(snvc.G))   
        ## run test
        snvc._minimise_edges(nodes)
        
        # check
        neighbours_post = 0
        for node in snvc.G:
            neighbours_post = neighbours_post + len(snvc.G.adj[node])
            
        cc_post = list(nx.connected_components(snvc.G))   
            
        self.assertEqual(cc_post, cc_pre)      # still all connected
        self.assertTrue(neighbours_post < neighbours_pre)       # edges minimised
        
class test_snvc_add_sample_0(unittest.TestCase):
    """ tests allocation of new clusterids """
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12)
        
        n1 = snvc._new_cluster_id()
        snvc.G.add_node('n1', cluster_id=[n1])
        n2 = snvc._new_cluster_id()
        snvc.G.add_node('n2', cluster_id=[n2])
        n3 = snvc._new_cluster_id()
        snvc.G.add_node('n3', cluster_id=[n3])
        
        self.assertEqual(n1,1)
        self.assertEqual(n2,2)
        self.assertEqual(n3,3)

        self.assertEqual(snvc._new_cluster_id(),4)
        snvc.G.remove_node('n2')    # number 2 should become available.
        self.assertEqual(snvc._new_cluster_id(),2)

class test_snvc_update_clusterid_1(unittest.TestCase):
    """ tests updating clusterid to that of the largest cluster"""
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12)

        snvc.G.add_node('n1', cluster_id=[1])
        snvc.G.add_node('n2', cluster_id=[1])
        snvc.G.add_node('n3', cluster_id=[2])
        snvc.G.add_node('n4', cluster_id=[3])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        snvc._update_clusterid_to_largest_cluster(['n1','n2','n3'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[1],'n3':[1],'n4':[3]})

class test_traverse_1(unittest.TestCase):
    """ tests updating clusterid to that of the largest cluster"""
    def runTest(self):
        # ignore; nothing mixed      
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='ignore')

        snvc.add_sample('n1', [])
        snvc.add_sample('n2', ['n1'])
        snvc.add_sample('n3', ['n1','n2'])      # this is going to be mixed
        snvc.add_sample('n4', ['n3'])
        snvc.add_sample('n5', ['n4','n3'])
    
        res = snvc._traverse_from('n1')
        self.assertEqual(set(['n1','n2','n3','n4','n5']), res)

class test_traverse_2(unittest.TestCase):
    """ tests updating clusterid to that of the largest cluster"""
    def runTest(self):
        # ignore; n3 mixed 
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='ignore')

        snvc.add_sample('n1', [])
        snvc.add_sample('n2', ['n1'])
        snvc.add_sample('n3', ['n1','n2'])      # this is going to be mixed
        snvc.add_sample('n4', ['n3'])
        snvc.add_sample('n5', ['n4','n3'])
        snvc.set_mixed('n3', ['n1','n2','n4','n5'])
        res = snvc._traverse_from('n1')
        self.assertEqual(set(['n1','n2','n3','n4','n5']), res)


class test_traverse_3(unittest.TestCase):
    """ tests updating clusterid to that of the largest cluster"""
    def runTest(self):
        # include; n3 mixed 
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')

        snvc.add_sample('n1', [])
        snvc.add_sample('n2', ['n1'])
        snvc.add_sample('n3', ['n1','n2'])      # this is going to be mixed
        snvc.add_sample('n4', ['n3'])
        snvc.add_sample('n5', ['n4','n3'])
        snvc.set_mixed('n3', ['n1','n2','n4','n5'])
        res = snvc._traverse_from('n1')
        self.assertEqual(set(['n1','n2','n3']), res)
        res = snvc._traverse_from('n5')
        self.assertEqual(set(['n5','n4','n3']), res)

class test_traverse_4(unittest.TestCase):
    """ tests updating clusterid to that of the largest cluster"""
    def runTest(self):
        # exclude; n3 mixed 
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='exclude')

        snvc.add_sample('n1', [])
        snvc.add_sample('n2', ['n1'])
        snvc.add_sample('n3', ['n1','n2'])      # this is going to be mixed
        snvc.add_sample('n4', ['n3'])
        snvc.add_sample('n5', ['n4','n3'])
        snvc.set_mixed('n3', ['n1','n2','n4','n5'])
        res = snvc._traverse_from('n1')
        self.assertEqual(set(['n1','n2']), res)
        res = snvc._traverse_from('n5')
        self.assertEqual(set(['n5','n4']), res)

       
class test_add_sample_2(unittest.TestCase):
    """ tests insertion of new samples, where all form one cluster. """
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12)
        
        # add three
        snvc.add_sample('n1')
        self.assertTrue(len(snvc.G.nodes),1)
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1]})
        snvc.add_sample('n2')
        self.assertTrue(len(snvc.G.nodes),2)
        cluster_ids = nx.get_node_attributes(snvc.G,'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[2]})
        snvc.add_sample('n3', ['n1','n2'])
        self.assertTrue(len(snvc.G.nodes),3)
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[1], 'n3':[1]})

class test_add_sample_3(unittest.TestCase):
    """ tests insertion of new samples, where more than one cluster is present.  """
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12)
        
        # add three
        snvc.add_sample('n1')
        self.assertTrue(len(snvc.G.nodes),1)
        cluster_ids = nx.get_node_attributes(snvc.G,'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1]})
        snvc.add_sample('n2')
        self.assertTrue(len(snvc.G.nodes),2)
        cluster_ids = nx.get_node_attributes(snvc.G,'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[2]})
        snvc.add_sample('n3')
        self.assertTrue(len(snvc.G.nodes),3)
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[2],'n3':[3]})
        snvc.add_sample('n4'),
        self.assertTrue(len(snvc.G.nodes),4)
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[2],'n3':[3],'n4':[4]})
        
        snvc.add_sample('n5', ['n1','n2','n3'])
        self.assertTrue(len(snvc.G.nodes),5)
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids, {'n1':[1],'n2':[1], 'n3':[1],'n4':[4], 'n5':[1]})
        

class test_add_sample_4(unittest.TestCase):
    """ tests combining clusters, one of which is mixed.  using exclude. """
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='exclude')
        
        # add two clusters
        snvc.add_sample('n1_1')      # cluster 1, member 1

        snvc.add_sample('n2_1')     # cluster 2, members 1-3
        snvc.add_sample('n2_2', ['n2_1'])
        snvc.add_sample('n2_3', ['n2_2'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')

        self.assertEqual(cluster_ids, {'n1_1':[1],'n2_1':[2],'n2_2':[2],'n2_3':[2]})
  
        # they are linked by a mixed sample
        snvc.set_mixed('n2_2',['n2_1','n2_3'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')

        self.assertFalse(cluster_ids['n1_1']==cluster_ids['n2_1'])
        self.assertFalse(cluster_ids['n2_1']==cluster_ids['n2_3'])
        self.assertEqual(len(set(cluster_ids['n2_2']+cluster_ids['n2_1']+cluster_ids['n2_3'])),3) # all in different clusters

        snvc.add_sample('n3')
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertFalse(cluster_ids['n1_1']==cluster_ids['n2_1'])
        self.assertFalse(cluster_ids['n2_1']==cluster_ids['n2_3'])
        self.assertEqual(len(set(cluster_ids['n2_2']+cluster_ids['n2_1']+cluster_ids['n2_3'])),3) # all in different clusters

        snvc.add_sample('n4', ['n1_1','n2_1'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids['n1_1'], cluster_ids['n2_1'])
        self.assertEqual(cluster_ids['n1_1'], cluster_ids['n4'])
        self.assertFalse(cluster_ids['n2_1']==cluster_ids['n2_3'])
        self.assertEqual(len(set(cluster_ids['n2_2']+cluster_ids['n2_1']+cluster_ids['n2_3'])),3) # all in different clusters

class test_add_sample_5(unittest.TestCase):
    """ tests combining clusters, one of which is mixed.  using ignore (the mixed sample)"""
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='ignore')
        
        # add two clusters
        snvc.add_sample('n1_1')      # cluster 1, member 1

        snvc.add_sample('n2_1')     # cluster 2, members 1-3
        snvc.add_sample('n2_2', ['n2_1'])
        snvc.add_sample('n2_3', ['n2_2'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids, {'n1_1':[1],'n2_1':[2],'n2_2':[2],'n2_3':[2]})
  
        # they are linked by a mixed sample
        snvc.set_mixed('n2_2',['n2_1','n2_3'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertFalse(cluster_ids['n1_1']==cluster_ids['n2_1'])
        self.assertTrue(cluster_ids['n2_1']==cluster_ids['n2_3'])
        self.assertTrue(cluster_ids['n2_2']==cluster_ids['n2_3'])

        snvc.add_sample('n3')
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertFalse(cluster_ids['n1_1']==cluster_ids['n2_1'])
        self.assertTrue(cluster_ids['n2_1']==cluster_ids['n2_3'])

        snvc.add_sample('n4', ['n1_1','n2_1'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        self.assertEqual(cluster_ids['n1_1'], cluster_ids['n2_1'])
        self.assertEqual(cluster_ids['n1_1'], cluster_ids['n4'])
        self.assertTrue(cluster_ids['n2_1']==cluster_ids['n2_3'])
        self.assertTrue(cluster_ids['n1_1']==cluster_ids['n4'])
        self.assertEqual(len(set(cluster_ids['n2_2']+cluster_ids['n2_1']+cluster_ids['n2_3'])),1) # all in same clusters

class test_add_sample_6(unittest.TestCase):
    """ tests combining clusters, one of which is mixed.  using include (the mixed sample).
    What should happen is that we should get two sets [2_1,2_2], [2-2,2_3]"""
    def runTest(self):
               
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')
        
        # add two clusters
        snvc.add_sample('n1_1')      # cluster 1, member 1

        snvc.add_sample('n2_1')     # cluster 2, members 1-3
        snvc.add_sample('n2_2', ['n2_1'])
        snvc.add_sample('n2_3', ['n2_2'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
         
        self.assertEqual(cluster_ids, {'n1_1':[1],'n2_1':[2],'n2_2':[2],'n2_3':[2]})
  
        # they are linked by a mixed sample
        snvc.set_mixed('n2_2',['n2_1','n2_3'])
        cluster_ids = nx.get_node_attributes(snvc.G, 'cluster_id')
        
        self.assertFalse(cluster_ids['n1_1']==cluster_ids['n2_1'])
        self.assertFalse(cluster_ids['n2_1']==cluster_ids['n2_3'])
        self.assertTrue(cluster_ids['n2_3'][0] in cluster_ids['n2_2'])
        self.assertTrue(cluster_ids['n2_1'][0] in cluster_ids['n2_2'])
      
class test_guids(unittest.TestCase):
    """ tests recovery of list of guids """
    def runTest(self):
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')
        self.assertEqual(snvc.guids(),set([]))
            
        # add two samples
        snvc.add_sample('n1')      
        snvc.add_sample('n2')      
        snvc.add_sample('n3')      
        self.assertEqual(snvc.guids(),set(['n1','n2','n3']))

class test_guid2clusters(unittest.TestCase):
    """ tests recovery of list of guids """
    def runTest(self):
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')
        self.assertEqual(snvc.guids(),set([]))
            
        # add two samples
        snvc.add_sample('n1')      
        snvc.add_sample('n2')      
        snvc.add_sample('n3')      
        self.assertEqual(snvc.guids(),set(['n1','n2','n3']))

        self.assertEqual(snvc.guid2clusters('n1'), [1])

class test_guid2clusters(unittest.TestCase):
    """ tests recovery of list of guids """
    def runTest(self):
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')
        self.assertEqual(snvc.guids(),set([]))
            
        # add two samples
        snvc.add_sample('n1')      
        snvc.add_sample('n2', ['n1'])      
        snvc.add_sample('n3')      
        self.assertEqual(snvc.guids(),set(['n1','n2','n3']))

        self.assertEqual(snvc.guid2clusters('n1'), [1])
        self.assertEqual(snvc.guid2clusters('n3'), [2])
        self.assertEqual(snvc.clusters2guid(), {1:['n1','n2'], 2:['n3']})
        

class test_clusters2guidmeta(unittest.TestCase):
    """ tests recovery of list of guids """
    def runTest(self):
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')
        self.assertEqual(snvc.guids(),set([]))
            
        # add three samples
        snvc.add_sample('n1')      
        snvc.add_sample('n2', ['n1'])      
        snvc.add_sample('n3', ['n2'])      
       
        res = snvc.clusters2guidmeta()
        self.assertEqual(res,
                        [{'guid': 'n1', 'cluster_id': 1, 'change_id': 1, 'is_mixed': False},
                         {'guid': 'n2', 'cluster_id': 1, 'change_id': 2, 'is_mixed': False},
                         {'guid': 'n3', 'cluster_id': 1, 'change_id': 3, 'is_mixed': False}]
                        )

        res2 = snvc.clusters2guidmeta(after_change_id = 2)     
        self.assertEqual(res2,
                        [{'guid': 'n3', 'cluster_id': 1, 'change_id': 3, 'is_mixed': False}]
                        )
        
        snvc.set_mixed('n2', neighbours=['n1','n3'])
        
        res3 = snvc.clusters2guidmeta()
        df = pd.DataFrame.from_records(res3)
        self.assertEqual(len(df.index), 5)
        self.assertEqual(len(df.query('guid=="n2"').index),3)
        self.assertEqual(df['change_id'].tolist(), [4,4,4,4,4])

class test_Raise_error(unittest.TestCase):
    """ tests raise_error"""
    def runTest(self):
        snvc = snv_clustering(snv_threshold=12, mixed_sample_management='include')
                      
        with self.assertRaises(ZeroDivisionError):
            snvc.raise_error("token")
        
        