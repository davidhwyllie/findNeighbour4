#!/usr/bin/env python
""" 
A component of a findNeighbour4 server which provides relatedness information for bacterial genomes.
It does so using PCA, and supports PCA based cluster generation.

The associated classes compute a variation model for samples in a findNeighbour4 server.
Computation uses data in MongoDb, and is not memory intensive, using configuration information in a 
config file. 

If no config file is provided, it will run in  'testing' mode with the  parameters
in default_test_config.json.  This expects a mongodb database to be running on
the default port on local host. 

Functionality is provided in following classes:
* VariationModel - stores results of producing variant matrix and running PCA
* VariantMatrix - computes sample x variant matrix (requires: PERSIST object for mongodb access; server configuration file) 
* PCARunner - runs PCA on VariantMatrix

These classes are not yet properly unit tested.

Example usage:
    see main()

Example of extra functionality: see '_snippet' variable

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.  See see <https://www.gnu.org/licenses/>.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

"""
 
# import libraries
import os
import io
import sys
import json
import logging
import warnings
import datetime
import glob
import time
import random
from pathlib import Path
from typing import List, Tuple, Set
from collections import defaultdict
import hashlib

import pandas as pd
import numpy as np
from Bio import Phylo

import argparse
import progressbar

import sqlalchemy

from scipy.stats import binom_test, median_abs_deviation
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, ClusterWarning
from scipy.cluster.hierarchy import ClusterWarning

from sklearn import linear_model
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import euclidean_distances

from sentry_sdk import capture_message, capture_exception
from sentry_sdk.integrations.flask import FlaskIntegration

# ignore clusterwarnings
warnings.simplefilter("ignore", ClusterWarning)     # nonsensical cluster warnings get issued

# logging
from logging.config import dictConfig

# reference based compression storage and clustering modules
from findn.mongoStore import fn3persistence


class VariationModel():
    """ Stores a VariantMatrix, the output of a PCA of the matrix, and (optionally) a clustering of the principal components.
        You should not normally have to call this class directly to create a VariationModel - the VariantMatrix class would do this for you.
        
        - You might wish to instantiate this class directly if you are restoring a previously serialised VariationModel - see constructor""" 
        
    def __init__(self,serialised_representation = None):
        """ 
           Inputs:
           serialised_representation:  either None (a new VariationModel which is empty is created) or a dictionary previously created by .serialise()


           This class can created in 3 ways:
            * an empty, if serialised representation is None, example:
            vm = VariationModel()
            - use `vm['attribute'] = value` to add properties to the model, and .finish() when the model is built.

            * a new model, generated from data.  Such models would have been created (typically with the ModelBuilder class), e.g. as follows:
            builder = ModelBuilder(CONFIG, PERSIST)
            vm = builder.build(train_on=500, n_components=15, pca_parameters= {'n_jobs':-1})       # method returns a VariationModel object
            to_serialise = json.dumps(vm.serialise())       # to_serialise is a string which can be stored, e.g. in a database.  .serialise() returns a dictionary which is jsonable.
            
            * restore  a previously build model
            deserialise = json.loads(to_serialise)          # deserialise is a dictionary
            vm2 = VariationModel(serialised_representation= deserialise)    # create a new VariationModel           

        """
        if serialised_representation is None:
            self.model = {'built':False}
        else:
            self.model = {}
            self._deserialise(serialised_representation)    
        return

    def __getitem__(self, key):
        if key not in self.model:
            raise KeyError(f"Key {key} not found")
        return self.model[key]

    def __setitem__(self, key, value):
        """ adds a key-value pair to the model """
        if key in self.model.keys():
            raise KeyError(f"Cannot replace key {key}")
        else:
            self.model[key] = value

    def _coefficients_hash(self):
        """ computes a hash on the coefficients in the variant model.
        This is useful for version tracking & storing patterns of masking. """
        h = hashlib.md5()
        h.update(self.model['eigenvectors'].to_csv().encode('utf-8'))
        md5_l = h.hexdigest()
        return("{0}".format(md5_l))

    def finish(self):
        """ completes construction of the VariationModel """
        self.model['build_time'] = datetime.datetime.now().isoformat()
        self.model['coefficients_hash'] = self._coefficients_hash()
        self.model['built'] = True

   
    def to_sqlite(self, persistdir = '', analysis_name = 'pca_output', rebuild_databases_if_present = True):
        """ write output to sqlite database 
        
        Inputs
        =======
        persistdir                         the directory the SQLite database goes into
        analysis_name                      name of the analysis.  Will become the first part of the file name
        rebuild_databases_if_present       delete any existing SQLite database 
        
        Returns
        =======
        path to sqlite database
        """

        # configure sqlite file for output.
        sqlite_file = os.path.join( persistdir,'{0}.sqlite'.format(analysis_name))
        engine = sqlalchemy.create_engine("sqlite:///{0}".format( sqlite_file, echo=True))

        # run checks on sqlite file
        if rebuild_databases_if_present:
            try:
                os.unlink( sqlite_file)
            except FileNotFoundError:
                pass

        # open connection
        conn =  engine.connect()

        metadata = []
        for key in self.model:
            print("Writing {0}".format(key))
            if not (key == 'variant_matrix' or isinstance(self.model[key], PCA)):   # we don't serialise these; one is massive and the other can't be serialised
                native_type = type(self.model[key])
                if native_type in [bool,int,  float, str]:
                    metadata.append({'variable':key, 'value': str(self.model[key]), 'native_type':str(native_type)})
                elif type(self.model[key]) in [set, list]:
                    if type(self.model[key]) == set:
                        list_data = sorted(list(self.model[key]))
                    else:
                        list_data = self.model[key]
                    tmp = pd.DataFrame(list_data, columns=[key])
                    
                    tmp.to_sql(key, conn, if_exists='fail')                # we don't serialise these at present
                elif type(self.model[key]) in [dict, defaultdict]:
                    output_records = []
                    for this_key in self.model[key].keys():
                        item = self.model[key][this_key]
                        if type(item) in [float,bool,int,str]:
                            item = [item]
                        if not type(item) == list:
                            raise TypeError("Can only export dictionaries which are of key:list or key:scalar format; the list element is of type {0} : {1}".format(type(item),item))
                        for list_element in item:
                            output_records.append({key:this_key, 'value':list_element})
                    tmp = pd.DataFrame.from_records(output_records)       

                elif type( self.model[key]) == np.int64:
                    metadata.append({'variable':key, 'value': str(int( self.model[key]))})
                elif type( self.model[key]) == pd.core.frame.DataFrame:
                    self.model[key].to_sql(key, conn, if_exists='fail')
                else:
                    warnings.warn("Not handled {0} with class {1}".format(key, type( self.model[key])))

        metadata_df = pd.DataFrame.from_records(metadata)
        metadata_df.to_sql('Metadata', conn, if_exists='fail')
        conn.close()
        return sqlite_file

    def serialise(self):
        """ serialise to model to a dictionary compatible with json serialisation
        
        Note that this may generate MemoryErrors for massive analyses.
        Persistence to database is probably a better option - see persist() method. """
      
        target_classes = dict()
        serialised = {}
        for key in vm.model:
            if not key == 'variant_matrix':                 # contains all variant sites for all isolates - really huge
                if not isinstance(vm.model[key], PCA):   # we don't serialise these
                    target_classes[key]= str(type(self.model[key]))

                    if type(self.model[key]) in (bool,int, dict, float, list, str):
                        serialised[key]=self.model[key]
                    elif type(self.model[key]) == set:
                        serialised[key]= list(self.model[key])
                    elif type(self.model[key]) == np.int64:
                        serialised[key]= int(self.model[key])
                    elif type(self.model[key]) == pd.core.frame.DataFrame:
                        serialised[key]= self.model[key].to_dict(orient='index')
                    else:
                        warnings.warn("Not handled {0} with class {1}".format(key, type(self.model[key])))

        serialised['_target_classes'] = target_classes

        return serialised


    def _deserialise(self, deserialise):
        """ transform to native format based on dictionary produced by serialise
        Called by the constructor; do not class directly.

        deserialise: a dictionary.  Should have been converted from json string if appropriate, e.g.         
                        deserialise = json.loads(to_serialise)
        """

        if not type(deserialise) == dict:
            raise TypeError("must be passed a dictionary, not a {0}.  You may need to json.loads() the string.  ".format(type(deserialise)))
        if not '_target_classes' in deserialise.keys():
            raise KeyError("The dictionary passed doesn't contain _target_classes.  It wasn't produced by .serialise, and can't be analysed")

        for key in deserialise['_target_classes'].keys():
            if not key == '_target_classes':
                if not key in deserialise:
                    raise KeyError("Expected a key {0} in the dictionary to deserialise.  It's not there.  dictionary cannot be deserialised.  Keys are {1}".format(key, deserialise.keys()))

                # restore to native types
                target_class  = deserialise['_target_classes'][key]
                if target_class in ("<class 'dict'>",  "<class 'list'>", "<class 'str'>"):
                    # no conversion needed
                    self.model[key] = deserialise[key]
                elif target_class in ("<class 'bool'>"):
                    self.model[key] = bool(deserialise[key])
                elif target_class in ("<class 'int'>"):
                    self.model[key] = int(deserialise[key])
                elif target_class in ("<class 'float'>"):
                    self.model[key] = float(deserialise[key])
                elif target_class in ("<class 'set'>"):
                    self.model[key] = set(deserialise[key])
                elif target_class in ("<class 'pandas.core.frame.DataFrame'>"):
                    self.model[key] = pd.DataFrame.from_dict(deserialise[key], orient='index')
                else:
                    raise TypeError("Unable to deserialise, don't know how to deal with {0}/{1}".format(key,deserialise[key]))       
            
    def _getNewick(self, node, newick, parentdist, leaf_names):
        """ generates a Newick format tree file 

        Inputs:
        node: a scipy hierarchy converted to a tree, e.g. below
            # dist is a distance matrix
            dist= euclidean_distances(cluster_centres, cluster_centres)
            Z = linkage(dist, 'ward')
            t = hierarchy.to_tree(Z, False)
            nwk = getNewick(t, "", t.dist, list_of_node_names)
        newick: the existing newick string
        parentdist:  the distance matrix of the parent
        leaf_names:  leaf names of the tree

        Note: code is modified from user jfn https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
        this method uses recursion.
        """

        if node.is_leaf():
            return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = self._getNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = self._getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick

    def cluster(self, eps=0.5, drop_first_pcs=0, build_tree=True):
        """ clusters the eigenvalues in the sparse PCA output using DBSCAN. 

        This technique finds regions of high point density (i.e. where there are lots of sequences)
        in the n-dimensional space represented by the eigenvalues using the DBSCAN algorithm.  
        This produces bigger and more biologically sensible clusters than the related OPTICS algorithm. 

        Ester, M., H. P. Kriegel, J. Sander, and X. Xu, “A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise”. In: 
            Proceedings of the 2nd International Conference on Knowledge Discovery and Data Mining, Portland, OR, AAAI Press, pp. 226-231. 1996
        Schubert, E., Sander, J., Ester, M., Kriegel, H. P., & Xu, X. (2017). DBSCAN revisited, revisited: why and how you should (still) use DBSCAN. 
            ACM Transactions on Database Systems (TODS), 42(3), 19.


        We use the scikit-learn implementation
        https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html

        For further discussion and benchmarking against other implementations, see
        dbscan: Fast Density-Based Clustering with R.  Michael Hahsler et al.
        Journal of Statistical Software  October 2019, Volume 91, Issue 1. doi: 10.18637/jss.v091.i0
        The scikit-learn implementation is not the fastest available, but it has good performance.

        Inputs:
            eps - a parameter controlling the minimum distance to a neighbouring point.  Values of 0.5 (default) seem appropriate.  Higher values give bigger clusters.
            The number of clusters is quite sensitive to the parameter and increases ~ monotonically for 0.005 < eps < 0.5 (not tested outside these ranges)

            drop_first_pcs:  omit the first pcs (up to this number) from the alignment.  Set to None or 0 to drop none.  Might be relevant if initial pcs were predominantly of technical origin
            build_tree - if True, generates an NJ tree from the eigenvalues of the cluster centres

        Outputs:
            Number of clusters generated

        Sets:
            self.cluster_designations: a data frame containing cluster names for each cluster
            self.cluster_centres: a data frame containging centroids of each cluster
            self.cluster_tree: a newick tree, generated by NJ method, including all the groups in the clustering

        NOTE:  this clustering method could be further parameterised to remove principal components of technical origin, if required.

        """

        # check there is a model
        if not self.model['built']:
            raise NotImplementedError("You must build a model before calling .cluster()")

        # prepare data for clustering
        if drop_first_pcs is None:
            drop_first_pcs = 0


        print("Clustering - running DBSCAN")
        drop_columns = [x for x in range(drop_first_pcs)]
        ev = self.model['eigenvalues'].copy()        # eigenvalues.  option to drop PCs of technical origin could be dropped. 
        ev.drop(ev.columns[drop_columns], axis = 1, inplace=True)

        fit = DBSCAN(metric= 'l2', min_samples=2, eps=eps).fit(ev)



        ev['cluster']=fit.labels_
        clustered = ev[ev['cluster']>=0].copy()

        # relabel clusters to include sizes, as this helps understand it later
        first_column = clustered.columns.to_list()[0]
        cluster_sizes = dict(clustered.groupby('cluster')[first_column].count())
        clustered['cluster'] = ["Cluster={0}|Size={1}".format(x, cluster_sizes[x]) for x in clustered['cluster']]

        # relabel individual samples which are numbered with negative numbers (all are labelled -1) as these are singletons
        not_clustered = ev[ev['cluster']<0].copy()
        not_clustered['cluster'] = ["Sample={0}|Size={1}".format(x,1) for x in not_clustered.index]   
        self.model['cluster_designations']= pd.concat([clustered, not_clustered])

        self.model['cluster_centres'] = self.model['cluster_designations'].groupby('cluster').mean()                  # compute centres of dimension for each cluster
        self.model['eigenvalue_clusters'] = ev
        
        # TODO:  categorise all eigenvalues, against using dbscan.
        # a contingency table of each categorised PC against all others
        # will indicate whether the PC is phylogenetically associated or not
        
        if build_tree:
            print("Building nj tree from cluster centres")
        
            dist= euclidean_distances(self.model['cluster_centres'], self.model['cluster_centres'])     # compute distance matrix
            if len(dist) > 1:       # can't build a tree with one itemprint(dist)
                # use scipy's C++ implementation; BioPython is v slow in comparison.
                print("Performing hierarchical clustering")
                Z = linkage(dist, 'ward')
                t = hierarchy.to_tree(Z, False)
                self.model['cluster_tree'] = self._getNewick(t, "", t.dist, self.model['cluster_centres'].index.to_list())
            else:
                warnings.warn('Only one cluster - no tree built')
                self.model['cluster_tree']=''

        return(len(self.model['cluster_centres'].index))
    
    def multicluster(self, drop_first_pcs=0, eps_cutoffs = [0.01,0.03,0.1,0.3,0.5,0.7,1]):
        """ clusters the eigenvalues at multiple levels, using DBSCAN.
        For more details, see .cluster()

        Inputs: 
        eps_cutoffs: a series of different values of eps to use 

        Outputs:
        self.cluster_designation: a data frame containing cluster name for each cluster
           
        """ 
        eps_cutoffs = sorted(eps_cutoffs)
        eps_cutoff_columns = []
        for eps in eps_cutoffs:
            ncl = self.cluster(drop_first_pcs = drop_first_pcs, eps = eps)
            print("Identified {0} clusters".format(ncl))
            cl = self.model['eigenvalue_clusters']['cluster'].to_frame()
            new_column_name = 'cl'+str(eps)
            eps_cutoff_columns.append(new_column_name)
            cl.rename(columns = {'cluster':new_column_name}, inplace=True)
            if eps == eps_cutoffs[0]:   # first one
                multicluster = cl
            else: 
                multicluster = pd.concat([multicluster, cl], axis=1)

        multicluster['cluster'] = ''
        composite = []
        for ix in multicluster.index: 
            composite.append('|'.join([str(x) for x in multicluster.loc[ix]]))
        
        multicluster['cluster'] = composite
        self.model['cluster_designations']=multicluster

    def annotate_tree_with_cluster_output(self, newick_tree_string):
        """ Marks the result of the most recent clustering done (either by cluster() or multicluster() on a phylogenetic tree.
        This is useful to consider whether the results of the PCA-based clustering performed here is similar to 
        that of (for example) Maximal Likelihood trees generated from the original fasta files loaded.

        For example, suppose a fasta file 'micro.fas' containing multiple genomes had been loaded.
        One could generate a ML tree using a tree building software, e.g. iqtree:

        ./iqtree2 -s micro.fas -T 6 -m GTR+I+G

        If this gave Newick format output in 'micro.fas.treefile', we could make a new version of the newick output in which
        - The cluster membership is marked in the node names
        - Whether the sample is considered mixed is marked

        Input:
            a newick string consisting of the ML tree, or any other kind of tree
            multiclusters : if True, uses the output of the most recent multicluster() call.  if False, uses the output of the most recent cluster call.

        Output:
            a modified cluster string in which any nodes clustered are marked with cluster numbers, and any ? mixed samples are marked. 
        """

        # check that clustering has occurred.
        if not 'cluster_designations' in self.model.keys():
            raise KeyError("Cannot annotate tree as not clustering is stored.  Call .cluster() on VariantModel object.")

        # load tree
        with io.StringIO(newick_tree_string) as f:
            Tree=Phylo.read(f, format='newick')

        guids = set(self.model['cluster_designations'].index)
        for i,node in enumerate(Tree.get_terminals()):
            if node.name in guids:        
                node.name  = "{0} cl {1}".format(node.name, self.model['cluster_designations'].at[node.name,'cluster'])

        badguids = set(self.model['suspect_quality_seqs'].index)
        for i,node in enumerate(Tree.get_terminals()):
            if node.name in badguids:        
                node.name  = "{0} cl {1}".format(node.name, '?mixed')


        with io.StringIO(initial_value = '') as f:
            Phylo.write(Tree, f, format='newick')
            retVal = f.getvalue()

        return retVal

    def to_excel(self, filename):
        """" exports model to excel
        
        Does so by coercion of data to Pandas.DataFrame objects and export using pd.to_excel method.
        
        Inputs:
        filename:   the name of the file to be written
        
        Returns:    nothing
        
        Raises:     user defined warning if an object which cannot be coerced to a pandas dataframe is present. 
        
        Note:       export can be very slow.  Using sqlite export as a source for excel is a better solution."""
        
        with pd.ExcelWriter(filename) as writer:
            meta_data = {}
            for key in self.model.keys():
                if not key == 'cluster_tree':       # too large for export
                    if type(self.model[key]) in [pd.DataFrame]:
                        self.model[key].to_excel(writer, sheet_name=key)
                    elif  type(self.model[key]) in [int,float,str,bool]:
                        meta_data[key]=self.model[key]
                    elif  type(self.model[key]) in [list]:
                        pd.DataFrame.from_dict({'value':self.model[key]}, orient='columns').to_excel(writer, sheet_name=key)
                    elif  type(self.model[key]) in [set]:
                        pd.DataFrame.from_dict({'value':sorted(self.model[key])}, orient='columns').to_excel(writer, sheet_name=key)
                    elif  type(self.model[key]) in [dict]:
                        pd.DataFrame.from_dict(self.model[key], orient='index').to_excel(writer, sheet_name=key)
                    elif    isinstance(vm.model[key], PCA):
                        pass    # don't export these
                    else:
                        warnings.warn("Unhandled type, cannot be written to excel: {0} {1}".format(key, type(self.model[key])))
            pd.DataFrame.from_dict(meta_data, orient='index').to_excel(writer, sheet_name = 'Metadata')


class MNStats():
    """ computes the number of M and N bases in a reference compressed object """
    def __init__(self, select_positions, analysed_reference_length):
        """ input:
            select_positions are the positions contributing to the pca model, as generated by ModelBuilder.
            analysed_reference_length are the number of reference bases analysed. """
        self.select_positions = select_positions
        self.analysed_reference_length = analysed_reference_length

    def examine(self, obj): 
        """ examines the reference compressed object obj,
        reporting
        * number of Ns and Ms in the sequence, subdivided by whether they are in 
          select_positions
        * "Test 2" (binomial test, as per findNeighbour3 paper) testing whether the frequency of Ns/Ms in the selected_positions exceed those elsewhere,
                   indicative of a mixture. """

        missing= {'M_in_model':0,
        'N_in_model':0,                     'model_positions':len(self.select_positions),                       'reference_positions':self.analysed_reference_length}

        for base in ['M','N']:      # compute missingness

            # record total numbers of N and M per guid
            try:
                missing["{0}_total".format(base)] = len(obj[base])
            except KeyError:
                missing["{0}_total".format(base)] = 0

            # examine all missing (N/M) sites, adding to a missingness model
            try:
                for pos in obj[base]:
                    if pos in self.select_positions:
                        try:
                            missing["{0}_in_model".format(base)] +=1                            
                        except KeyError:
                            pass

            except KeyError:
                pass        # if there are no M,N then we can ignore these

            # compute an factor which could be used to adjust an eigenvalue 
            # based on the number of Ns assuming Ns are missing at random.
            try:
                missing["{0}_eigenvalue_inflation_parameter".format(base)]  = len(self.select_positions)/(len(self.select_positions)-missing["{0}_in_model".format(base)])
            except KeyError:        # no N/M
                missing["{0}_eigenvalue_inflation_parameter".format(base)]  = 1

            ## do binomial test
            not_model = self.analysed_reference_length- len(self.select_positions)
            p_expected = (missing["{0}_total".format(base)]-missing["{0}_in_model".format(base)])/not_model
            missing["{0}_expected_proportion".format(base)]=p_expected
            in_model = len(self.select_positions)
            p_observed = missing["{0}_in_model".format(base)]/len(self.select_positions)
            missing["{0}_observed_proportion".format(base)]=p_observed
            p_val = binom_test(
                missing["{0}_in_model".format(base)], 
                len(self.select_positions), 
                p_expected, 
                alternative='greater')
 
            missing["{0}_p_value".format(base)]=p_val
            
        return missing


class VariantMatrix:
    """In charge of producing a sample x SNP matrix"""
    def __init__(self,CONFIG, PERSIST):
        """ Using values in CONFIG, estimates based on distances in the server with CONFIG['NAME'] on port CONFIG['PORT'].
        
        Note that the findNeighbour4 server itself does not have to be operational for this class to work, but the underlying database does have to be accessible.

        related to error handling:
        SENTRY_URL (optional)
        Note: if a FN_SENTRY URL environment variable is present, then the value of this will take precedence over any values in the config file.
        This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
        PERSIST is a storage object needs to be supplied.  The fn3Persistence class in mongoStore is one suitable object.
        PERSIST=fn4persistence(connString=CONFIG['FNPERSISTENCE_CONNSTRING'])

        """
        # store the persistence object as part of the object
        self.PERSIST=PERSIST
        
        from common_utils import validate_server_config
        required_keys=set(['IP','INPUTREF','EXCLUDEFILE','DEBUGMODE','SERVERNAME',
                           'FNPERSISTENCE_CONNSTRING', 'MAXN_STORAGE',
                            "SNPCEILING", 'MAXN_PROP_DEFAULT', 'REST_PORT',
                           'LOGFILE','LOGLEVEL','CLUSTERING'])
        validate_server_config(CONFIG, required_keys)

        # load reference
        cfg = self.PERSIST.config_read('config')
        
        # set easy to read properties from the config
        self.analysed_reference_length = len( cfg['reference']) - len(set(cfg['excludePositions']))

        # we start without any variation model
        self._reset()
        
    def _reset(self):
        """ clears existing variation model and pca result """ 
        self.vm = VariationModel()      
        self._invalid = set()   # invalid guids for which we can't compute pcs
        self.model = {'built':False, 'built_with_guids':[]}
        self.validation_data = None

    def guids(self):
        """ returns list of guids currently in the findNeighbour4 database"""
        return sorted(self.PERSIST.refcompressedsequence_guids())   # sorting is not essential, but makes more deterministic for testing
        
    def _column_name(self,pos,base):
        """ given a base at a position, returns a position:base string suitable for use as a pandas column name """
        return  f"{pos}:{base}"

    @staticmethod
    def get_position_counts(guids, num_train_on, PERSIST) -> Tuple[set, dict,dict]:
        vmodel = defaultdict(int)     # variants
        mmodel = defaultdict(int)     # missingness
        guids_analysed = set()
        bar = progressbar.ProgressBar(max_value = num_train_on)
        for num_loaded,guid in enumerate(guids):
            if num_loaded>=num_train_on:       # trained on enough samples
                break
            bar.update(num_loaded)
            refcompressed_sample = PERSIST.refcompressedsequence_read(guid) # ref compressed sequence

            if refcompressed_sample['invalid'] != 1:
                guids_analysed.add(guid)
                # for definite calls, compute variation at each position
                for base in ['A','C','G','T']:
                    var_positions = refcompressed_sample.get(base, [])
                    for var_pos in var_positions:
                        vmodel[var_pos]+=1
                # compute missingness/gaps if it's mixed (M) or N
                for base in ['M','N']:      # compute missingness/gaps if it's mixed or N
                    missingness_positions = refcompressed_sample.get(base, [])
                    for missingness_pos in missingness_positions:
                        mmodel[missingness_pos]+=1

        bar.finish()
        return guids_analysed, vmodel, mmodel


    @staticmethod
    def get_missingness_cutoff(positions: Set[int], mmodel: dict) -> int:
        missingness = map(lambda pos: mmodel.get(pos, 0), positions)
        mean_missingness = np.mean(list(missingness))
        sd_missingness = np.sqrt(mean_missingness)              # poisson assumption
        upper_cutoff = int(2*sd_missingness + mean_missingness) # normal approximation
        return upper_cutoff

    @staticmethod
    def remove_high_missingness_positions(positions: Set[int], mmodel: dict, cutoff: float) -> int:
        """Remove positions with missingness above cutoff.
        
        Returns: the number of removed positions
        """
        num_removed = 0
        for pos in mmodel.keys():
            if mmodel[pos]>cutoff and pos in positions:
                num_removed += 1
                positions.remove(pos)
        return num_removed

            
    def build(self, min_variant_freq=None, num_train_on = None, deterministic = True):
        """
            input: 
                min_variant_freq: the minimum proportion of samples with variation at that site for the site to be included.  If none, is set to 10/train_on, i.e. each variant has to appear 10 times to be considered
                num_train_on: only compute PCA on a subset of train_on samples.  Set to None for all samples.
        """
        # determine guids there in the database
        guids = self.guids()
        
        # randomise order for model training purposes if required
        if not deterministic:
            random.shuffle(guids)       # makes pipeline non-deterministic if not all samples are analysed
        else:
            guids = sorted(guids)       # keep order constant between runs
        
        if num_train_on is None:            # if we are not told how many to use then we use
            num_train_on = len(guids)       # all samples 
        
        # persist parameters used
        self.vm['num_train_on'] = num_train_on
        
        # if minimum variation is not set, only analyse variants seen at least 10 times. 
        if min_variant_freq is None:
            min_variant_freq = 10/num_train_on
        
        print(">>Determining variant sites, from a derivation set of up to {0} samples ".format(num_train_on))
        guids_analysed_stage1, vmodel, mmodel = self.get_position_counts(guids, num_train_on, self.PERSIST)

        # store variant model
        self.vm['variant_frequencies'] = vmodel
        self.vm['min_variant_freq'] = min_variant_freq
        self.vm['analysed_reference_length'] = self.analysed_reference_length
        
        # from the totality of the variation, select positions with > cutoff % variation
        cutoff_variant_number = num_train_on*min_variant_freq
        self.vm['cutoff_variant_number'] = cutoff_variant_number

        select_positions = set()
        for pos, variant_count in vmodel.items():
            if variant_count>=cutoff_variant_number:
                select_positions.add(pos)
        print("Found {0} positions which vary at frequencies more than {1}.".format(len(select_positions),min_variant_freq))

        if len(select_positions) == 0:
            raise ValueError("No variation found above cutoff. normally this is because you ran the PCA operation against an empty database; this is what happens if you omit a config file parameter, when a test database is examined by default.   Cannot continue")

        self.model['variant_positions_gt_cutoff_variant_number'] = len(select_positions)

        upper_cutoff = self.get_missingness_cutoff(select_positions, mmodel)
        print("Identified max acceptable missingness per site as {0} positions ({1}%)".format(upper_cutoff,int(100*upper_cutoff/num_train_on)))
        self.vm['max_ok_missingness'] = upper_cutoff
        self.vm['max_ok_missingness_pc'] = int(100*upper_cutoff/num_train_on)
    
        num_removed = self.remove_high_missingness_positions(select_positions, mmodel, upper_cutoff)
        self.select_positions = select_positions
        print(f"Removed {num_removed} positions with missingness > cutoff [note: we are analysing non-missing data here, which will ignore deletions]")

        print(f"Found {len(select_positions)} positions which vary at frequencies more than {min_variant_freq} and pass missingness cutoff.")
        self.vm['variant_positions_ok_missingness'] = len(select_positions)
        
        # find any samples which have a high number of missing bases
        self.mns = MNStats(select_positions, self.vm.model['analysed_reference_length'] )
        print(">> scanning for samples with unexpectedly high missingness (N), or likely to be mixed (M)")
        bar = progressbar.ProgressBar(max_value = len(guids_analysed_stage1))
        
        guid2missing = {}
        for nLoaded,guid in enumerate(guids_analysed_stage1):
            bar.update(nLoaded)
            obj = self.PERSIST.refcompressedsequence_read(guid) # ref compressed sequence
            for base in ['M','N']:      # compute how many bases in this position are either M or N

                # examine all missing (N/M) sites, adding to a missingness model
                try:
                    for pos in obj[base]:
                        if pos in select_positions:
                            try:
                                mmodel[pos]=mmodel[pos]+1
                            except KeyError:
                                if not pos in vmodel.keys():
                                    mmodel[pos]={}
                                mmodel[pos]=1

                except KeyError:
                    pass        # if there are no M,N then we can ignore these

            ## do binomial test and n count
            guid2missing[guid] = self.mns.examine(obj)

        bar.finish()


        # collate mixture quality information, and identify low quality (mixed, \
        # as judged by high Ns or Ms in the variant sites)
        mix_quality_info = pd.DataFrame.from_dict(guid2missing, orient='index')
        self.vm['mix_quality_info'] =mix_quality_info

        # identify any mixed samples.  we don't build the model from these.
        # mixed are defined as having significantly more N or M in the variant
        # positions than in other bases.

        ### CAUTION - not sure whether this is applicable to SARS-COV-2
        mix_quality_cutoff = 0.05 / len(mix_quality_info.index)     # 0.05 divided by the number of mixed samples?  Bonferroni adj; maybe better to use FDR      
        suspect_quality = mix_quality_info.query('M_p_value < {0} or N_p_value < {0}'.format(mix_quality_cutoff))
        self.vm['suspect_quality_seqs'] = suspect_quality
        print(">> Identified {0} sequences based on their having unexpectedly higher Ns in variant vs. non-variant bases; excluded from model as may be mixed.".format(len(set(suspect_quality.index.to_list()))))

        guids_analysed_stage2= guids_analysed_stage1 - set(suspect_quality.index.to_list())

        # build a variation matrix for variant sites
        vmodel = {}
        nLoaded = 0
        print(">>Gathering variation for matrix construction from {0} unmixed samples into dictionary ".format(len(guids_analysed_stage2)))

        bar = progressbar.ProgressBar(max_value = len(guids_analysed_stage2))
        
        self.model['built_with_guids'] = []
        for guid in guids_analysed_stage2:
            nLoaded+=1
            if nLoaded>=num_train_on:       # if we're told only to use a proportion, and we've analysed enough,
                break
                pass
            bar.update(nLoaded)

            obj = self.PERSIST.refcompressedsequence_read(guid) # ref compressed sequence
            # for definite calls, compute variation at each position
            
            # for invalid samples, compute nothing
            if obj['invalid'] == 1:
                self._invalid.add(guid)
            else:
                # compute a variation model - a list of bases and variants where variation occurs
                variants={}     # variation for this guid

                # for definite calls, compute variation at each position
                # positions of variation where a call was made
                for base in set(['A','C','G','T']).intersection(obj.keys()):
                    target_positions = select_positions.intersection(obj[base])
                    called_positions = dict((self._column_name(pos,base),1) for pos in target_positions) 
                    
                    variants = {**variants, **called_positions}                   
                    
                vmodel[guid]=variants
        bar.finish()

        print(">>Building variant matrix from dictionary (pandas - no progression information)")  
        vmodel = pd.DataFrame.from_dict(vmodel, orient='index')
        vmodel.fillna(value=0, inplace=True)    # if not completed, then it's reference
                                                # unless it's null, which we are ignoring at present- we have preselected sites as non-null
        # print(">>Deduplicating the original {0} sequences".format(len(vmodel.index)))
        # vmodel = vmodel.drop_duplicates()      # deduplicate by row - removes weighting by repeated isolates, but retains diversity 
        #self.vm['null_elements_in_vmodel'] = vmodel.isna().sum().sum()
        print(">>Post-deduplication there are {0} sequences".format(len(vmodel.index)))
        self.vm['variant_matrix'] = vmodel


class PCARunner():
    """ Performs PCA on a VariantMatrix """ 
        
    def __init__(self,snp_matrix: VariantMatrix):
        self.vm = snp_matrix.vm

        
    def run(self, n_components, pca_parameters={}, deterministic = True) -> VariationModel:
        """conducts pca on a snp_matrix, storing the results in the snp_matrix's VariantModel object.

            input: 
                n_components: the maximum number of components to extract.  
                pca_parameters: a dictionary of parameters passed to the PCA command 
                                Example: {'n_jobs':-1} to use all processors.  The contents of the dictionary are passed as-is to the PCA command, without any checking.
        """
        print(">>Performing pca extracting {0} components".format(n_components))
        self.vm['n_pca_components'] = n_components
            
        # if necessary, can perform incremental PCA see https://stackoverflow.com/questions/31428581/incremental-pca-on-big-data
        pca = PCA(n_components=n_components, **pca_parameters)     
        variant_matrix = self.vm["variant_matrix"]
        pca.fit(variant_matrix)
        contribs = []
        nz_columns = {}

        # summarise the positions and variants responsible for each pc
        pc2contributing_pos = {}
        contributing_basepos = set()
        contributing_pos = set()
        for i,row in enumerate(pca.components_,0):
            # mark values far from the median, which is close to zero
            row_median = np.median(row)
            row_mad = median_abs_deviation(row)
            row_upper_ci = row_median + 3*row_mad
            row_lower_ci = row_median - 3*row_mad

            pc2contributing_pos[i] = set()
            for j,cell in enumerate(row,0):
                if cell > row_upper_ci or cell < row_lower_ci:
                    outside_3mad = True
                else:
                    outside_3mad = False
                pos = int(variant_matrix.columns[j].split(':')[0])
                allele = variant_matrix.columns[j].split(':')[1]

                # indicate whether positions are strongly weighted
                if outside_3mad:
                    pc2contributing_pos[i].add(pos)
                    contributing_basepos.add(variant_matrix.columns[j])
                    contributing_pos.add(pos)
                contribs.append({'pc':i, 'pos':pos, 'allele': allele,'col':variant_matrix.columns[j], 'weight':cell, 'outside_3mad':outside_3mad})
            pc2contributing_pos[i] = sorted(list(pc2contributing_pos[i]))       # can be json serialised, unlike set
                
        # report eigenvectors which are different from median +- 3 median absolute deviations
        self.eigenvectors = pd.DataFrame.from_records(contribs)

        if len(self.eigenvectors.index)==0:
            raise KeyError("PCA problem.  No eigenvectors found.  Contributions found are as follows: {0}.  This usually means there is insufficient data to build PCs.  Try increasing the sample number".format(contribs))
        
        # compute eigenvalues for the data on which the fit was performed.
        print(">>Computing eigenvalues for {variant_matrix.index.size} unmixed samples ")
        eigenvalues_dict = {}
        for guid, evs in zip(variant_matrix.index, pca.transform(variant_matrix)):
            eigenvalues_dict[guid] = evs
        eigenvalues = pd.DataFrame.from_dict(eigenvalues_dict, orient= 'index')
        eigenvalues.columns = range(n_components)

        self.vm['pca'] = pca
        self.vm['eigenvalues'] = eigenvalues
        self.vm['eigenvectors'] = self.eigenvectors
        self.vm['explained_variance_ratio'] = list(pca.explained_variance_ratio_)
        self.vm['n_contributing_positions'] = len(contributing_pos)
        self.vm['pc2_contributing_positions'] = pc2contributing_pos
        self.vm['n_contributing_variants'] = len(contributing_basepos)
        self.vm['contributing_basepos'] = contributing_basepos
        self.vm['contributing_pos'] = contributing_pos
        self.vm['built_with_guids'] = variant_matrix.index.tolist()
        self.vm['pos_per_pc'] = [len(x) for x in self.vm.model['pc2_contributing_positions'].values()]

        print("PCA completed, identified {0} strongly contributing base/positions".format(len(contributing_basepos)))
        self.vm.finish()
        
        return self.vm
            



def main():
    # command line usage.  Pass the location of a config file as a single argument.
    parser = argparse.ArgumentParser(
        formatter_class= argparse.RawTextHelpFormatter,
        description="""Runs findNeighbour3-varmod, a service for bacterial relatedness monitoring.
                                     

Example usage: 
============== 
# show command line options 
python findNeighbour3-varmod.py --help  

# run with debug settings; only do this for unit testing.
python findNeighbour3-varmod.py     

# run using settings in myConfigFile.json.  
python findNeighbour3-varmod.py ../config/myConfigFile.json     

""")
    parser.add_argument('path_to_config_file', type=str, action='store', nargs='?',
                        help='the path to the configuration file', default=None)
    
    args = parser.parse_args()
    

    ############################ LOAD CONFIG ######################################
    print("findNeighbour4 PCA modelling .. reading configuration file.")

    from findn import DEFAULT_CONFIG_FILE
    from common_utils import read_server_config

    configFile = args.path_to_config_file
    if configFile is None:
        configFile = DEFAULT_CONFIG_FILE
        warnings.warn(f"No config file name supplied; using configuration in {DEFAULT_CONFIG_FILE}, suitable only for testing, not for production.")

    required_keys=set(['IP', 'REST_PORT', 'DEBUGMODE', 'LOGFILE', 'MAXN_PROP_DEFAULT'])
    CONFIG = read_server_config(configFile, required_keys=required_keys)

    # determine whether a FNPERSISTENCE_CONNSTRING environment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FNPERSISTENCE_CONNSTRING") is not None:
        CONFIG["FNPERSISTENCE_CONNSTRING"] = os.environ.get("FNPERSISTENCE_CONNSTRING")
        print("Set mongodb connection string  from environment variable")
    else:
        print("Using mongodb connection string from configuration file.")

    # determine whether a FN_SENTRY_URLenvironment variable is present,
    # if so, the value of this will take precedence over any values in the config file.
    # This allows 'secret' connstrings involving passwords etc to be specified without the values going into a configuraton file.
    if os.environ.get("FN_SENTRY_URL") is not None:
        CONFIG["SENTRY_URL"] = os.environ.get("FN_SENTRY_URL")
        print("Set Sentry connection string from environment variable")
    else:
        print("Using Sentry connection string from configuration file.")
        
    ########################### SET UP LOGGING #####################################  
    # create a log file if it does not exist.
    print("Starting logging")
    logdir = os.path.dirname(CONFIG['LOGFILE'])
    Path(os.path.dirname(CONFIG['LOGFILE'])).mkdir(parents=True, exist_ok=True)

    # set up logger
    loglevel=logging.INFO
    if 'LOGLEVEL' in CONFIG.keys():
            if CONFIG['LOGLEVEL']=='WARN':
                    loglevel=logging.WARN
            elif CONFIG['LOGLEVEL']=='DEBUG':
                    loglevel=logging.DEBUG

    # prepare to connection
    print("Connecting to backend data store")
    try:
            PERSIST=fn3persistence( dbname = CONFIG['SERVERNAME'],
                                    connString=CONFIG['FNPERSISTENCE_CONNSTRING'],
                                    debug=CONFIG['DEBUGMODE'],
                                    server_monitoring_min_interval_msec = 0)
    except Exception as e:
            print("Error raised on creating persistence object")
            raise


    # instantiate builder for PCA object
    rebuild = True
    print("Rebuild is {0}".format(rebuild))
    if rebuild:
        try:
            var_matrix = VariantMatrix(CONFIG, PERSIST)
        except Exception as e:
            print("Error raised on instantiating findNeighbour3 distance estimator object")
            raise

        print("Building snp matrix")
        var_matrix.build()
        print("Running PCA on snp matrix")
        pca_runner = PCARunner(var_matrix)
        pca_runner.run(n_components=200, pca_parameters={}) ## for MiniBatchSparsePCA (which gives unstable results) {'batch_size':10, 'random_state':0, 'n_jobs':-1}
        vm = pca_runner.vm
        
        print("Exporting sqlite")
        fn = vm.to_sqlite()
        
# startup
if __name__ == '__main__':
    main()

_snippet = """

    print("Exporting excel - not run as very slow.  Use sqlite as source for excel")
    #vm.to_excel('output2.xlsx')

    print("Serialising")
    to_serialise = json.dumps(vm.serialise())

    with open('test_serialisation.json', 'w') as f:
        f.write(to_serialise)

# deserialise from stored results if not built result on this run
if not rebuild:
    print('Recovering from serialised version.')
    with open('test_serialisation.json', 'r') as f:
        serialised = f.read()
        serialised = json.loads(serialised)
        vm2 = VariationModel(serialised_representation= serialised)  
else:
    vm2 = vm         


# cluster - an experimental process - may not be very useful
# too slow at large numbers with current algorithm
print("Clustering Eigenvalues")
#vm2.cluster(eps = 0.3, drop_first_pcs=4)       # just one eps
vm2.multicluster(drop_first_pcs=0)           # multiple thresholds


# write tree file
filename = '/srv/data/mixfiles/covid/milk_micro.fas.treefile'       # ML tree of a 2500 sample subset
with (open(filename,'r')) as f:
    ml_tree_string = f.read()

annotated_ml_tree_string = vm2.annotate_tree_with_cluster_output(ml_tree_string)

filename = 'output.treefile'
with (open(filename,'w')) as f:
    f.write(annotated_ml_tree_string) 

print('Finished - output in output.treefile')
exit()

## todo: VALIDATE PCS VS. tREE- BUILD SMALL TREES, SEE WHETHER PHYLOGENETICALLY COHERENT
## CONSIDER TREE BASED pc EXAMINATION

"""
