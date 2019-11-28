#!/usr/bin/env python3

import unittest
import copy
import pycw_client
import json
import logging
class preComparer():
    """ compares reference compressed sequences.  

    It is designed to provide either an approximate or exact method of distance estimation.

    In the exact method, it will take into account all Ms and Ns (mixed or uncertain bases).  method which does not use or store Ns, although it can do so. 
    In the approximate method, Ns are typically omitted.  This overestimates SNV distances.
    However, Omitting Ns is useful thing to do if there are a lot of Ns, as they take up a great deal of RAM and may or may not contribute to the computation of SNV distances very much.

    The key parameter altering this behaviour is 'uncertain_base'. 
    If this is 'M', then only mixed bases are considered and the distances computed are approximate.  At least for TB, this is a helpful setting.  This approach is much faster and uses much less RAM than previous alternatives.    However, it requires that parameters for the function are carefully set and validated to ensure the approximate computation detects all pairs of samples which require more detailed computation.
    if this is 'N', then only uncertian basesa re considered and the distances computed are approximate.  This not usually a sensible setting to use.
    if this is 'N_or_M', then the distances are exact and comparisons with M or N bases are considered to count nothing to the SNV distance.

    This class can use a catWalk (compiled, high speed) relatedness engine if available.  For details of catWalk, see https://gitea.mmmoxford.uk/dvolk/catWalk.


    """
    def __init__(self,                    
                    selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5,
                    uncertain_base = 'M',  
                    catWalk_url = "http://127.0.0.1:5000",
                    **kwargs
                ):

        """ instantiates the preComparer, an object which manages in-memory reference compressed sequences.  It does not manage persistence, nor does it automatically load sequences, but it does provide methods by which external applications can load data into it.
        
        Parameters:
        ===========
        selection_cutoff: (int)
            snp distances more than this are not of interest epidemiologically and are not reported

        uncertain_base (str): one of 'M', 'N', or 'N_or_M'
            which bases to regard as uncertain in the computation.  Default to M.
            if 'N_or_M', will store all us and Ns and will give the same result as the seqComparer module.

        over_selection_cutoff_ignore_factor:
            SNP distances more than over_selection_cutoff_ignore_factor * selection_cutoff do not need to be further analysed.  For example, if a SNP cutoff was 20, and over_selection_cutoff_ignore_factor is 5, we can safely consider with SNV distances > 100 (=20*5) as being unrelated.

    	for optimal performance, over_selection_cutoff_ignore_factor cutoffs need to be tuned.  Code doing so is available in preComparer_calibrator.py.

	> if over_selection_cutoff_ignore_factor is very large ( > genome_length/selection_cutoff )
	then preComparer will report all distances as requirng additional computation.
	> if over_selection_cutoff_ignore_factor is 1, then preComparer will use the estimated snp distance computed as the basis for its decision as to whether additional testing is needed.

        catWalk_url: either None, or the url of a running catwalk instance.  If None, uses python computation.  

        
        Returns:
        ========
        None

        Note:
        - to run unit tests, do
        python3 -m unittest preComparer
        """
        self.set_operating_parameters( selection_cutoff,
                    over_selection_cutoff_ignore_factor, uncertain_base)
        self.catWalk=catWalk_url
        self.uncertain_base = uncertain_base

        # initialise data structures
        self._refresh()

        # if catWalk, startup the catWalk server
        if catWalk_url is not None:
            # startup catwalk
            # TODO if it fails to start or returns inappropriately, trap the error, log it and self self.catWalk to None;
            # log a warning that we are falling back to using a python comparison engine.
            #pycw_client.init(name="catWalk_pc",
            #   reference_name="ref",
            #   reference_sequence=self.catWalk,
            #   mask_name="no_mask",
            #   mask_str = "0")

         #   pycw_client.start(cw_binary_filepath="/home/phe.gov.uk/david.wyllie/catwalk/src/cw_server",
         #                     instance_name="test",
         #                     reference_filepath="../reference/TB-ref.fasta",
         #                     mask_filepath="../reference/TB-exclude-adaptive.txt",
         #                     max_distance=20)

            print("Using catwalk server at {0}".format(catWalk_url))
            pass
        else:
            print("No catwalk server; using python based computation.")
            pass

    def set_operating_parameters(self, 
                                 selection_cutoff,
                                 over_selection_cutoff_ignore_factor, 
                                 uncertain_base
                    ):
        """ sets the parameters used by the algorithm """
        self.selection_cutoff = selection_cutoff 
        self.over_selection_cutoff_ignore_factor = over_selection_cutoff_ignore_factor
        self.uncertain_base = uncertain_base

    def check_operating_parameters(self, 
                    selection_cutoff,
                    over_selection_cutoff_ignore_factor, 
                    uncertain_base, 
                    **kwargs):
        """ returns true if the existing operating parameters are those passed, otherwise returns false """
  
        if  (self.selection_cutoff == selection_cutoff  and 
            self.over_selection_cutoff_ignore_factor == over_selection_cutoff_ignore_factor and 
            self.uncertain_base == uncertain_base ):
            return True
        else:
            return False
      
    def _refresh(self):
        """ initialise in ramn stores """
        self.seqProfile={}      	# if catWalk is not used, then the reference compressed sequence (missing Ns, depending on self.uncertain_base) are stored in here
                                    # if catWalk is used, then only the invalid key is stored.
        self.composition = {}  		# composition statistics for each sequence. Computed and stored on loading for both catWalk and python comparison engines.

    def raise_error(self,token):
        """ raises a ZeroDivisionError, with token as the message.
            useful for unit tests of error logging """
        raise ZeroDivisionError(token)
    
    def persist(self, obj, guid):
        """ keeps a reference compressed object in the preComparer.

            if the catWalk engine is being used, then the sequences are loaded into catWalk.

            Parameters:
                obj: 
                a reference compressed object, which is a dictionary of the form
                {'A':set([1,2,3,4]} which would represent non-reference A at positions 1-4.
                The following keys are possible:
                A,C,G,T,N,M, invalid.
                A,C,G,T,N keys map to sets representing the positions where the reference differs by the respective keys 
                M is a dictionary e.g. {1:'y',2:'R'} where the keys of the dictionary represent the positions where the base is mixed and the value is the IUPAC code for the mixture.
                invalid is 0 of the sequence is valid (i.e. comparisons should be made) and 1 if it is not (no comparison are made with the sequence).
                This format is that used and stored by findNeighbour4.

                guid:
                a unique identifier for a sequence

            Returns:
                invalid: 0 if valid, 1 if invalid
            """

        # check if the guid exists.  Does not allow overwriting
        if guid in self.seqProfile.keys():
            return self.seqProfile[guid]['invalid']
            #raise KeyError("Duplicate guid supplied; already exists in preComparer: {0}".format(guid))


        # check it is not invalid
        isinvalid = False
        try:
            isinvalid = (obj['invalid']==1)
        except KeyError:
            # no invalid tag
            obj['invalid']=0

        # construct an object to be stored either directly in the self.seqProfile or in catWalk

        # determine the object's composition, i.e. numbers of bases differing from reference.
        # this is not required for SNV computation, but it is required for composition-based mixture testing (mixPORE etc).
        # computing this here is very lightweight computationally
        obj_composition= {'A':0,'C':0,'G':0,'T':0,'M':0,'N':0, 'invalid':0}
        if isinvalid:
            obj_composition['invalid']=1
  
        # construct an object to store
        for key in set(['A','C','G','T']) - set(obj.keys()):        # what is missing
                obj[key]=set()      # add empty set if no key exists.

        # create a smaller object to store.
        smaller_obj = {}
        for item in ['A','C','G','T']:
            try:
                smaller_obj[item] = copy.deepcopy(obj[item])
            except KeyError:
                pass        # if it doesn't exist, that's OK	

        # store uncertain bases - either us, Ns, or both
        # if there are no M or N keys, add them
        if not 'N' in obj.keys():
                obj['N']=set()
        if not 'M' in obj.keys():
                obj['M']=set()
        else:
                obj['M']=set(obj['M'].keys())       # make a set of the positions of us

        if self.uncertain_base == 'M':
            obj['U']=obj['M']
        elif self.uncertain_base == 'N':
            obj['U']=obj['N']
        elif self.uncertain_base == 'N_or_M':
            obj['U']=obj['N'].union(obj['M'])
        else:
            raise KeyError("Invalid uncertain_base: got {0}".format(self.uncertain_base))

        # add the uncertain bases to the smaller object for storage
        smaller_obj['U']=copy.deepcopy(obj['U'])
        key_mapping= {'A':'A','C':'C','G':'G','T':'T','U':'N'} # catwalk uses N, not U, for unknown.
        if self.catWalk is not None:
            if not isinvalid:       # we only store valid sequences in catWalk
                to_catwalk = {}
                for key in smaller_obj.keys():                    
                    to_catwalk[key_mapping[key]]=list(smaller_obj[key])
                     
                # TODO  trap all errors and raise appropropriate error 
                rcs_json = json.dumps(to_catwalk)
                pycw_client.add_sample_from_refcomp(guid, rcs_json)             
                #logging.debug("CATWALK INPUT|{0}|{1}".format(guid,rcs_json))          # log what we wrote to catwalk 
            self.seqProfile[guid]={'invalid':obj['invalid']}      # that's all we store in python if catWalk is in use
 
        else:           # cache in python dictionary
            smaller_obj['invalid']=obj['invalid']
            smaller_obj['U']=copy.deepcopy(obj['U'])
            self.seqProfile[guid]=smaller_obj           # store in python dictionary for comparison

        # set composition
        for key in set(obj.keys()).intersection(['A','C','G','T','M','N','U']):
                obj_composition[key] = len(obj[key])                    
        self.composition[guid] = obj_composition

        return obj['invalid']

    def remove(self, guid):
        """ removes a reference compressed object into RAM.
            If compression relative to other sequences has been carried out post-hoc in ram,
            only the sequence is removed; any consensus linked to it (and potentially to other sequences)
            remain unaltered.
            """

        try:
               del self.seqProfile[guid]
               del self.composition[guid]
        except KeyError:
               pass 	# we permit attempts to delete things which don't exist

        if self.catWalk is not None:
            # we cannot currently delete from catwalk
            pass

    def guids(self):
        """ returns the sample identifiers in ram in the preComparer"""
        return set(self.seqProfile.keys())
    def _classify(self, key1, key2, dEstim):
        """ returns a dictionary classifying our reponse to the observed dEstim between key1 and key2"""
        res= {
        	'guid1_invalid':self.seqProfile[key1]['invalid'],
			'guid2_invalid':self.seqProfile[key2]['invalid'],
			'guid1_Ms':self.composition[key1]['M'],
			'guid2_Ms':self.composition[key2]['M'],
			'guid1_Ns':self.composition[key1]['N'],
			'guid2_Ns':self.composition[key2]['N'],
			'guid1_Us':self.composition[key1]['U'],
			'guid2_Us':self.composition[key2]['U'],
            'no_retest':False,
			'reported_category':'not assigned'
        }

        if self.uncertain_base == 'N_or_M':
            # then distances are exact
            res['dexact'] = dEstim
            res['no_retest'] = True
            return res


        if dEstim is None:          # invalid
            res['reported_category'] = 'No retest - invalid'
            res['no_retest'] = True

        elif dEstim > self.over_selection_cutoff_ignore_factor * self.selection_cutoff:		# very large          
                res['reported_category'] = 'No retest - High estimated distance'
                res['no_retest'] = True

        else:
                res['reported_category'] = 'Retest - Estimated distance small'
                res['no_retest'] = False

        return res

    def mcompare(self, guid, guids=None):
        """ performs comparison of one guid with 
        all guids, which are all stored samples stored with .persist().
        """

        # if guids are not specified, we do all vs all
        if guids is None:
            guids = set(self.seqProfile.keys())     # guids are all guids which exist in the engine.  TODO consider if catWalk would want to do this differnetly
        
        if not guid in self.seqProfile.keys():      # seqProfile keys do include all the samples stored in catWalk + the invalid ones; what we're asked for should be in here
            raise KeyError("Asked to compare {0}  but guid requested has not been stored.  call .persist() on the sample to be added before using mcompare.")
        
        if self.seqProfile[guid]['invalid']==1:     # sequence is invalid
            return []                               # no neighbours reported

        # now do the guid vs. guids comparison
        neighbours = []
        if self.catWalk is None:
            # use python comparison engine
            guids = list(set(guids))       
            sampleCount = len(guids)
            
            for key2 in guids:
                if not guid==key2:
                    res=self.compare(guid,key2)

                    res.update(self._classify(guid,key2, res['destim']))       
                    neighbours.append(res)     
        else:

            # use catwalk
            sample_neighbours= pycw_client.neighbours(guid)
            for (neighbour, dist) in sample_neighbours:
                if neighbour in guids:
                    res = {'guid1':guid, 'guid2':neighbour, 'destim':dist}
                    res.update(self._classify(guid, neighbour, dist))       
                    neighbours.append(res)

        return neighbours

    def summarise_stored_items(self):
        """ counts how many sequences exist of various types """
        retVal = {}
        retVal['server|pcstat|nSeqs'] = len(self.seqProfile.keys())
  
        # call the catWalk server, and return the dictionary including status information in key-value format.
        # the keys must be pipe delimited and should be for the format 
        # 'server|catWalk|{key}'.  The values must be scalars, not lists or sets.
        # the dictionary thus manipulated should be added to retVal.
        if self.catWalk is not None:
            cw_status=pycw_client.info()
            for item in cw_status:
                key = "server|catwalk|{0}".format(item)
                retVal[key]=cw_status[item]

        return(retVal)

    def iscachedinram(self,guid):
        """ returns true or false depending whether we have a local copy of the refCompressed representation of a sequence (name=guid) in the seqComparer """
        if guid in self.seqProfile.keys():
            return(True)
        else:
            return(False)

    def guidscachedinram(self):
        """ returns all guids with sequence profiles currently in this machine """
        retVal=set()
        for item in self.seqProfile.keys():
            retVal.add(item)
        return(retVal)


    def compare(self,key1, key2):
        """ compares the sequences identified by key1, key2 pairwise, using python. 
        """

        # determine whether we can 'exit fast' - stop computation at an early stage
        # optimisation discussed with Denis Volk
        # assign default return values
        res = {'guid1':key1, 
			'guid2':key2, 
			'destim':None}

        # if it is invalid, we can exit fast
        if self.seqProfile[key1]['invalid']+self.seqProfile[key2]['invalid']>0:
            dEstim=None
        else:
            # otherwise we do the computation
            dEstim = self._compare_python(key1,key2)        # estimated distance
        res['destim']=dEstim
        return res


  
    def _compare_python(self,key1, key2):
        """ compares two seqProfiles.  Called when self.catWalk is None.
            Parameters:

            key1: a guid
            key2: a guid
            
            Returns: pairwise estimated SNV distance (integer)

        """

        # check keys are present
        if not key1 in self.seqProfile.keys():
            raise KeyError("{0} not present in preComparer".format(key1))
        if not key2 in self.seqProfile.keys():
            raise KeyError("{0} not present in preComparer".format(key2))

	    # compute us in non-reference positions
        try:	
            u1 = set(list(self.seqProfile[key1]['U']))
        except KeyError:
            u1 = set()
        try:
    	    u2 = set(list(self.seqProfile[key2]['U']))
        except KeyError:
            u2 = set()

        us = u1 | u2
	    
        # compute positions which differ;
        nDiff=0

        differing_positions = set()
        for nucleotide in ['C','G','A','T']:

            # we do not consider differences relative to the reference if the other nucleotide is an M
            nonU_seq1=self.seqProfile[key1][nucleotide]-us
            nonU_seq2=self.seqProfile[key2][nucleotide]-us
            differing_positions = differing_positions | (nonU_seq1 ^ nonU_seq2)
            
        us= differing_positions.intersection(us)
        n_us = len(us)

        destim = len(differing_positions)       # number of positions varying
        return destim

 

class test_preComparer_1(unittest.TestCase):
    """ tests __init__ method"""
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor= 5)

class test_preComparer_1a(unittest.TestCase):
    """ tests __init__ method with catwalk"""
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor= 5,
                    catWalk_url = "AAAAAAAA")


class test_preComparer_1b(unittest.TestCase):
    """ tests check_operating_parameters method"""
    def runTest(self):
        # initialise comparer
        sc=preComparer( selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5, uncertain_base='M')

        sc.set_operating_parameters(selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5, uncertain_base='M')

        self.assertTrue(sc.check_operating_parameters(selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5, uncertain_base='M'))

        sc.set_operating_parameters(  selection_cutoff = 10,
                    over_selection_cutoff_ignore_factor = 5, uncertain_base='M'
                    )

        self.assertFalse(sc.check_operating_parameters(selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5, uncertain_base='M'))

class test_preComparer_2(unittest.TestCase):
    """ tests storage """
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20)
        obj = {'A':set([1,2,3,4])}
        sc.persist(obj,'guid1')
        self.assertEqual(sc.composition['guid1']['A'],4)
        self.assertEqual(sc.composition['guid1']['C'],0)
        self.assertEqual(sc.composition['guid1']['T'],0)
        self.assertEqual(sc.composition['guid1']['G'],0)
        self.assertEqual(sc.composition['guid1']['N'],0)
        self.assertEqual(sc.composition['guid1']['M'],0)
        self.assertEqual(sc.composition['guid1']['invalid'],0)

class test_preComparer_2a(unittest.TestCase):
    """ tests storage """
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20, catWalk_url="AAAA")
        obj = {'A':set([1,2,3,4])}
        sc.persist(obj,'guid1')
        self.assertEqual(sc.composition['guid1']['A'],4)
        self.assertEqual(sc.composition['guid1']['C'],0)
        self.assertEqual(sc.composition['guid1']['T'],0)
        self.assertEqual(sc.composition['guid1']['G'],0)
        self.assertEqual(sc.composition['guid1']['N'],0)
        self.assertEqual(sc.composition['guid1']['M'],0)
        self.assertEqual(sc.composition['guid1']['invalid'],0)

class test_preComparer_3(unittest.TestCase):
    """ tests storage of invalid samples """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20)
        obj = {'A':set([1,2,3,4]), 'invalid':1}
        sc.persist(obj,'guid1')
        self.assertEqual(sc.composition['guid1']['A'],4)
        self.assertEqual(sc.composition['guid1']['C'],0)
        self.assertEqual(sc.composition['guid1']['T'],0)
        self.assertEqual(sc.composition['guid1']['G'],0)
        self.assertEqual(sc.composition['guid1']['N'],0)
        self.assertEqual(sc.composition['guid1']['M'],0)
        self.assertEqual(sc.composition['guid1']['invalid'],1)

        self.assertEqual(set(sc.seqProfile['guid1'].keys()), set(['invalid','A','T','C','U','G']))
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
        self.assertEqual(sc.composition['guid2']['A'],4)
        self.assertEqual(sc.composition['guid2']['C'],0)
        self.assertEqual(sc.composition['guid2']['T'],0)
        self.assertEqual(sc.composition['guid2']['G'],0)
        self.assertEqual(sc.composition['guid2']['N'],0)
        self.assertEqual(sc.composition['guid2']['M'],0)
        self.assertEqual(sc.composition['guid2']['invalid'],0)

        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')
        self.assertEqual(sc.composition['guid3']['A'],4)
        self.assertEqual(sc.composition['guid3']['C'],0)
        self.assertEqual(sc.composition['guid3']['T'],0)
        self.assertEqual(sc.composition['guid3']['G'],0)
        self.assertEqual(sc.composition['guid3']['N'],2)
        self.assertEqual(sc.composition['guid3']['M'],0)
        self.assertEqual(sc.composition['guid3']['invalid'],0)

        obj = {'A':set([1,2,3,4]), 'N':set([10,11,12]), 'M':{13:'Y'}, 'invalid':0}
        sc.persist(obj,'guid4')
        self.assertEqual(sc.composition['guid4']['A'],4)
        self.assertEqual(sc.composition['guid4']['C'],0)
        self.assertEqual(sc.composition['guid4']['T'],0)
        self.assertEqual(sc.composition['guid4']['G'],0)
        self.assertEqual(sc.composition['guid4']['N'],3)
        self.assertEqual(sc.composition['guid4']['M'],1)
        self.assertEqual(sc.composition['guid4']['invalid'],0)

        self.assertEqual(sc.guidscachedinram(),set(['guid1','guid2','guid3','guid4']))

        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['A','C','T','G','U','invalid']))
        self.assertEqual(set(sc.seqProfile['guid1'].keys()), set(['A','C','T','G','U','invalid']))

class test_preComparer_3a(unittest.TestCase):
    """ tests storage of invalid samples """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20, catWalk_url = "AAAA")
        obj = {'A':set([1,2,3,4]), 'invalid':1}
        sc.persist(obj,'guid1')
        self.assertEqual(sc.composition['guid1']['A'],4)
        self.assertEqual(sc.composition['guid1']['C'],0)
        self.assertEqual(sc.composition['guid1']['T'],0)
        self.assertEqual(sc.composition['guid1']['G'],0)
        self.assertEqual(sc.composition['guid1']['N'],0)
        self.assertEqual(sc.composition['guid1']['M'],0)
        self.assertEqual(sc.composition['guid1']['invalid'],1)

        self.assertEqual(set(sc.seqProfile['guid1'].keys()), set(['invalid']))

        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
        self.assertEqual(sc.composition['guid2']['A'],4)
        self.assertEqual(sc.composition['guid2']['C'],0)
        self.assertEqual(sc.composition['guid2']['T'],0)
        self.assertEqual(sc.composition['guid2']['G'],0)
        self.assertEqual(sc.composition['guid2']['N'],0)
        self.assertEqual(sc.composition['guid2']['M'],0)
        self.assertEqual(sc.composition['guid2']['invalid'],0)

        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')
        self.assertEqual(sc.composition['guid3']['A'],4)
        self.assertEqual(sc.composition['guid3']['C'],0)
        self.assertEqual(sc.composition['guid3']['T'],0)
        self.assertEqual(sc.composition['guid3']['G'],0)
        self.assertEqual(sc.composition['guid3']['N'],2)
        self.assertEqual(sc.composition['guid3']['M'],0)
        self.assertEqual(sc.composition['guid3']['invalid'],0)

        obj = {'A':set([1,2,3,4]), 'N':set([10,11,12]), 'M':{13:'Y'}, 'invalid':0}
        sc.persist(obj,'guid4')
        self.assertEqual(sc.composition['guid4']['A'],4)
        self.assertEqual(sc.composition['guid4']['C'],0)
        self.assertEqual(sc.composition['guid4']['T'],0)
        self.assertEqual(sc.composition['guid4']['G'],0)
        self.assertEqual(sc.composition['guid4']['N'],3)
        self.assertEqual(sc.composition['guid4']['M'],1)
        self.assertEqual(sc.composition['guid4']['invalid'],0)

        self.assertEqual(sc.guidscachedinram(),set(['guid1','guid2','guid3','guid4']))
        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['invalid']))
        self.assertEqual(set(sc.seqProfile['guid1'].keys()), set(['invalid']))
        
class test_preComparer_4(unittest.TestCase):
    """ tests reporting of server status """
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5)
        res1 = sc.summarise_stored_items()
        self.assertEqual({'server|pcstat|nSeqs': 0 }, res1)

        obj = {'A':set([1,2,3,4])}
        sc.persist(obj,'guid1')
        res2 = sc.summarise_stored_items()
        self.assertEqual({'server|pcstat|nSeqs': 1}, res2)

class test_preComparer_4a(unittest.TestCase):
    """ tests reporting of server status """
    def runTest(self):
        # initialise comparer
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5, catWalk_url = "AAAA")
        res1 = sc.summarise_stored_items()
        self.assertEqual(res1['server|pcstat|nSeqs'], 0 )
        self.assertEqual(res1['server|catwalk|mask_name'], 'no_mask' )

class test_preComparer_5(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')

        res = sc.compare('guid2','guid3')


class test_preComparer_6(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([3,4,5,6]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')

        res = sc.compare('guid2','guid3')

class test_preComparer_7(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5   
                    )
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([3,4,5,6]), 'N':set([10,11]), 'invalid':1}
        sc.persist(obj,'guid3')
 
        # check invalid set
        self.assertEqual(sc.seqProfile['guid3']['invalid'], 1)
        self.assertEqual(sc.seqProfile['guid2']['invalid'], 0)
              
        res = sc.compare('guid2','guid3')
        self.assertIsNone(res['destim'])

class test_preComparer_8(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        # check all keys are present
        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['U','A','C','G','T','invalid']))
       
        obj = {'A':set([1,2,3,4]), 'invalid':1}
        sc.persist(obj,'guid3')
      
        # check all keys are present
        self.assertEqual(set(sc.seqProfile['guid3'].keys()), set(['U','A','C','G','T','invalid']))

        # check invalid set
        self.assertEqual(sc.seqProfile['guid3']['invalid'], 1)
 
class test_preComparer_9(unittest.TestCase):
    """ tests mcompare """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5)
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        # check all keys are present
        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['U','A','C','G','T','invalid']))
       
        obj = {'A':set([1,2,3,4,5]), 'invalid':0}
        sc.persist(obj,'guid3')
      
        res = sc.mcompare('guid2')
        self.assertEqual(len(res),1)

class test_preComparer_9a(unittest.TestCase):
    """ tests mcompare """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5, catWalk_url = "TTTTTT")
       
        obj = {'A':set([1,2,3,4]), 'invalid':0}
        sc.persist(obj,'guid2')
      
        # check only invalid key is present
        self.assertEqual(set(sc.seqProfile['guid2'].keys()), set(['invalid']))
       
        obj = {'A':set([1,2,3,4,5]), 'invalid':0}
        sc.persist(obj,'guid3')
      
        res = sc.mcompare('guid2')
        self.assertEqual(len(res),1)

class test_preComparer_10(unittest.TestCase):
    """ tests comparison """
    def runTest(self):
        
        sc=preComparer(  selection_cutoff = 20,
                    over_selection_cutoff_ignore_factor = 5)
       
        obj = {'invalid':1}
        sc.persist(obj,'guid2')
      
        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        sc.persist(obj,'guid3')

        obj = {'A':set([1,2,3,4]), 'N':set([10,11]), 'invalid':0}
        res = sc.persist(obj,'guid2')     # should return results for the previously stored guid2

        self.assertEqual(res, 1)

