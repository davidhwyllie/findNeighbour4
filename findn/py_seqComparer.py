#!/usr/bin/env python3
""" performs reference based compression and comparisons of SNV distance using these 

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@phe.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.

 
"""
# python3 code to compare fasta sequences
import hashlib
import json
import numpy as np
from scipy.stats import binom_test
import datetime
import pandas as pd
from collections import Counter
from findn.msa import MSAResult


class py_seqComparer:
    def __init__(
        self,
        reference,
        maxNs,
        snpCeiling,
        excludePositions=set(),
        return_none_if_high_snp_distance=True,
    ):

        """instantiates the sequence comparer, an object which manages reference compressed sequences.

        It does not manage persistence, nor does it automatically load sequences.

        reference is a string consisting of the reference sequence.
        This is required because, as a data compression technique,
        only differences from the reference are stored.

        excludePositions contains a zero indexed set of bases which should not be considered at all in the sequence comparisons.
        Any bases which are always N should be added to this set.
        Not doing so will substantially degrade the algorithm's performance.

        self.return_none_if_high_snp_distance: if True (default), does not report high SNP distances; returns None if the distance is higher than a SNP cutoff

        If the number of Ns are more than maxNs, no data from the sequence is stored.


        Results > snpCeiling are not returned or stored.

        unknown_base_type is either N or M, and is used for computation of mixture statistics
        David Wyllie, Nov 2018

        - to run unit tests, do
        python3 -m unittest py_seqComparer
        """

        # we support three kinds of sequences.
        # sequence in strings;
        # reference based compression relative to reference 'compressed_sequence';
        self.return_none_if_high_snp_distance = return_none_if_high_snp_distance
        self.compressed_sequence_keys = set(
            ["invalid", "A", "C", "G", "T", "N", "M", "U"]
        )
        self.snpCeiling = snpCeiling

        # sequences with more than maxNs Ns will be considered invalid and their details (apart from their invalidity) will not be stored.
        self.maxNs = maxNs

        # check composition of the reference.
        self.reference = str(reference)  # if passed a Bio.Seq object, coerce to string.
        self.reference_list = list(self.reference)  # a list, one element per nt
        letters = Counter(self.reference)
        if len(set(letters.keys()) - set(["A", "C", "G", "T"])) > 0:
            raise TypeError(
                "Reference sequence supplied contains characters other than ACTG: {0}".format(
                    letters
                )
            )

        # load the excluded bases
        self.excluded = excludePositions

        # define what is included
        self.included = set(range(len(self.reference))) - self.excluded

        # initialise pairwise sequences for comparison.
        self._refresh()

    def persist(self, object, guid):
        """keeps a reference compressed object into RAM.
        Note: the sequences are stored on disc/db relative to the reference.
        Compression relative to each other is carried out post-hoc in ram
        """

        # older databases don't store the M/N combination positions in the 'U' key
        # we create it on load into RAM
        if "U" not in object.keys():
            object["U"] = set()
            if "N" in object.keys():
                object["U"] = object[
                    "N"
                ]  # if no Ms, then this is a reference taking ~ no memory
            if "M" in object.keys():
                if len(object["M"]) > 0:
                    object["U"] = object["N"].intersection(object["M"])

        self.seqProfile[guid] = object

    def remove(self, guid):
        """removes a reference compressed object into RAM.
        If compression relative to other sequences has been carried out post-hoc in ram,
        only the sequence is removed; any consensus linked to it (and potentially to other sequences)
        remain unaltered.
        """
        try:
            del self.seqProfile[guid]

        except KeyError:
            pass  # we permit attempts to delete things which don't exist

    def load(self, guid):
        """returns a reference compressed object into RAM.
        Note: this function loads stored on disc/db relative to the reference.
        Compression relative to each other is carried out post-hoc in ram
        """
        return self.seqProfile[guid]

    def remove_all(self):
        """empties any sequence data from ram"""
        self._refresh()

    def _refresh(self):
        """empties any sequence data from ram"""
        self.seqProfile = {}

    def mcompare(self, guid, guids=None):
        """performs comparison of one guid with
        all guids, which are also stored samples.
        """

        # if guids are not specified, we do all vs all
        if guids is None:
            guids = set(self.seqProfile.keys())

        if guid not in self.seqProfile.keys():
            raise KeyError(
                "Asked to compare {0}  but guid requested has not been stored.  call .persist() on the sample to be added before using mcompare.".format(
                    guid
                )
            )

        guids = list(set(guids))
        neighbours = []

        for key2 in guids:
            if not guid == key2:
                (
                    guid1,
                    guid2,
                    dist,
                    n1,
                    n2,
                    nboth,
                    N1pos,
                    N2pos,
                    Nbothpos,
                ) = self.countDifferences_byKey(
                    keyPair=(guid, key2), cutoff=self.snpCeiling
                )
                neighbours.append(
                    [guid1, guid2, dist, n1, n2, nboth, N1pos, N2pos, Nbothpos]
                )

        return neighbours

    def distmat(self, half=True, diagonal=False):
        """returns a distance matrix.
        parameters
        If half=True, returns only half the matrix.
        If diagonal=True, return the diagonal (all of which are zeros)

        returns:
        A generator yielding neighbours, in format [guid1,guid2,dist,n1,n2,nboth,N1pos, N2pos, Nbothpos]
        """

        for guid1 in self.seqProfile.keys():
            for guid2 in self.seqProfile.keys():
                include = False
                if not half and not guid1 == guid2:
                    include = True
                if guid1 == guid2 and diagonal:
                    include = True
                if guid1 < guid2 and half:
                    include = True

                if include:
                    yield self.countDifferences_byKey(keyPair=(guid1, guid2))

    def summarise_stored_items(self):
        """counts how many sequences exist of various types"""
        retVal = {}
        retVal["server|scstat|nSeqs"] = len(self.seqProfile.keys())

        return retVal

    def iscachedinram(self, guid):
        """returns true or false depending whether we have a local copy of the refCompressed representation of a sequence (name=guid) in this machine"""
        if guid in self.seqProfile.keys():
            return True
        else:
            return False

    def guidscachedinram(self):
        """returns all guids with sequence profiles currently in this machine"""
        retVal = set()
        for item in self.seqProfile.keys():
            retVal.add(item)
        return retVal

    def _delta(self, x):
        """returns the difference between two numbers in a tuple x"""
        return x[1] - x[0]

    def excluded_hash(self):
        """returns a string containing the number of nt excluded, and a hash of their positions.
        This is useful for version tracking & storing patterns of masking."""
        list_excluded = sorted(list(self.excluded))
        len_l = len(list_excluded)
        h = hashlib.md5()
        h.update(json.dumps(list_excluded).encode("utf-8"))
        md5_l = h.hexdigest()
        return "Excl {0} nt [{1}]".format(len_l, md5_l)

    def uncompress(self, compressed_sequence):
        """returns a sequence from a compressed_sequence"""
        if "invalid" in compressed_sequence.keys():
            if compressed_sequence["invalid"] == 1:
                raise ValueError(
                    "Cannot uncompress an invalid sequence, because the sequence it is not stored {0}".format(
                        compressed_sequence.keys()
                    )
                )

        seq = list(self.reference)

        # mark all positions excluded as N
        for x in self.excluded:
            seq[x] = "N"

        # mark any bases definitively called as whatever they are called as
        for item in ["A", "C", "T", "G", "N"]:
            for x in compressed_sequence[item]:  # these are different from reference;
                seq[x] = item

        # add any mixed bases
        for x in compressed_sequence[
            "M"
        ].keys():  # the 'M' option includes iupac coding
            seq[x] = compressed_sequence["M"][x]
        return "".join(seq)

    def compress(self, sequence):
        """reads a string sequence and extracts position - genome information from it.
        returns a dictionary consisting of zero-indexed positions of non-reference bases.

        """
        if not len(sequence) == len(self.reference):
            raise TypeError(
                "sequence must of the same length as reference; seq is {0} and ref is {1}".format(
                    len(sequence), len(self.reference)
                )
            )
        if len(self.reference) == 0:
            raise TypeError("reference cannot be of zero length")

        # we consider - characters to be the same as N
        sequence = sequence.replace("-", "N")

        # we only record differences relative to to refSeq.
        # anything the same as the refSeq is not recorded.
        # a dictionary, M, records the mixed base calls.
        # we also store a set containing either N or M.  This is wasteful of RAM, and the need for it could be removed but
        # it is included here as it massively increases computational speed; the process of combining sets of Ns and Ms is very expensive.

        diffDict = {
            "A": set([]),
            "C": set([]),
            "T": set([]),
            "G": set([]),
            "N": set([]),
            "M": {},
            "U": set([]),
        }

        for i in self.included:  # for the bases we need to compress
            if not sequence[i] == self.reference[i]:  # if it's not reference
                if sequence[i] in ["A", "C", "T", "G", "N"]:
                    diffDict[sequence[i]].add(i)  # if it's a definitively called base
                else:
                    # we regard it as a code representing a mixed base.  we store the results in a dictionary
                    diffDict["M"][i] = sequence[i]
                    diffDict["U"].add(i)
            if sequence[i] == "N":
                diffDict["U"].add(i)

        # check how many Ns or Ms
        if len(diffDict["U"]) > self.maxNs:
            # we store it, but not with sequence details if is invalid
            diffDict = {"invalid": 1}
        else:
            diffDict["invalid"] = 0

        return diffDict

    def _setStats(self, i1, i2):
        """compares either:
        * two sets (if i1 or i2 is a set)
        OR
        * the keys of the dictionaries i1 or i2 is a dictionary

        returns
        * the number of elements in i1
        * the number of elements in i2
        * the number of elements in the union of i1 and i2
        * i1
        * i2
        * the union of set1 and set2

        """
        if isinstance(i1, set):
            retVal1 = i1
        elif isinstance(i1, dict):
            retVal1 = set(i1.keys())
        if isinstance(i2, set):
            retVal2 = i2
        elif isinstance(i2, dict):
            retVal2 = set(i2.keys())

        retVal = retVal2 | retVal1

        return (len(retVal1), len(retVal2), len(retVal), retVal1, retVal2, retVal)

    def countDifferences_byKey(self, keyPair, cutoff=None):
        """compares the in memory refCompressed sequences at
        self.seqProfile[key1] and self.seqProfile[key2]

        Returns the number of SNPs between self._seq1 and self._seq2, and,
        if the pairwise SNP distance is less than cutoff,
        the number of N and Ms in the two sequences and the union of their positions.
        """

        if not type(keyPair) is tuple:
            raise TypeError(
                "Wanted tuple keyPair, but got keyPair={0} with type {1}".format(
                    keyPair, type(keyPair)
                )
            )
        if not len(keyPair) == 2:
            raise TypeError(
                "Wanted a keyPair with two elements, but got {0}".format(keyPair)
            )

        ## test the keys exist
        (key1, key2) = keyPair
        if key1 not in self.seqProfile.keys():
            raise KeyError(
                "Key1={0} does not exist in the in-memory store.".format(key1)
            )
        if key2 not in self.seqProfile.keys():
            raise KeyError(
                "Key1={0} does not exist in the in-memory store.".format(key1)
            )

        # if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling

        ## do the computation
        nDiff = self.countDifferences(
            self.seqProfile[key1], self.seqProfile[key2], cutoff=cutoff
        )

        if nDiff is None:
            return (key1, key2, nDiff, None, None, None, None, None, None)
        elif nDiff <= cutoff:
            seq1_uncertain = self.seqProfile[key1]["N"]
            seq2_uncertain = self.seqProfile[key2]["N"]
            (n1, n2, nboth, N1pos, N2pos, Nbothpos) = self._setStats(
                seq1_uncertain, seq2_uncertain
            )
            return (key1, key2, nDiff, n1, n2, nboth, N1pos, N2pos, Nbothpos)
        else:
            return (key1, key2, nDiff, None, None, None, None, None, None)

    def countDifferences(self, seq1, seq2, cutoff=None):
        """compares seq1 with seq2.

        Ns and Ms (uncertain bases) are ignored in snp computations.

        """
        #  if cutoff is not specified, we use snpCeiling
        if cutoff is None:
            cutoff = self.snpCeiling

        nDiff = 0
        if seq1 is None or seq2 is None:
            return None

        if seq1["invalid"] == 1 or seq2["invalid"] == 1:
            return None

        # compute positions which differ;
        differing_positions = set()
        for nucleotide in ["C", "G", "A", "T"]:

            # we do not consider differences relative to the reference if the other nucleotide is an N or M
            nonN_seq1 = seq1[nucleotide] - seq2["U"]
            nonN_seq2 = seq2[nucleotide] - seq1["U"]
            differing_positions = differing_positions | (nonN_seq1 ^ nonN_seq2)

            # if the number of differences is already larger than the cutoff,
            # then we do not need to perform additional computations; we return None,
            # which is indicated there is no link less than or equal to cutoff
            # unless we are told to return exact differences
            if (
                self.return_none_if_high_snp_distance
                and len(differing_positions) > cutoff
            ):
                return None

        nDiff = len(differing_positions)

        if self.return_none_if_high_snp_distance and nDiff > cutoff:
            return None
        else:
            return nDiff

    def compare(self, key1, key2):
        """ Provides an exact distance between two guids, guid1 and guid2
        
        Will raise KeyError if either guid does not exist. """
       
        return self.countDifferences(
            self.seqProfile[key1], 
            self.seqProfile[key2]
        )

    def compressed_sequence_hash(self, compressed_sequence):
        """returns a string containing a hash of a compressed object.
        Used for identifying compressed objects, including consensus sequences.
        """
        keys = sorted(compressed_sequence.keys())
        serialised_compressed_sequence = ""
        for key in keys:
            if isinstance(compressed_sequence[key], set):
                list_compressed = sorted(list(compressed_sequence[key]))
            else:
                list_compressed = compressed_sequence[key]
            serialised_compressed_sequence = (
                serialised_compressed_sequence + key + ":" + str(list_compressed) + ";"
            )
        h = hashlib.md5()
        h.update(serialised_compressed_sequence.encode("utf-8"))
        md5 = h.hexdigest()
        return md5

    def estimate_expected_proportion(self, seqs):
        """computes the median Ns for seqs, a list.
        Returns None if
        * the length of the sequences in seqs are 0, or
        * there are <= 3 seqs
        """
        if len(seqs) < 3:
            return None
        if len(seqs[0]) == 0:
            return None
        Ns = []
        for seq in seqs:
            Ns.append(seq.count("N"))
        return np.median(Ns) / len(seqs[0])

    def estimate_expected_unk(self, sample_size=30, exclude_guids=set(), unk_type="N"):
        """computes the median unk_type for sample_size guids, randomly selected from all guids except for exclude_guids.
        Used to estimate the expected number of unk_type bases in an alignment.
        unk_type can be one of 'N' 'M' 'N_or_M'.
        """

        guids = list(set(self.seqProfile.keys()) - set(exclude_guids))
        np.random.shuffle(list(guids))

        unks = []
        for guid in guids:
            this_unk = 0
            try:

                if unk_type in ["N", "N_or_M"]:
                    this_unk = this_unk + len(self.seqProfile[guid]["N"])
                if unk_type in ["M", "N_or_M"]:
                    this_unk = this_unk + len(self.seqProfile[guid]["M"])
                unks.append(this_unk)
            except KeyError:
                # it is invalid
                pass
            if len(unks) >= sample_size:
                break
        if len(unks) >= sample_size:
            return np.median(unks)
        else:
            return None

    def estimate_expected_unk_sites(
        self, sample_size=30, sites=set(), exclude_guids=set(), unk_type="N"
    ):
        """computes the median unk_type for sample_size guids, randomly selected from all guids except for exclude_guids.
        Only reports unk_type bases at sites().
        Used to estimate the expected number of unk_type (one of 'N','M',or 'N_or_M' in an alignment"""

        guids = list(set(self.seqProfile.keys()) - set(exclude_guids))
        np.random.shuffle(list(guids))

        unks = []
        for guid in guids:
            try:
                this_unk = 0
                if unk_type in ["N", "N_or_M"]:
                    this_unk = this_unk + len(
                        self.seqProfile[guid]["N"].intersection(sites)
                    )
                if unk_type in ["M", "N_or_M"]:
                    this_unk = this_unk + len(
                        set(self.seqProfile[guid]["M"].keys()).intersection(sites)
                    )
                unks.append(this_unk)
            except ValueError:
                # it is invalid
                pass
            if len(unks) >= sample_size:
                break
        if len(unks) >= sample_size:
            return np.median(unks)
        else:
            return None

    def multi_sequence_alignment(
        self,
        guids,
        output="dict",
        sample_size=30,
        expected_p1=None,
        uncertain_base_type="N",
    ):
        """computes a multiple sequence alignment containing only sites which vary between guids.

        sample_size is the number of samples to randomly sample to estimate the expected number of Ns in
        the population of sequences currently in the server.  From this, the routine computes expected_p1,
        which is expected_expected_N/ the length of sequence.
        if expected_p1 is supplied, then such sampling does not occur.

        output can be either
        'dict', in which case the output is presented as dictionaries mapping guid to results; or
        'df' in which case the results is a pandas data frame like the below, where the index consists of the
        guids identifying the sequences, or

            (index)      aligned_seq  allN  alignN   p_value
            AAACGN-1        AAAC     1       0  0.250000
            CCCCGN-2        CCCC     1       0  0.250000
            TTTCGN-3        TTTC     1       0  0.250000
            GGGGGN-4        GGGG     1       0  0.250000
            NNNCGN-5        NNNC     4       3  0.003906
            ACTCGN-6        ACTC     1       0  0.250000
            TCTNGN-7        TCTN     2       1  0.062500
            AAACGN-8        AAAC     1       0  0.250000

        'df_dict'.  This is a serialisation of the above, which correctly json serialised.  It can be turned back into a
        pandas DataFrame as follows:

        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')     # make the dictionary, see unit test _47
        df = pd.DataFrame.from_dict(res,orient='index')                         # turn it back.

        The p values reported are derived from exact, one-sided binomial tests as implemented in python's scipy.stats.binom_test().

        TEST 1:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences.

        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by
        i) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the genome.  The estimate_expected_unk() function performs this.
        ii) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the relevant  genome.  The estimate_expected_unk() function performs this.

        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.

        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.

        TEST 2:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences
        *at the bases examined in the alignment*.
        This might be relevant if these particular bases are generally hard to call.

        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base at the relevant sites.  The estimate_expected_unk_sites() function performs this.

        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.

        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.

        TEST 3: tests whether the proportion of Ns in the alignment is greater
        than in the bases not in the alignment, for this sequence.

        ## TEST 4:
        tests whether the proportion of Ns in the alignment  for this sequence
        is greater than the median proportion of Ns in the alignment for all other sequences.
        This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
        aligned sequences differs markedly by cluster.

        This test is computed if there are four or more samples in the cluster, and if the alignment length is non-zero.

        """

        # -1 validate input
        if expected_p1 is not None:
            if expected_p1 < 0 or expected_p1 > 1:
                raise ValueError("Expected_p1 must lie between 0 and 1")
        if sample_size is None:
            sample_size = 30

        # step 0: find all valid guids

        valid_guids = []
        invalid_guids = []

        comparatorSeq = {}
        for guid in guids:
            try:
                comparatorSeq[guid] = self.seqProfile[guid]
                valid_guids.append(guid)
            except ValueError:
                invalid_guids.append(guid)
        # Estimate expected N as median(observed Ns),
        # which is a valid thing to do if the proportion of mixed samples is low.

        if expected_p1 is None:
            expected_N1 = self.estimate_expected_unk(
                sample_size=sample_size, exclude_guids=invalid_guids
            )
            if expected_N1 is None:
                expected_p1 = None
            else:
                expected_p1 = expected_N1 / len(self.reference)
        else:
            expected_N1 = np.floor(expected_p1 * len(self.reference))

        return self._msa(
            valid_guids,
            invalid_guids,
            expected_p1,
            output,
            sample_size,
            uncertain_base_type,
        )

    def _msa(
        self,
        valid_guids,
        invalid_guids,
        expected_p1,
        output,
        sample_size,
        uncertain_base_type="N",
    ):
        """perform multisequence alignment and significance tests.

        Parameters:
        valid_guids:  the guids in valid_guids are those on which the msa is computed.

        Additionally, the software performs four significance tests
        on the number of uncertain_bases in an alignment containing only the variant bases between variant_guids.

        If uncertain_base_type == 'N' analyses Ns.
        If uncertain_base_type == 'M' analyses mixed bases.
        If uncertain_base_type == 'N_or_M' analyses both combined

        expected_p1:  The expected proportion of Ns or Ms (as specified by uncertain_base_type) is expected_p1.

        invalid_guids: these are not analysed when computed expected numbers of uncertain_base_type.

        output:   can be either
        'dict', in which case the output is presented as dictionaries mapping guid to results; or
        'df' in which case the results is a pandas data frame like the below, where the index consists of the
        guids identifying the sequences, or

            (index)      aligned_seq  allN  alignN   p_value
            AAACGN-1        AAAC     1       0  0.250000
            CCCCGN-2        CCCC     1       0  0.250000
            TTTCGN-3        TTTC     1       0  0.250000
            GGGGGN-4        GGGG     1       0  0.250000
            NNNCGN-5        NNNC     4       3  0.003906
            ACTCGN-6        ACTC     1       0  0.250000
            TCTNGN-7        TCTN     2       1  0.062500
            AAACGN-8        AAAC     1       0  0.250000

        'df_dict'.  This is a serialisation of the above, which correctly json serialised.  It can be turned back into a
        pandas DataFrame as follows:

        res= sc.multi_sequence_alignment(guid_names[0:8], output='df_dict')     # make the dictionary, see unit test _47
        df = pd.DataFrame.from_dict(res,orient='index')                         # turn it back.

        The p values reported are derived from exact, one-sided binomial tests as implemented in pythons scipy.stats.binom_test().

        TEST 1:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences.

        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by
        i) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the genome.  The estimate_expected_unk() function performs this.
        ii) randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base across the relevant  genome.  The estimate_expected_unk() function performs this.

        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.

        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.

        TEST 2:
        This tests the hypothesis that the number of Ns in the *alignment*
        is GREATER than those expected from the expected_N in the population of whole sequences
        *at the bases examined in the alignment*.
        This might be relevant if these particular bases are generally hard to call.

        Does so by comparing the observed number of Ns in the alignment (alignN),
        given the alignment length (4 in the above case) and an expectation of the proportion of bases which will be N.
        The expected number of Ns is estimated by randomly sampling sample_size guids from those stored in the server and
        observing the number of Ns per base at the relevant sites.  The estimate_expected_unk_sites() function performs this.

        This approach determines the median number of Ns in valid sequences, which (if bad samples with large Ns are rare)
        is a relatively unbiased estimate of the median number of Ns in the good quality samples.

        If there  are not enough samples in the server to obtain an estimate, p_value is not computed, being
        reported as None.

        TEST 3: tests whether the proportion of Ns in the alignment is greater
        than in the bases not in the alignment, for this sequence.

        TEST 4: tests whether the proportion of Ns in the alignment  for this sequence
        is greater than the median proportion of Ns in the alignment for all other sequences.
        This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
        aligned sequences differs markedly by cluster.

        This test is computed if there are two or more samples in the cluster.

        """

        # step 1: define all positions in these sequences for which there is a non-reference position
        # we ignore Ms and Ns in this analysis apart from
        # computing the total number of Ms or Ns in each guid
        nrps = {}
        guid2all = {"N": {}, "M": {}, "N_or_M": {}}
        guid2align = {"N": {}, "M": {}, "N_or_M": {}}

        for guid in valid_guids:
            if self.seqProfile[guid]["invalid"] == 1:
                raise TypeError(
                    "Invalid sequence {0} passed in valid_guids".format(guid)
                )

            for unk_type in ["N", "M"]:

                unks = len(self.seqProfile[guid][unk_type])
                guid2all[unk_type][guid] = unks

            guid2all["N_or_M"][guid] = guid2all["N"][guid] + guid2all["M"][guid]

            for base in ["A", "C", "T", "G"]:
                positions = self.seqProfile[guid][base]
                for position in positions:
                    if (
                        position not in nrps.keys()
                    ):  # if it's non-reference, and we've got no record of this position
                        nrps[
                            position
                        ] = set()  # then we generate a set of bases at this position
                    nrps[position].add(
                        base
                    )  # either way add the current non-reference base there

        # step 2: for the non-reference called positions, check if there's a reference base there.
        for guid in valid_guids:
            for position in nrps.keys():
                psn_accounted_for = 0
                for base in ["A", "C", "T", "G"]:
                    if position in self.seqProfile[guid][base]:
                        psn_accounted_for = 1
                if psn_accounted_for == 0:
                    # it is reference; this guid has no record of a variant base at this position, so it must be reference.
                    nrps[position].add(self.reference[position])

        # step 3: find those which have multiple bases at a position
        variant_positions = set()
        for position in nrps.keys():
            if len(nrps[position]) > 1:
                variant_positions.add(position)

        # step 3b: define constant base frequencies - i.e. frequencies outside the region of variation.  This is used by iqtree, see http://www.iqtree.org/doc/Command-Reference
        constant_positions = list(self.included - variant_positions)
        constant_nt = [self.reference_list[x] for x in constant_positions]
        fconst_dict = Counter(constant_nt)
        fconst = [0, 0, 0, 0]
        for i, base in enumerate(["A", "C", "G", "T"]):
            if base in fconst_dict.keys():
                fconst[i] = fconst_dict[base]

        # step 4: determine the sequences of all variant positions.
        ordered_variant_positions = sorted(list(variant_positions))
        guid2seq = {}
        guid2msa_seq = {}
        for guid in valid_guids:
            guid2seq[guid] = []
            for position in ordered_variant_positions:
                this_base = self.reference[position]
                for base in ["A", "C", "T", "G", "N", "M"]:
                    if not base == "M":
                        positions = self.seqProfile[guid][base]
                    else:
                        positions = self.seqProfile[guid][base].keys()

                    if position in positions:
                        this_base = base
                guid2seq[guid].append(this_base)
            guid2msa_seq[guid] = "".join(guid2seq[guid])

        # step 5: determine the expected_p2 at the ordered_variant_positions:
        expected_N2 = self.estimate_expected_unk_sites(
            sample_size=sample_size,
            exclude_guids=invalid_guids,
            sites=set(ordered_variant_positions),
            unk_type=uncertain_base_type,
        )
        if expected_N2 is None:
            expected_p2 = None
        elif len(ordered_variant_positions) == 0:
            expected_p2 = None
        else:
            expected_p2 = expected_N2 / len(ordered_variant_positions)

        # step 6: perform Binomial tests on all samples
        if len(valid_guids) > 0:
            guid2pvalue1 = {}
            guid2pvalue2 = {}
            guid2pvalue3 = {}
            guid2pvalue4 = {}
            guid2observed_p = {}
            guid2expected_p1 = {}
            guid2expected_p2 = {}
            guid2expected_p3 = {}
            guid2expected_p4 = {}
            for guid in valid_guids:

                # compute p value 1.  This tests the hypothesis that the number of Ns in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences.
                for unk_type in ["N", "M"]:
                    guid2align[unk_type][guid] = guid2msa_seq[guid].count(unk_type)
                guid2align["N_or_M"][guid] = (
                    guid2align["N"][guid] + guid2align["M"][guid]
                )

                if (
                    expected_p1 is None
                ):  # we don't have an expectation, so we can't assess the first binomial test;
                    p_value1 = None
                    observed_p = None
                elif (
                    len(guid2msa_seq[guid]) == 0
                ):  # we don't have any information to work with
                    p_value1 = None
                    observed_p = None
                else:
                    observed_p = guid2align[uncertain_base_type][guid] / len(
                        guid2msa_seq[guid]
                    )
                    p_value1 = binom_test(
                        guid2align[uncertain_base_type][guid],
                        len(guid2msa_seq[guid]),
                        expected_p1,
                        alternative="greater",
                    )
                    # print(uncertain_base_type, guid2msa_seq[guid],guid2align[uncertain_base_type][guid],len(guid2msa_seq[guid]), observed_p, expected_p1, p_value1)

                guid2pvalue1[guid] = p_value1
                guid2observed_p[guid] = observed_p
                guid2expected_p1[guid] = expected_p1

                # compute p value 2.  This tests the hypothesis that the number of Ns in the *alignment*
                # is GREATER than those expected from the expected_N in the population of whole sequences
                # at these sites.
                if (
                    expected_p2 is None
                ):  # we don't have an expectation, so we can't assess the binomial test;
                    p_value2 = None
                elif (
                    len(guid2msa_seq[guid]) == 0
                ):  # we don't have any information to work with
                    p_value2 = None
                else:
                    p_value2 = binom_test(
                        guid2align[uncertain_base_type][guid],
                        len(guid2msa_seq[guid]),
                        expected_p2,
                        alternative="greater",
                    )
                guid2pvalue2[guid] = p_value2
                guid2expected_p2[guid] = expected_p2

                # compute p value 3.  This tests the hypothesis that the number of Ns in the alignment of THIS SEQUENCE
                # is GREATER than the number of Ns not in the alignment  IN THIS SEQUENCE
                # based on sequences not in the alignment

                expected_p3 = (
                    guid2all[uncertain_base_type][guid]
                    - guid2align[uncertain_base_type][guid]
                ) / (len(self.reference) - len(guid2msa_seq[guid]))
                p_value = binom_test(
                    guid2align[uncertain_base_type][guid],
                    len(guid2msa_seq[guid]),
                    expected_p3,
                    alternative="greater",
                )
                guid2pvalue3[guid] = p_value
                guid2expected_p3[guid] = expected_p3

                # compute p value 4.
                # tests whether the proportion of Ns in the alignment  for this sequence
                # is greater than the median proportion of Ns in the alignment for all other sequences.
                # This is a sensible option if the test is performed per-cluster, and the proportion of Ns in the
                # aligned sequences differs markedly by cluster.
                # This test is computed if there are two or more samples in the cluster.

                # get the sequences which are not guid;
                other_seqs = []
                for this_guid in guid2msa_seq.keys():
                    if not guid == this_guid:
                        other_seqs.append(guid2msa_seq[this_guid])
                expected_p4 = self.estimate_expected_proportion(other_seqs)
                if expected_p4 is None:
                    p_value = None
                else:
                    p_value = binom_test(
                        guid2align[uncertain_base_type][guid],
                        len(guid2msa_seq[guid]),
                        expected_p4,
                        alternative="greater",
                    )
                guid2pvalue4[guid] = p_value
                guid2expected_p4[guid] = expected_p4

            # assemble dataframe
            df1 = pd.DataFrame.from_dict(guid2msa_seq, orient="index")
            df1.columns = ["aligned_seq"]
            df1["aligned_seq_len"] = len(ordered_variant_positions)

            df2n = pd.DataFrame.from_dict(guid2all["N"], orient="index")
            df2n.columns = ["allN"]
            df3n = pd.DataFrame.from_dict(guid2align["N"], orient="index")
            df3n.columns = ["alignN"]

            df2m = pd.DataFrame.from_dict(guid2all["M"], orient="index")
            df2m.columns = ["allM"]
            df3m = pd.DataFrame.from_dict(guid2align["M"], orient="index")
            df3m.columns = ["alignM"]

            df2mn = pd.DataFrame.from_dict(guid2all["N_or_M"], orient="index")
            df2mn.columns = ["allN_or_M"]
            df3mn = pd.DataFrame.from_dict(guid2align["N_or_M"], orient="index")
            df3mn.columns = ["alignN_or_M"]

            df4 = pd.DataFrame.from_dict(guid2pvalue1, orient="index")
            df4.columns = ["p_value1"]
            df5 = pd.DataFrame.from_dict(guid2pvalue2, orient="index")
            df5.columns = ["p_value2"]
            df6 = pd.DataFrame.from_dict(guid2pvalue3, orient="index")
            df6.columns = ["p_value3"]
            df7 = pd.DataFrame.from_dict(guid2pvalue4, orient="index")
            df7.columns = ["p_value4"]
            df8 = pd.DataFrame.from_dict(guid2observed_p, orient="index")
            df8.columns = ["observed_proportion"]
            df9 = pd.DataFrame.from_dict(guid2expected_p1, orient="index")
            df9.columns = ["expected_proportion1"]
            df10 = pd.DataFrame.from_dict(guid2expected_p3, orient="index")
            df10.columns = ["expected_proportion2"]
            df11 = pd.DataFrame.from_dict(guid2expected_p3, orient="index")
            df11.columns = ["expected_proportion3"]
            df12 = pd.DataFrame.from_dict(guid2expected_p4, orient="index")
            df12.columns = ["expected_proportion4"]
            df12["what_tested"] = uncertain_base_type

            df = df1.merge(df2n, left_index=True, right_index=True)
            df = df.merge(df3n, left_index=True, right_index=True)
            df = df.merge(df2m, left_index=True, right_index=True)
            df = df.merge(df3m, left_index=True, right_index=True)
            df = df.merge(df2mn, left_index=True, right_index=True)
            df = df.merge(df3mn, left_index=True, right_index=True)
            df = df.merge(df4, left_index=True, right_index=True)
            df = df.merge(df5, left_index=True, right_index=True)
            df = df.merge(df6, left_index=True, right_index=True)
            df = df.merge(df7, left_index=True, right_index=True)
            df = df.merge(df8, left_index=True, right_index=True)
            df = df.merge(df9, left_index=True, right_index=True)
            df = df.merge(df10, left_index=True, right_index=True)
            df = df.merge(df11, left_index=True, right_index=True)
            df = df.merge(df12, left_index=True, right_index=True)

            retDict = {
                "variant_positions": ordered_variant_positions,
                "invalid_guids": invalid_guids,
                "valid_guids": valid_guids,
                "expected_p1": expected_p1,
                "sample_size": sample_size,
                "df_dict": df.to_dict(orient="index"),
                "what_tested": uncertain_base_type,
                "outgroup": "",
                "creation_time": datetime.datetime.now().isoformat(),
                "fconst": fconst_dict,
            }

            return MSAResult(**retDict)

        else:
            return None

    def raise_error(self, token):
        """raises a ZeroDivisionError, with token as the message.
        useful for unit tests of error logging"""
        raise ZeroDivisionError(token)