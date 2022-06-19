#!/usr/bin/env python
""" 
Stores reference compressed sequence data, in a sparse-array format suitable for 
rapid construction of on-hot distance matrices, in a local directory.

A component of the findNeighbour4 system for bacterial relatedness monitoring
Copyright (C) 2021 David Wyllie david.wyllie@ukhsa.gov.uk
repo: https://github.com/davidhwyllie/findNeighbour4

This program is free software: you can redistribute it and/or modify
it under the terms of the MIT License as published
by the Free Software Foundation.  See <https://opensource.org/licenses/MIT>, and the LICENSE file.
bu
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without tcen the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  S
"""

# import libraries
import os
import logging
import math
import numpy as np
from scipy import sparse as sp
import progressbar
from localstore.localstoreutils import LocalStore


class PCASequenceStore:
    """stores sequences in a folder, in a format optimised for fast
    matrix production for PCA

    Supports the following activities:
    * Storing samples, either from strings containing sequence or from reference compressed representations
    * Recovering samples
    """

    def __init__(self, reference_sequence, persistdir, show_bar=True):
        """initiates a PCASequenceStore object

        Parameters:
        reference_sequence: a string containing the reference sequence.
        persistdir: a directory to store the output in.  If the directory does not exist, then the directory is created
        show_bar: display a progressbar during long operations
        """

        self.show_bar = show_bar
        # create an output directory if it does not exist
        os.makedirs(persistdir, exist_ok=True)

        # note what it is
        self.persistdir = persistdir

        # sanity check
        if not isinstance(reference_sequence, str):
            raise TypeError("Reference sequence must be a string")

        # ensure reference is upper case
        self.reference_sequence = reference_sequence.upper()
        self.reference_sequence_list = list(self.reference_sequence)
        self.reference_sequence_length = len(self.reference_sequence)

        # compute offsets for mapping sequence data to integers
        self.offsets = dict(
            A=0,
            C=self.reference_sequence_length,
            T=self.reference_sequence_length * 2,
            G=self.reference_sequence_length * 3,
        )
        self.rcs_array_length = 4 * self.reference_sequence_length
        self.N_array_length = self.reference_sequence_length

        # create a localstore object to keep the data in
        self.tarfile_name = os.path.join(self.persistdir, "sparse.tar")
        self.localstore = LocalStore(self.tarfile_name)
        self.sequence_ids = set(self.localstore.sequence_ids())
        logging.info(
            "Sparse array storage file {0} contains {1} samples".format(
                self.tarfile_name, len(self.sequence_ids)
            )
        )

        # to generate a table containing the mapping of variant_int_id to variant positions,
        # call _generate_mapping_file()
        self.mapping = None  # none created

    def _variant_int_id(self, nt, pos):
        """computes a variant_int_id from a base / nt combination

        Paramters:
        nt: a nucleotide, either A,C,G,T
        pos: the position the variant occurs

        """

        nt = nt.upper()
        return self.offsets[nt] + pos

    def allelemodel2vmodel(self, allelemodel):
        """converts an allele model to a variant model array

        An allelemodel is a (possibly sparse) array of length 4x genome length.  It contains all position/variant combinations.
        By contrast, a variant model is a (possibly sparse) array of length 1x genome length.
        It contains the positions in the genome of variation.
        """
        to_slice = allelemodel[0]
        return np.sum(
            [
                to_slice[
                    self.offsets["A"]: (
                        self.offsets["A"] + self.reference_sequence_length
                    )
                ],
                to_slice[
                    self.offsets["C"]: (
                        self.offsets["C"] + self.reference_sequence_length
                    )
                ],
                to_slice[
                    self.offsets["G"]: (
                        self.offsets["G"] + self.reference_sequence_length
                    )
                ],
                to_slice[
                    self.offsets["T"]: (
                        self.offsets["T"] + self.reference_sequence_length
                    )
                ],
            ],
            axis=0,
        )

    def identify_variant(self, variant_int_id):
        """identifies a variant from a variant_int_id

        Parameters: variant_int_id is an internally computed integer representing a variant

        returns:
        a dictionary containing the position and variant with keys
        variant_int_id
        base
        pos

        """

        nt_int = math.floor(variant_int_id / self.reference_sequence_length)
        nt = ["A", "C", "T", "G"][nt_int]
        offset = self.offsets[nt]
        pos = variant_int_id - offset
        return dict(pos=pos, variant_int_id=variant_int_id, nt=nt)

    def _generate_mapping_file(self):
        """generates and saves a numpy array relating the
        variant_int_id to position and allele"""

        mapping_file_name = os.path.join(self.persistdir, "mapping.csv")

        if not os.path.exists(mapping_file_name):
            logging.info("Creating mapping file [one-off operation]")
            mapping_file_records = []

            nt_int = 0
            for nt in ["A", "C", "T", "G"]:
                nt_int += 1
                for pos in list(range(self.reference_sequence_length)):

                    variant_int_id = self._variant_int_id(nt, pos)
                    mapping_file_records.append([variant_int_id, nt_int, pos])
            self.mapping = np.array(mapping_file_records)
            np.save(mapping_file_name, self.mapping)

        else:
            # read the mapping from disc
            self.mapping = np.load(mapping_file_name)

    def _encode_rcs(self, obj):
        """encodes a reference compressed sequence obj obj,
        which is a dictionary comprising a dictionary with
        * zero indexed non reference positions of A,C,G,T
        * zero indexed positions of uncertain bases as 'U'

        Example, in which there is no variation w.r.t the reference:
        obj = {'A': {}, 'C': {}, 'T': {}, 'G': {}, 'N': {}, 'M': {}, 'invalid': 0, 'U': []}

        Parameters:
        obj: a dictionary containing reference compressed sequence data as above

        Returns:
        a dictionary with two keys:
        variants: a scipy.sparse.coo_array with dimensions 1 x self.rcs_array_length
        ns : a scipy.sparse.coo_array with dimensions 1 x self.reference_length
        """

        rcs_array = []
        n_array = []
        for char in ["A", "C", "G", "T"]:
            for pos in obj[char]:
                rcs_array.append(self._variant_int_id(char, pos))
        n_array = obj["U"]

        # note that int8 is sufficient to encode the single bit here, but if we add up these arrays,
        # e.g. when counting over millions of samples, a larger data type is needed
        # int32 allow up to 2 billion samples
        rcs_sparse_array = sp.coo_array(
            ([1] * len(rcs_array), ([0] * len(rcs_array), rcs_array)),
            shape=(1, self.rcs_array_length),
            dtype=np.int32,
        )

        n_sparse_array = sp.coo_array(
            ([1] * len(n_array), ([0] * len(n_array), n_array)),
            shape=(1, self.N_array_length),
            dtype=np.int32,
        )

        return dict(variants=rcs_sparse_array, ns=n_sparse_array)

    def _encode_sequence(self, seq):
        """encodes a sequence seq, which is a string,
        as a dictionary to two sparse matrices:
        one covers non-reference variation, and the other covers Ns (unknown characters)

        A,C,T,G are encoded (irrespective of case)
        All other characters are encoded as N (Unknown)

        Parameters:
        seq: a string, the same length as the reference, to encode

        Returns:
        a dictionary with two keys:
        variants: a scipy.sparse.coo_array with dimensions 1 x self.rcs_array_length
        ns : a scipy.sparse.coo_array with dimensions 1 x self.reference_length
        """

        rcs_array = []
        n_array = []
        seq = seq.upper()

        # sanity check
        if not len(seq) == self.reference_sequence_length:
            raise ValueError(
                "The supplied sequence is the wrong length.  Should be {0} but is {1}".format(
                    self.reference_sequence_length, len(seq)
                )
            )

        for i, char in enumerate(list(seq)):
            if char in ["A", "C", "G", "T"]:
                # encode this if non reference
                if not self.reference_sequence_list[i] == char:
                    rcs_array.append(self._variant_int_id(char, i))
            else:
                # mark it as an N
                n_array.append(i)

        rcs_sparse_array = sp.coo_array(
            ([1] * len(rcs_array), ([0] * len(rcs_array), rcs_array)),
            shape=(1, self.rcs_array_length),
            dtype=np.int32,
        )

        n_sparse_array = sp.coo_array(
            ([1] * len(n_array), ([0] * len(n_array), n_array)),
            shape=(1, self.N_array_length),
            dtype=np.int32,
        )

        return dict(variants=rcs_sparse_array, ns=n_sparse_array)

    def sequence_ids(self):
        """returns a list of sequence_ids stored in the tarfile"""
        return self.localstore.sequence_ids()

    def store_sequence(self, sequence_id, seq):
        """stores a sequence, identified by sequence_id
        whose sequence is a string in seq

        to force the changes to disc, call flush().
        It is more efficient to call flush() when multiple samples have been stored; by default,
        flushing occurs after every 1,000 samples."""
        if sequence_id not in self.sequence_ids:
            obj = self._encode_sequence(seq)
            self._store(sequence_id, obj)

    def store_rcs(self, sequence_id, rcs):
        """store a reference compressed object rcs, identified by sequence_id, in the tar file.
        to force the changes to disc, call flush().
        It is more efficient to call flush() when multiple samples have been stored; by default,
        flushing occurs after every 1,000 samples."""
        if sequence_id not in self.sequence_ids:
            obj = self._encode_rcs(rcs)
            self._store(sequence_id, obj)

    def _store(self, sequence_id, obj):
        """stores an compressed version of the sequence object obj
        identified by sequence_id in a tar file

        Parameters:
            obj: an object for storage

        Returns:
            True

        Note:
            1. sequence_id will be part of a filename, so alphanumeric, _, - only should be used
            2. writing to the tar file may be delayed;
            to increase write speed,s batches of files are written;
            this method stores obj for writing. Writing should occur transparently
            without further action.
            3. However, to force the changes to disc, call flush()"""
        self.localstore.store(sequence_id, obj)
        self.sequence_ids.add(sequence_id)

    def flush(self):
        """ensures data is written to disc"""
        self.localstore.flush()

    def read(self, sequence_id):
        """reads the object stored as sequence_id
        returns a tuple
            (sequence_id, reference compressed object)"""
        return self.localstore.read(sequence_id)

    def read_all(self):
        """reads all items in the tar file
        Returns:
        a generator which provides tuples
            (sequence_id, reference compressed object)
        """

        # iterate over the tarfile, returning all objects
        for sequence_id, obj in self.localstore.read_all():
            if sequence_id is None:
                return
            yield sequence_id, obj

    def read_many(self, select_sequence_ids=None):
        """reads selected items in the tar file

        Parameters:
        select_sequence_ids: either a set of sequence_ids to find, or None (in which case all samples are loaded)

        Returns:
        a generator which provides tuples
            (sequence_id, reference compressed object)

        Note:
        if select_sequence_ids is None, will read all samples.
        read_all is a much faster method of doing this."""

        for sequence_id, obj in self.localstore.read_many(select_sequence_ids):
            if sequence_id is not None:
                yield sequence_id, obj

    def summarise_all(self):
        """summarises variant positions all items in the tar file
        Returns:
        a dictionary, containing two arrays with keys
        allelemodel
        mmodel

        These are 4x and 1x the genome length, respectively and contain the number of
        sequences with variants and the number of sequences with Ns at each position

        """

        # create arrays
        encode_reference = self._encode_sequence(self.reference_sequence)
        allelemodel = encode_reference["variants"].toarray()  # will be empty
        mmodel = encode_reference["ns"].toarray()
        sequence_ids = []

        # iterate over the tarfile, returning all objects
        if self.show_bar:
            bar = progressbar.ProgressBar()

        # iterate over the tarfile, returning all objects
        num_loaded = 0
        for sequence_id, obj in self.localstore.read_all():
            num_loaded += 1
            if sequence_id is None:
                break
            sequence_ids.append(sequence_id)
            allelemodel = allelemodel + obj["variants"]
            mmodel = mmodel + obj["ns"]
            if self.show_bar:
                bar.update(num_loaded)

        if self.show_bar:
            bar.finish()

        return dict(
            allelemodel=allelemodel,
            mmodel=mmodel,
            vmodel=[self.allelemodel2vmodel(allelemodel)],
            sequence_ids=set(sequence_ids),
        )

    def summarise_many(self, select_sequence_ids=None, starting_result=None):

        """summarises variant positions in selected items in the tar file

        Parameters:
        select_sequence_ids: either a set of sequence_ids to find, or None (in which case all samples are loaded)
        starting_result:   if starting_result is supplied, it should the result of a previous call to summarise_all or summarise_many:
        a dictionary with keys as above.

        Returns:
        a dictionary, containing
        allelemodel: array
        mmodel: array
        sequence_ids : set

        These are 4x and 1x the genome length, respectively and contain the number of
        sequences with variants and the number of sequences with Ns at each position

        Note:
        if select_sequence_ids is None, will read all samples.
        read_all is a much faster method of doing this."""

        # create arrays
        if starting_result is not None:
            # allelemodel and mmodel have to be in the keys

            allelemodel = starting_result["allelemodel"]
            mmodel = starting_result["mmodel"]
            sequence_ids = starting_result["sequence_ids"]
        else:
            encode_reference = self._encode_sequence(self.reference_sequence)
            allelemodel = encode_reference["variants"].toarray()  # will be empty
            mmodel = encode_reference["ns"].toarray()
            sequence_ids = set()

        # iterate over the tarfile, returning all objects
        if self.show_bar:
            if select_sequence_ids is not None:
                max_value=len(select_sequence_ids)

            else:
                max_value=len(self.sequence_ids)

            bar = progressbar.ProgressBar(max_value=max_value)

        num_loaded = 0
        for sequence_id, obj in self.localstore.read_many(select_sequence_ids):
            if sequence_id is None:
                break
            num_loaded += 1
            if self.show_bar:
                if num_loaded <= max_value:
                    bar.update(num_loaded)
                else:
                    # this is an error condition, one explanation for which is that the .tar file has been updated
                    logging.warning("Unexpectedly high number of samples encountered {0} vs {1}".format(num_loaded, max_value))

            if sequence_id not in sequence_ids:
                sequence_ids.add(sequence_id)
                allelemodel = allelemodel + obj["variants"]
                mmodel = mmodel + obj["ns"]

        if self.show_bar:
            bar.finish()

        return dict(
            allelemodel=allelemodel,
            mmodel=mmodel,
            vmodel=[self.allelemodel2vmodel(allelemodel)],
            sequence_ids=sequence_ids,
        )

    def update_summary(self, starting_result=None):
        """summarises variant positions in all items in the tar file,
        considering only those which have not previously been considered.

        Parameters:
        starting_result:   if starting_result is supplied, it should the result of a previous call to summarise_all or summarise_many:
        a dictionary with keys as above.  If it not supplied, this is the equivalent of summarise_all()

        Returns:
        a dictionary, containing
        allelemodel: array
        mmodel: array
        sequence_ids : set

        These are 4x and 1x the genome length, respectively and contain the number of
        sequences with variants and the number of sequences with Ns at each position

        Note:
        """

        if starting_result is None:
            return self.summarise_all()

        # starting result is supplied
        assessed_sequence_ids = starting_result["sequence_ids"]
        total_sequence_ids = self.sequence_ids()

        new_sequence_ids = total_sequence_ids - assessed_sequence_ids

        if len(new_sequence_ids) == 0:
            return starting_result
        else:
            return self.summarise_many(
                select_sequence_ids=new_sequence_ids, starting_result=starting_result
            )
