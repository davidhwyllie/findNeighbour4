""" utilities to support local storage of sequence compressed and fasta data """

import io
import os
import time
import tarfile
import lzma
import gzip
import pickle
from findn.seq2json import SeqDictConverter


class LocalStore:
    """methods to store json sequence objects, identified by a guid, in a
    tar file.

    Transparent compression of the json object occurs.
    """

    def __init__(
        self, tarfile_name, write_batch_size=1000, compression_method="pickle"
    ):
        """creates a LocalStore object, providing access to a tarfile

        Parameters:
        tarfile_name:
                    the tarfile to use.  Tarfile will be created if it does
                    not exist
                    
        write_batch_size:
                    how many sequences are held in ram before writing to tar file
                    (higher values speed writes but use more RAM)

        compression_method:
                    the compression method used.
                    options are: 'lzma' 'gzip' 'pickle'

                    benchmarking indicates the fastest method is pickling but that lzma makes much smaller files.

                    If the potential security risk associated with unpickling (if the tar file in manipulated,
                    unpickling can be a route to remote code execution) the pickle is fastest
                    Otherwise gzip is a good compromise

                    Using a synthetic 5,000 sample set including a consecutive series of bases
                    (note: this is found in some types of sequences, typically with runs of Ns)
                    which may favour lzma more than real samples.

                    5000 synthetic compressed sequences

                          ADD       READ      FILE SIZE
                    lzma  82.448591 11.517303  15923200
                    gzip  26.173684  9.019711 115210240
                    pickle 1.735323  3.132463 153610240


        Returns:
                    Nothing
        """
        self.sjc = SeqDictConverter()
        self.tarfile_name = tarfile_name
        self._to_store = []
        self.write_batch_size = write_batch_size
        self.compression_method = compression_method
        # create and close tarfile if it does not exist
        try:
            with tarfile.open(self.tarfile_name, "x"):
                pass  # just create the file
        except FileExistsError:
            pass  # no action

    def __del__(self):
        """destructor"""
        self._write()  # ensure everything is written

    def sequence_ids(self):
        """returns a list of sequence_ids stored in the tarfile"""
        with tarfile.TarFile(self.tarfile_name, "r") as tar:
            sequence_ids = tar.getnames()
        return sequence_ids

    def store(self, sequence_id, obj):
        """stores a compressed version of the sequence object obj
        identified by sequence_id in a tar file

        Parameters:
            obj: a reference compressed object

        Returns:
            True

        Note:
            1. sequence_id will be part of a filename, so alphanumeric, _, - only should be used
            2. writing to the tar file may be delayed;
            to increase write speed,s batches of files are written;
            this method stores obj for writing. Writing should occur transparently
            without further action."""
        self._to_store.append((sequence_id, self._compress(obj)))

        if len(self._to_store) > self.write_batch_size:
            self._write()

    def flush(self):
        """writes any stored items to the tarfile"""
        self._write()

    def read(self, sequence_id):
        """reads the object stored as sequence_id
        returns a tuple
            (sequence_id, reference compressed object)"""
        with tarfile.TarFile(self.tarfile_name, "r") as tar:
            buf = tar.extractfile(sequence_id)
            compressed_bytes = buf.read()
            return (sequence_id, self._decompress(compressed_bytes))

    def read_many(self, select_sequence_ids=None):
        """reads items in the tar file

        Parameters:
        select_sequence_ids: either a set of sequence_ids to find, or None (in which case all samples are loaded)

        Returns:
        a generator which provides tuples
            (sequence_id, reference compressed object)
            
        Note:
        if select_sequence_ids is None, will read all samples """

        all_sequence_ids = self.sequence_ids()

        # find the ones we need, in the order they are in the tar file
        if select_sequence_ids is not None:
            select_sequence_ids = select_sequence_ids.intersection(set(all_sequence_ids))
            sequence_ids = [x for x in all_sequence_ids if x in select_sequence_ids]
        else:
            sequence_ids = all_sequence_ids

        # iterate over the tarfile, returning all objects
        with tarfile.TarFile(self.tarfile_name, "r") as tar:
            for sequence_id in sequence_ids:
                buf = tar.extractfile(sequence_id)
                compressed_bytes = buf.read()
                yield (sequence_id, self._decompress(compressed_bytes))

    def _write(self):
        """writes the contents of self._to_store, a list of reference compressed sequences,
        to the tarfile at self.tarfile_name"""
        with tarfile.TarFile(self.tarfile_name, "a") as tar:
            for sequence_id, compressed_bytes in self._to_store:
                info = tarfile.TarInfo(
                    name=sequence_id
                )
                info.size = len(compressed_bytes)
                uid = os.getuid()
                gid = os.getgid()
                info.uid = uid
                info.gid = gid
                info.mtime = int(time.time())

                tar.addfile(info, io.BytesIO(compressed_bytes))
        self._to_store = []

    def _compress(self, obj):
        """converts a python object obj to json, and compresses it
        returning a bytes object

        Parameters:
        obj: a python object, serialisable to json

        Returns:
        bytes object, a json serialised, LZMA compressed string"""

        if self.compression_method == "pickle":
            to_compress = pickle.dumps(obj)
            return to_compress
        elif self.compression_method == "lzma":
            to_compress = self.sjc.to_json(obj)
            return lzma.compress(to_compress.encode("utf-8"))
        elif self.compression_method == "gzip":
            to_compress = self.sjc.to_json(obj)
            return gzip.compress(to_compress.encode("utf-8"))
        else:
            raise ValueError(
                "{0} is not a valid compression_method, options are lzma, gzip and pickle"
            )

    def _decompress(self, compressedobj):
        """converts compressedobj, a bytes object containing reference compressed sequence data,
        generated by _compress,
        back to the original python object

        Parameters:

        compressedobj: a bytes object from _compress()

        Returns:
        a python object"""

        if self.compression_method == "pickle":
            return pickle.loads(compressedobj)
        elif self.compression_method == "lzma":
            to_decompress = lzma.decompress(compressedobj)
            return self.sjc.from_json(to_decompress.decode("utf-8"))
        elif self.compression_method == "gzip":
            to_decompress = gzip.decompress(compressedobj)
            return self.sjc.from_json(to_decompress.decode("utf-8"))
        else:
            raise ValueError(
                "{0} is not a valid compression_method, options are lzma, gzip and pickle"
            )
