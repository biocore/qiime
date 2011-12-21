#!/usr/bin/env python

from datetime import datetime
from json import dumps
from types import NoneType
from operator import itemgetter, xor
from itertools import izip
from collections import defaultdict, Hashable
from pysparse.spmatrix import LLMatType, ll_mat
from numpy import ndarray, asarray, array, newaxis, zeros
from cogent.util.misc import unzip, flatten
from qiime.format import get_biom_format_version_string,\
        get_biom_format_url_string
from qiime.util import get_qiime_library_version

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2007-2011, QIIME"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.3.0dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Prototype"

"""add in support methods for filtering, etc"""

"""rows and cols might be backwards to represent samples/otus"""

"""got atleast one 1 bug... keys() method of ll_mat is unreliable following a
slice:
    >>> from pysparse.spmatrix import ll_mat
    >>> foo = ll_mat(2,2)
    >>> foo[0,0] = 5
    >>> foo[0,1] = 6
    >>> foo[1,0] = 7
    >>> foo[1,1] = 8
    >>> col = foo[:,0]
    >>> print col
    ll_mat(general, [2,1]):
         5.000000 
          7.000000 

          >>> col.shape
          (2, 1)
          >>> col.items()
          [((0, 0), 5.0), ((1, 0), 7.0)]
          >>> col.keys()
          [[0, 1], [0, 0]]

   ...i don't have internet right now so cannot report but this should be done
   soon
   """

class TableException(Exception):
    pass

def to_ll_mat(values, transpose=False):
    """Tries to returns a populated ll_mat object
    
    NOTE: assumes the max value observed in row and col defines the size of the
    matrix
    """
    ### SCREAMING decomposition here... this is ugly

    # if it is a vector
    if isinstance(values, ndarray) and len(values.shape) == 1:
        if transpose:
            mat = ll_mat(len(values), 1)
        else:
            mat = ll_mat(1, len(values))

        for i,v in enumerate(values):
            if transpose:
                mat[i,0] = v
            else:
                # pysparse whines about numpy.int64
                mat[0,i] = float(v)

        return mat
    # list of dicts, each representing a row in row order
    elif isinstance(values, list) and isinstance(values[0], dict):
        n_rows = len(values)
        all_keys = flatten(map(lambda x:x.keys(),values))
        n_cols = max(all_keys, key=itemgetter(1))[1] # index 1 == column
       
        if transpose:
            mat = ll_mat(n_cols+1, n_rows)
        else:
            mat = ll_mat(n_rows, n_cols+1) # adjust cols for 0-based index

        # we ignore the given row as these _should_ be dicts of vectors
        for row,d in enumerate(values):
            for (ignore, column),v in d.items():
                if transpose:
                    mat[column, row] = v
                else:
                    mat[row, column] = v
        return mat
    # list of ll_mat objects (EXPECTING VECTORS!!)
    # this is unfortunately necessary because ll_mat decided that the
    # 'keys' method should behave like dict.keys. thanks. implement the
    # other dict methods like items() and values() but LIE about keys.
    elif isinstance(values, list) and isinstance(values[0], LLMatType):
        n_rows = len(values)
        all_keys = []

        # god damn you pysparse
        for d in values:
            all_keys.extend([k for k,v in d.items()])
        n_cols = max(all_keys, key=itemgetter(1))[1] # index 1 == column
        
        if transpose:
            mat = ll_mat(n_cols+1, n_rows)
        else:
            mat = ll_mat(n_rows, n_cols+1) # adjust cols for 0-based index

        # we ignore the given row as these _should_ be dicts of vectors
        for row,d in enumerate(values):
            for (ignore, column),v in d.items():
                if transpose:
                    mat[column, row] = v
                else:
                    mat[row, column] = v
        return mat

    else:
        # if values does not appear dict-like, try a cast
        if not hasattr(values, 'items'): 
            try:
                values = dict(values)
            except IndexError:
                raise TableException, "Unable to cast to known type"
            except TypeError:
                raise TableException, "Unable to cast to known type"
        
        rows, cols = unzip([(r,c) for (r,c),v in values.items()])
        n_rows = sorted(rows)[-1]
        n_cols = sorted(cols)[-1]

        if transpose:
            mat = ll_mat(n_cols+1, n_rows+1)
        else:
            mat = ll_mat(n_rows+1, n_cols+1)

        for (row,col),val in values.items():
            if transpose:
                mat[col,row] = val
            else:
                mat[row,col] = val

        return mat
    
class Table(object):
    """ """
    _biom_type = None
    _biom_matrix_type = None

    def __init__(self, Data, SampleIds, ObservationIds, SampleMetadata=None, 
                 ObservationMetadata=None, TableId=None, **kwargs):
        self.TableId = TableId
        self._data = Data

        ### DO WE WANT IMMUTABLE TYPES? or some programitic lie to that effect?
        self.SampleIds = SampleIds
        self.ObservationIds = ObservationIds
        self.SampleMetadata = SampleMetadata
        self.ObservationMetadata = ObservationMetadata

        self._verify_metadata()
        self._cast_metadata()

    def _conv_to_self_type(self, vals, transpose=False):
        """For converting vectors to a compatible self type"""
        raise NotImplementedError

    def _verify_metadata(self):
        """Obtain some notion of sanity on object construction with inputs"""
        try:
            n_obs, n_samp = self._data.shape
        except:
            n_obs = n_samp = 0

        if n_obs != len(self.ObservationIds):
            raise TableException, \
                    "Number of ObservationIds differs from matrix size!"

        if n_obs != len(set(self.ObservationIds)):
            raise TableException, "Duplicate ObservationIds"

        if n_samp != len(self.SampleIds):
            raise TableException, "Number of SampleIds differs from matrix size!"

        if n_samp != len(set(self.SampleIds)):
            raise TableException, "Duplicate SampleIds"

        if self.SampleMetadata is not None and \
           n_samp != len(self.SampleMetadata):
            raise TableException, "SampleMetadata not in a compatible shape \
                                   with data matrix!"

        if self.ObservationMetadata is not None and \
           n_obs != len(self.ObservationMetadata):
            raise TableException, "ObservationMetadata not in a compatible \
                                   shape with data matrix!"

    def _cast_metadata(self):
        """Casts all metadata to defaultdict to support default values"""
        default_samp_md = []
        default_obs_md = []
   
        if self.SampleMetadata is not None:
            for samp_md in self.SampleMetadata:
                d = defaultdict(lambda: None)
    
                if isinstance(samp_md, dict):
                    d.update(samp_md)
                elif samp_md is None:
                    pass
                else:
                    raise TableException, "Unable to cast metadata: %s" % \
                            repr(samp_md)

                default_samp_md.append(d)
            self.SampleMetadata = default_samp_md

        if self.ObservationMetadata is not None:
            for obs_md in self.ObservationMetadata:
                d = defaultdict(lambda: None)

                if isinstance(obs_md, dict):
                    d.update(obs_md)
                elif obs_md is None:
                    pass
                else:
                    raise TableException, "Unable to cast metadata: %s" % \
                            repr(obs_md)

                default_obs_md.append(d)
            self.ObservationMetadata = default_obs_md

    ### IS THIS WHAT WE WANT GETITEM/SETITEM AS?
    def __getitem__(self, args):
        """Passes through to internal matrix"""
        return self._data[args]

    def __setitem__(self, args, value):
        """Passes through to internal matrix"""
        self._data[args] = value
   
    ### DEFAULT OUTPUT SHOULD BE GTF
    def __str__(self):
        """Stringify self

        Default str output for a Table is just row/col ids and data values
        """
        return self.delimitedSelf()

    def delimitedSelf(self, delim='\t'):
        """Stringify self in a delimited form
        
        Default str output for the Table is just row/col ids and table data
        without any metadata
        """
        if self._biom_matrix_type is None:
            raise TableException, "Cannot delimit self if I don't have data..."

        samp_ids = delim.join(map(str, self.SampleIds))
        output = ['#RowIDs%s%s' % (delim, samp_ids)]
        
        for obs_id, obs_values in zip(self.ObservationIds, self._iter_obs()):
            str_obs_vals = delim.join(map(str, self._conv_to_np(obs_values)))
            output.append('%s%s%s' % (obs_id, delim, str_obs_vals))

        return '\n'.join(output)

    def __iter__(self):
        """Defined by subclass"""
        raise NotImplementedError

    def _iter_obs(self):
        """Defined by subclass"""
        raise NotImplementedError

    def _iter_samp(self):
        """Defined by subclass"""
        raise NotImplementedError

    def __eq__(self, other):
        """Equality is determined by the data matrix not metadata or IDs"""
        if self.ObservationIds != other.ObservationIds:
            return False
        if self.SampleIds != other.SampleIds:
            return False
        if self.ObservationMetadata != other.ObservationMetadata:
            return False
        if self.SampleMetadata != other.SampleMetadata:
            return False
        if not self._data_equality(other):
            return False

        return True

    def _data_equality(self,other):
        """Private method to determine equality of data

        The issue is that PySparse does not implement the __eq__ correctly
        """
        raise NotImplementedError

    def __ne__(self,other):
        return not (self == other)

    def _conv_to_np(self, v):
        """Convert values of v to numpy arrays"""
        raise NotImplementedError

    def copy(self):
        """Returns a copy of the Table"""
        #### NEEDS TO BE A DEEP COPY
        return self.__class__(self._data, self.SampleIds, self.ObservationIds,
                self.SampleMetadata, self.ObservationMetadata, self.TableId)

    def iterSamples(self, conv_to_np=True):
        """Yields (sample_values, sample_id, sample_metadata)

        NOTE: will return None in sample_metadata positions if 
        self.SampleMetadata is set to None
        """
        if self.SampleMetadata is None:
            samp_metadata = [None] * len(self.SampleIds)
        else:
            samp_metadata = self.SampleMetadata

        iterator = izip(self._iter_samp(), self.SampleIds, samp_metadata)
        for samp_v, samp_id, samp_md in iterator:
            if conv_to_np:
                yield (self._conv_to_np(samp_v), samp_id, samp_md)
            else:
                yield (samp_v, samp_id, samp_md)

    def iterObservations(self, conv_to_np=True):
        """Yields (observation_value, observation_id, observation_metadata)

        NOTE: will return None in observation_metadata positions if 
        self.ObservationMetadata is set to None
        """
        if self.ObservationMetadata is None:
            obs_metadata = [None] * len(self.ObservationIds)
        else:
            obs_metadata = self.ObservationMetadata
        
        iterator = izip(self._iter_obs(), self.ObservationIds, obs_metadata)
        for obs_v, obs_id, obs_md in iterator:
            if conv_to_np:
                yield (self._conv_to_np(obs_v), obs_id, obs_md)
            else:
                yield (obs_v, obs_id, obs_md)

    def filterSamples(self, f, invert=False):
        """Filter samples in self based on f
        
        f must accept three variables, the sample values, sample IDs and sample 
        metadata. The function must only return true or false.
        """
        samp_ids = []
        samp_vals = []
        samp_metadata = []

        # builtin filter puts all of this into memory and then return to the for
        # loop. This will impact memory substantially on large sparse matrices
        for s_val, s_id, s_md in self.iterSamples():
            if not xor(f(s_val, s_id, s_md), invert):
                continue

            # there is an implicit converstion to numpy types, want to make 
            # sure to convert back to underlying representation.
            samp_vals.append(self._conv_to_self_type(s_val))
            samp_metadata.append(s_md)
            samp_ids.append(s_id)
        
        # the additional call to _conv_to_self_type is to convert a list of 
        # vectors to a matrix
        # transpose is necessary as the underlying storage is sample == col
        return self.__class__(self._conv_to_self_type(samp_vals,transpose=True),
                samp_ids, self.ObservationIds, samp_metadata, 
                self.ObservationMetadata, self.TableId)

    def filterObservations(self, f, invert=False):
        """Filter observations in self based on f
        
        f must accept three variables, the observation values, observation ids
        and observation metadata. The function must only return true or false.
        """
        obs_ids = []
        obs_vals = []
        obs_metadata = []

        # builtin filter puts all of this into memory and then return to the for
        # loop. This will impact memory substantially on large sparse matrices
        for o_val, o_id, o_md in self.iterObservations():
            if not xor(f(o_val, o_id, o_md), invert):
                continue

            # there is an implicit converstion to numpy types, want to make 
            # sure to convert back to underlying representation.
            obs_vals.append(self._conv_to_self_type(o_val))
            obs_metadata.append(o_md)
            obs_ids.append(o_id)

        return self.__class__(self._conv_to_self_type(obs_vals),self.SampleIds,
                obs_ids, self.SampleMetadata, obs_metadata, self.TableId)

    def binSamplesByMetadata(self, f):
        """Yields tables by metadata
        
        f is given the sample metadata by row and must return what "bin" the
        sample is part of.
        """
        bins = {}
        # conversion of vector types is not necessary, vectors are not
        # being passed to an arbitrary function
        for samp_v, samp_id, samp_md in self.iterSamples(conv_to_np=False):
            bin = f(samp_md)

            # try to make it hashable...
            if not isinstance(bin, Hashable):
                bin = tuple(bin)

            if bin not in bins:
                bins[bin] = [[], [], []]

            bins[bin][0].append(samp_id)
            bins[bin][1].append(samp_v)
            bins[bin][2].append(samp_md)

        for bin, (samp_ids, samp_values, samp_md) in bins.iteritems():
            data = self._conv_to_self_type(samp_values, transpose=True)
            yield bin, self.__class__(data, samp_ids, self.ObservationIds, 
                            samp_md, self.ObservationMetadata, self.TableId)

    def binObservationsByMetadata(self, f):
        """Yields tables by metadata
        
        f is given the sample metadata by row and must return what "bin" the
        sample is part of.
        """
        bins = {}
        # conversion of vector types is not necessary, vectors are not
        # being passed to an arbitrary function
        for obs_v, obs_id, obs_md in self.iterObservations(conv_to_np=False):
            bin = f(obs_md)

            # try to make it hashable...
            if not isinstance(bin, Hashable):
                bin = tuple(bin)

            if bin not in bins:
                bins[bin] = [[], [], []]

            bins[bin][0].append(obs_id)
            bins[bin][1].append(obs_v)
            bins[bin][2].append(obs_md)

        for bin, (obs_ids, obs_values, obs_md) in bins.iteritems():
            yield bin, self.__class__(self._conv_to_self_type(obs_values), 
                            self.SampleIds, obs_ids, self.SampleMetadata,
                            obs_md, self.TableId)

    def transformSamples(self, f):
        """Apply a function to each sample
        
        f is passed a numpy vector and must return a vector
        """
        new_samp_v = []
        for samp_v, samp_id, samp_md in self.iterSamples():
            new_samp_v.append(self._conv_to_self_type(f(samp_v)))

        return self.__class__(self._conv_to_self_type(new_samp_v), 
                self.SampleIds, self.ObservationIds, self.SampleMetadata,
                self.ObservationMetadata, self.TableId)

    def transformObservations(self, f):
        """Apply a function to each observation

        f is passed a numpy vector and must return a vector
        """
        new_obs_v = []
        for obs_v, obs_id, obs_md in self.iterObservations():
            new_obs_v.append(self._conv_to_self_type(f(obs_v)))

        return self.__class__(self._conv_to_self_type(new_obs_v),
                self.SampleIds, self.ObservationIds, self.SampleMetadata,
                self.ObservationMetadata, self.TableId)

    #mergeTables, in place
        ### currently can only merge tables that do not have overlapping sample ids
        ### bail if overlapping

    def getBiomFormatObject(self):
        """Returns a dictionary representing the table in Biom format.

        This dictionary can then be easily converted into a JSON string for
        serialization.

        TODO: This method may be very inefficient in terms of memory usage, so
        it needs to be tested with several large tables to determine if
        optimizations are necessary or not (i.e. subclassing JSONEncoder, using
        generators, etc...).
        """
        if self._biom_type is None:
            raise TableException, "Unknown biom type"

        # Fill in top-level metadata.
        biom_format_obj = {}
        biom_format_obj["id"] = self.TableId
        biom_format_obj["format"] = get_biom_format_version_string()
        biom_format_obj["format_url"] =\
                get_biom_format_url_string()
        biom_format_obj["generated_by"] = "QIIME %s" %\
                get_qiime_library_version()
        biom_format_obj["date"] = "%s" % datetime.now().isoformat()

        # Determine what type of table we are (i.e. OTU table) and what type of
        # matrix we are storing (i.e. sparse or dense).
        try:
            num_rows, num_cols = self._data.shape
        except:
            num_rows = num_cols = 0
        hasData = True if num_rows > 0 and num_cols > 0 else False
        
        # Default the matrix element type to test to be an integer in case we
        # don't have any data in the matrix to test.
        test_element = 0
        if hasData:
            test_element = self[0,0]

        # Determine the type of elements the matrix is storing.
        if isinstance(test_element, int):
            dtype, matrix_element_type = int, "int"
        elif isinstance(test_element, float):
            dtype, matrix_element_type = float, "float"
        elif isinstance(test_element, str):
            dtype, matrix_element_type = str, "str"
        else:
            raise TableException("Unsupported matrix data type.")

        # Fill in details about the matrix.
        biom_format_obj["type"] = self._biom_type
        biom_format_obj["matrix_type"] = self._biom_matrix_type
        biom_format_obj["matrix_element_type"] = "%s" % matrix_element_type
        biom_format_obj["shape"] = [num_rows, num_cols]

        # Fill in details about the rows in the table and fill in the matrix's
        # data.
        biom_format_obj["rows"] = []
        biom_format_obj["data"] = []
        for obs_index, obs in enumerate(self.iterObservations()):
            biom_format_obj["rows"].append(
                    {"id" : "%s" % obs[1], "metadata" : obs[2]})
            # If the matrix is dense, simply convert the numpy array to a list
            # of data values. If the matrix is sparse, we need to store the
            # data in sparse format, as it is given to us in a numpy array in
            # dense format (i.e. includes zeroes) by iterObservations().
            if self._biom_matrix_type == "dense":
                # convert to python types, JSON doesn't like numpy types
                biom_format_obj["data"].append(map(dtype,obs[0]))
            elif self._biom_matrix_type == "sparse":
                dense_values = list(obs[0])
                sparse_values = []
                for col_index, val in enumerate(dense_values):
                    if int(val) != 0:
                        sparse_values.append([obs_index, col_index, val])
                biom_format_obj["data"].extend(sparse_values)

        # Fill in details about the columns in the table.
        biom_format_obj["columns"] = []
        for samp_index, samp in enumerate(self.iterSamples()):
            biom_format_obj["columns"].append(
                    {"id" : "%s" % samp[1], "metadata" : samp[2]})
        return biom_format_obj

    def getBiomFormatJsonString(self):
        """Returns a JSON string representing the table in Biom format."""
        return dumps(self.getBiomFormatObject())

class SparseTable(Table):
    _biom_matrix_type = "sparse"
    def __init__(self, *args, **kwargs):
        super(SparseTable, self).__init__(*args, **kwargs)
   
    def _data_equality(self, other):
        """Two pysparse matrices are equal if the items are equal"""
        return self._data.items() == other._data.items()

    def _conv_to_np(self, v):
        """Converts a vector to a numpy array

        Always returns a row vector for consistancy with numpy iteration over
        arrays
        """
        vals = v.items()

        num_rows, num_cols = v.shape
        if num_rows > num_cols:
            # i think pysparse matrices are always float...
            new_v = zeros(num_rows, dtype=float)
            for (row,col),val in vals:
                new_v[row] = val
        else:
            # i think pysparse matrices are always float...
            new_v = zeros(num_cols, dtype=float)
            for (row,col),val in vals:
                new_v[col] = val
        return new_v

    def _conv_to_self_type(self, vals, transpose=False):
        """For converting vectors to a compatible self type"""
        return to_ll_mat(vals, transpose)

    def __iter__(self):
        """Defined by subclass"""
        return self.iterSamples()

    def __eq__(self, other):
        """Equality is based on table data"""
        if self._data.items() == other._data.items():
            return True
        else:
            return False

    ### this method is type conversion heavy... but only fix if a burden when
    ### in use
    def _iter_samp(self):
        """Return sample vectors of data matrix vectors"""  
        rows, cols = self._data.shape
        for c in range(cols):
            # this pulls out col vectors but need to convert to the expected row
            # vector
            colvec = self._data[:,c]
            rowvec = ll_mat(1,rows)
            for (v_r, v_c), val in colvec.items():
                rowvec[0,v_r] = val
            yield rowvec

    def _iter_obs(self):
        """Return observation vectors of data matrix"""
        for r in range(self._data.shape[0]):
            yield self._data[r,:]

class DenseTable(Table):
    _biom_matrix_type = "dense"
    def __init__(self, *args, **kwargs):
        super(DenseTable, self).__init__(*args, **kwargs)

    def _data_equality(self, other):
        """Checks if the data matrices are equal"""
        return (self._data == other._data).all()

    def _conv_to_np(self, v):
        """Converts a vector to a numpy array"""
        return asarray(v)

    def _conv_to_self_type(self, vals, transpose=False):
        """For converting vectors to a compatible self type"""
        # expects row vector here...
        if transpose:
            return asarray(vals).T
        else:
            return asarray(vals)

    def __iter__(self):
        """Defined by subclass"""
        return self.iterSamples()

    def _iter_obs(self):
        """Return observations of data matrix"""
        for r in self._data:
            yield r

    def _iter_samp(self):
        """Return samples of data matrix in row vectors"""  
        for c in self._data.T:
            yield c

class OTUTable(object):
    _biom_type = "OTU table"
    pass

class DenseOTUTable(OTUTable, DenseTable):
    pass

class SparseOTUTable(OTUTable, SparseTable):
    pass

