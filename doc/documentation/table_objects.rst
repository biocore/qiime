.. _table_objects:

===========================================
QIIME Table Objects
===========================================

As of 1.4.0-dev, QIIME began providing rich ``Table`` objects to facilitate analyses. The objects encapsulate matrix data (such as OTU counts) and abstract the interaction away from the programmer. This provides the immediate benefit of the programmer not having to worry about what the underlying data object is, and in turn allows for different data representations to be supported. Currently, QIIME supports a ``dense`` object built off of ``numpy.array`` (`NumPy <http://http://numpy.scipy.org/>`_) and a ``sparse`` object built off of ``pysparse.spmatrix.ll_mat`` (`PySparse <http://pysparse.sourceforge.net/>`_). 

Motivation for the objects
==========================

Prior to QIIME 1.4.0-dev, all OTU tables used a dense representation in memory. OTU tables tend to be highly "sparse", or contain a significant number of entries with zero counts. These zero counts consume a dramatic amount of memory as the tables grow in size. Because of this, the QIIME team began investigating sparse matrix representations that only store the nonzero values. 

We decided early on to couple the support of the sparse data representation with rich objects that express the complex OTU table datatype. The refactoring comes with two specific benefits: a common API (application programming interface) that allows the developers easier access to common OTU table functionality and an abstraction from the underlying matrix data.

Description of available ``Table`` objects
==========================================

There are multiple objects available but some of them are unofficial abstract base classes (does not use the ``abc`` module for historical reasons). In practice, the objects used should be the derived Tables such as ``SparseOTUTable`` or ``DenseAbundanceTable``. 

Abstract base classes
---------------------

Abstract base classes establish standard interfaces for subclassed types and provide common functionality for derived types. 

``Table``
^^^^^^^^^

``` is a container object and an abstract base class that provides a common and required API for subclassed objects. Through the use of private interfaces, it is possible to create public methods that operate on the underlying datatype without having to implement each method in each subclass. For instance, ``Table.iterSamplesData`` will return a generator that always yields ``numpy.array`` vectors for each sample regardless of how the table data is actually stored. This functionality results from derived classes implementing private interfaces, such as ``Table._conv_to_np``.

``OTUTable``
^^^^^^^^^^^^

The ``OTUTable`` base class provides functionality specific for OTU tables. Currently, it only provides a static private member variable that describes its ``BIOM`` type. This object was stubbed out incase future methods are developed that do not make sense with the context of, say, an MG-RAST metagenomic abundance table. It is advised to always use an object that subclasses ``OTUTable`` if the analysis is on OTU data.

``AbundanceTable``
^^^^^^^^^^^^^^^^^^

The ``AbundanceTable`` base class provides functionality for abundance tables such as those returned by MG-RAST. Currently, the object only defines a single static private member variable that describes its ``BIOM`` type. The object was stubbed out in case abundance table specific functionality is developed in the future. It is advised to always use an object that subclasses ``AbundanceTable`` if the analysis is on abundance data.

Container classes
-----------------
The container classes implement required private member variable interfaces as defined by the ``Table`` abstract base class. Specifically, these objects define the ways in which data is moved into and out of the contained data object. These are fully functional and usable objects, however they do not implement table type specifc functionality.

``SparseTable``
^^^^^^^^^^^^^^^

The subclass ``SparseTable`` can be derived for use with table data. This object implemented all of the required private interfaces specified by the ``Table`` base class. The object contains a ``_data`` private member variable that is an instance of ``pysparse.spmatrix.ll_mat``. The ``ll_mat`` object is a linked-list structure that represensts only non-zero values. The object maintains the locations of non-zero values for quick lookups and returns zero when a element zero'd element is requested. It is advised to used derived objects of SparseTable if the data being operated on is sparse.

``DenseTable``
^^^^^^^^^^^^^^

The ``DenseTable`` object fulfills all private member methods stubbed out by the ``Table`` base class. The dense table contains a private member variable that is an instance of ``numpy.array``. The ``array`` object is a matrix that contains all values including zeros. It is advised to use this table only if the number of samples and observations is reasonble. Unfortunately, it isn't reasonable to define reasonable in this context. However, if either the number of observations or the number of samples is > 1000, it would probably be a good idea to rely on a ``SparseTable``.

Table type objects
------------------

The table type objects define variables and methods specific to a table type. Under the majority of situations, these are the objects that should be instantiated.

``DenseOTUTable``
^^^^^^^^^^^^^^^^^

A dense representation of an OTU table.

``SparseOTUTable``
^^^^^^^^^^^^^^^^^^

A sparse representation of an OTU table.

``DenseAbundanceTable``
^^^^^^^^^^^^^^^^^^^^^^^

A dense representation of an abundance table.

``SparseAbundanceTable``
^^^^^^^^^^^^^^^^^^^^^^^^

A sparse representation of an abundance table.

Description of ``Table`` API
============================

Public Variables
----------------

``ObservationIds``
^^^^^^^^^^^^^^^^^^^^^^^^

``ObservationIds`` is a list of strings containing observation identifiers. The order of this list matters: the index position of an id corresponds to its position within the underlying data matrix and its metadata. Observation ids are required and there must be one id for every observation within the data matrix. 

``SampleIds``
^^^^^^^^^^^^^^^^^^^

``SampleIds`` is a list of strings containing sample identifiers. The order of this list matters: the index position of an id corresponds to its position within the underlying data matrix and its metadata. Sample ids are required and there must be one id for every sample within the data matrix.

``ObservationMetadata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ObservationMetadata`` is a list of dicts. Each dict contains metadata associated to the observation id at the same index position. Observation metadata is optional. If metadata is not present, this variable will be set to None.

``SampleMetadata``
^^^^^^^^^^^^^^^^^^^^^^^^

``SampleMetadata`` is a list of dicts. Each dict contains metadata associated to the sample id at the same postion. Sample metadata is optional. If metadata is not present, this variable will be set to None.

``TableId``
^^^^^^^^^^^^^^^^^

``TableId`` is an optional variable that can be used to indentify a table. The value stored in this variable is written to the ``id`` field within a ``BIOM`` file. 

Private variables
-----------------

``_biom_matrix_type``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The BIOM matrix type can either be 'sparse' or 'dense'. This information is used by ``Table.getBiomFormatObject`` in order to figure out what ``matrix_type`` to put to file.

``_biom_type``
^^^^^^^^^^^^^^^^^^^^

The BIOM type defines the table "type" within a BIOM file.

``_data``
^^^^^^^^^^^^^^^

The underlying data object.

``_obs_index``
^^^^^^^^^^^^^^^^^^^^

A helper lookup dictionary that is {'observation_id': index_in_data}.

``_sample_index``
^^^^^^^^^^^^^^^^^^^^^^^

A helper lookup dictionary that is {'sample_id': index_in_data}.

Public Methods
--------------

``binObservationsByMetadata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yields tables by metadata. A user supplied function ``f`` is given the observation  metadata must return what "bin" the observation is part of.

``binSamplesByMetadata``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yields tables by metadata. A user supplied function ``f`` is given the sample metadata by row and must return what "bin" the sample is part of.

``copy``
^^^^^^^^^^^^^^

Returns a shallow copy of the Table

``delimitedSelf``
^^^^^^^^^^^^^^^^^^^^^^^

Stringify self in a delimited form. Default str output for the ``Table`` is just row/col ids and table data without any metadata. This is the "classic" table. If ``header_key`` is not None, it tries to pull out that key from observation metadata. If ``header_value`` is not None, use the ``header_value`` in the output as a column id. ``metadata_formatter`` is a function which takes a metadata entry and returns a formatted version that should be written to file.

``filterObservations``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Filter observations in self based on a user supplied function. ``f`` must accept three variables, the observation values, observation ids and observation metadata. The function must only return True or False.

``filterSamples``
^^^^^^^^^^^^^^^^^^^^^^^

Filter samples in self based on a user supplied function. ``f`` must accept three variables, the sample values, sample IDs and sample metadata. The function must only return true or false.
        
``getBiomFormatJsonString``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns a JSON string representing the table in BIOMformat.

``getBiomFormatObject``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns a dict representing the table in BIOM format. This dictionary can then be easily converted into a JSON string for serialization.

``getBiomFormatPrettyPrint``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns a 'pretty print' format of a BIOM file. WARNING: This method displays data values in a columnar format and can be misleading.

``getObservationIndex``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns the observation index position in _data.

``getSampleIndex``
^^^^^^^^^^^^^^^^^^^^^^^^

Returns the sample index position in _data.

``getValueByIds``
^^^^^^^^^^^^^^^^^^^^^^^

Return the value in the matrix corresponding to (obs_id, samp_id)

``isEmpty``
^^^^^^^^^^^^^^^^^

Returns True if ``Table._data`` is empty. It should not be possible for ``Table._data`` to be empty and have ``Table.SampleIds`` or ``Table.ObservationIds`` populated.

``iterObservationData``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yields vectors of observation data where the values at each index correspond to the ``Table.SampleId`` and ``Table.SampleMetadata`` at the same index. 

``iterObservations``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Yields (observation_value, observation_id, observation_metadata) NOTE: will return None in observation_metadata positions if ``Table.ObservationMetadata`` is set to None.
        
``iterSampleData``
^^^^^^^^^^^^^^^^^^^^^^^^

Yields vectors of sample sample data where the values at each index correspond to the ``Table.ObservationId`` and ``Table.ObservationMetadata`` at the same index.

``iterSamples``
^^^^^^^^^^^^^^^^^^^^^

Yields (sample_values, sample_id, sample_metadata)/ NOTE: will return None in sample_metadata positions if ``Table.SampleMetadata`` is set to None

``merge``
^^^^^^^^^^^^^^^

Merge two tables together. The axes, samples and observations, can be controlled independently and can both work on either 'union' or 'intersection'.  ``merge_f`` is a function that takes two arguments and returns a value. The method is parameterized so that the programmer can define how values are handled where there is overlap in (sample_id, observation_id) values between tables. ``sample_metadata_f`` and ``observation_metadata_f`` define how to merge metadata between tables. The default is to prefer the metadata associated to self if self has metadata otherwise take metadata from other. These functions are given both metadata dicts and must return a single metadata dict. NOTE: There is an implicit type conversion to float. Tables using strings as the type are not supported but no active check is in place. NOTE: The return type is always that of self

``nonzero``
^^^^^^^^^^^^^^^^^

Returns nonzero locations within the data matrix. The values returned are (observation_id, sample_id).

``normObservationBySample``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Return new table with relative abundance in each sample.

``normSampleByObservation``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Return new table with relative abundance in each observation.

``observationData``
^^^^^^^^^^^^^^^^^^^^^^^^^

Return a numpy vector with samples values associated to an observation id.

``observationExists``
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Returns True if observation exists, False otherwise.

``reduce``
^^^^^^^^^^^^^^^^

Reduce over axis with ``f``. Axis can be either 'sample' or 'observation'

``sampleData``
^^^^^^^^^^^^^^^^^^^^

Return a numpy vector with observation values associated to a sample id.

``sampleExists``
^^^^^^^^^^^^^^^^^^^^^^

Returns True if sample exists, False otherwise.

``setValueByIds``
^^^^^^^^^^^^^^^^^^^^^^^

Set the value in the matrix corresponding to (observation_id, sample_id).

``sortByObservationId``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Return a table with the observation ids sorted by a user supplied function. The default is natural sort.

``sortBySampleId``
^^^^^^^^^^^^^^^^^^^^^^^^

Return a table with the sample ids sorted by a user supplied function/

``sortObservationOrder``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Return a new table with the observation ids in the order given.

``sortSampleOrder``
^^^^^^^^^^^^^^^^^^^^^^^^^

Return a new table with the sample ids in the order given.

``sum``
^^^^^^^^^^^^^

Returns the sum by axis. The axis can either be 'whole', 'sample' or 'observation'. For the 'sample' and 'observation' axes, a vector is returned with the sum of the orthoganl vector. For example, if ``sum`` is called with 'sample', a vector is returned with a sum of the observations within each sample. Index 0 in that resulting vector would correspond to ``Table.SampleIds[0]``
        
``transformObservations``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Apply a user defined function to each observation. ``f`` is passed a numpy vector and must return a numpy vector of the same shape and datatype.
        
``transformSamples``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Apply a user defined function to each sample. ``f`` is passed a numpy vector and must return a numpy vector of the same shape and datatype.

Private methods
---------------

``_cast_metadata``
^^^^^^^^^^^^^^^^^^^^^^^^

Casts all metadata to defaultdict to support default values

``_conv_to_np``
^^^^^^^^^^^^^^^^^^^^^

Converts a vector to a numpy array. Always returns a row vector for consistancy with numpy iteration over arrays
        

``_conv_to_self_type``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For converting vectors to a compatible self type.

``_data_equality``
^^^^^^^^^^^^^^^^^^^^^^^^

A private method that defines how to test equality between Table._data variables.

``_index_ids``
^^^^^^^^^^^^^^^^^^^^

Sets lookups {id:index in _data}

``_intersect_id_order``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Determines the merge order for id lists A and B

``_iter_obs``
^^^^^^^^^^^^^^^^^^^

Return observation vectors of data matrix

``_iter_samp``
^^^^^^^^^^^^^^^^^^^^

Return sample vectors of data matrix vectors

``_union_id_order``
^^^^^^^^^^^^^^^^^^^^^^^^^

Determines merge order for id lists A and B

``_verify_metadata``
^^^^^^^^^^^^^^^^^^^^^^^^^^

Obtain some notion of sanity on object construction with inputs

Mini How-To
===========


