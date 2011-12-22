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

There are multiple objects available but some of them are unofficial abstract base classes (does not use the ``abc`` module for historical reasonas).

Table
-----

``Table`` is a container object and an abstract base class that provides a common and required API for subclassed objects. Through the use of private interfaces, it is possible to create public methods that operate on the underlying datatype without having to implement each method in each subclass. For instance, ``Table.iterSamplesData`` will return a generator that always yields ``numpy.array`` vectors for each sample regardless of how the table data is actually stored. This functionality results from derived classes implementing private interfaces, such as ``Table._conv_to_np``.

OTUTable
--------

The ``OTUTable`` base class provides functionality specific for OTU tables. Currently, it only provides a static private member variable that describes its ``BIOM`` type. This object was stubbed out incase future methods are developed that do not make sense with the context of, say, an MG-RAST metagenomic abundance table. It is advised to always an object that subclasses ``OTUTable`` if the analysis is on OTU data.

SparseTable
-----------

The subclass ``SparseTable`` can be derived for use with table data. This object implemented all of the required private interfaces specified by the ``Table`` base class. However, the ``SparseTable`` does not subclass
Description of the base class public API
======================================================

Description of OTUtable specific interfaces
===========================================

Mini How-To
===========


