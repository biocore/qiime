.. _table_objects:

===========================================
biom-format table objects
===========================================

As of 1.4.0-dev, QIIME began using rich ``Table`` objects from the `biom-format project <http://www.biom-format.org>`_ to store OTU tables. The objects simplify many programming tasks associated with OTU tables. Full documentation of the biom-format ``Table`` objects can be found `here <http://www.biom-format.org/documentation/table_objects.html>`_.

Motivation for the objects
==========================

Prior to QIIME 1.4.0-dev, all OTU tables used a dense representation in memory. OTU tables tend to be highly "sparse", or contain a significant number of entries with zero counts. These zero counts consume a dramatic amount of memory as the tables grow in size. Because of this, the QIIME team began investigating sparse matrix representations that only store the nonzero values. 

We decided early on to couple the support of the sparse data representation with rich objects that express the complex OTU table datatype. The refactoring comes with two specific benefits: a common API (application programming interface) that allows the developers easier access to common OTU table functionality, and an abstraction from the underlying matrix data.

