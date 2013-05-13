.. _working_with_biom_tables:

=================================
Working with BIOM tables in QIIME
=================================

The BIOM (canonically pronounced *biome*) table is the core data type for downstream analyses in QIIME. It is a matrix of counts of observations on a per-sample basis. Most commonly, the observations are OTUs or taxa, and the samples are the unit of sampling in a study. These tables are often referred to as OTU tables in QIIME. To learn more about the BIOM format, you should review the `main documentation for the BIOM format <http://biom-format.org/>`_. 

This document covers the tools commonly used in QIIME for working with BIOM tables. Some of this code is in QIIME, and some is in the BIOM project, but since the BIOM project is a core dependency of QIIME, users won't necessarily know which is which. This document is intended to point users toward information on a number of scripts in QIIME. See the script help pages or call a script with the ``-h`` parameter to get information about how to work with it.

Filtering
=========

even sampling (rarefaction)
---------------------------

The `single_rarefaction.py <../scripts/single_rarefaction.html>`_ script allows the user to 

filtering observations (otus)
-----------------------------

filtering samples
-----------------

filtering taxa
--------------

Splitting and merging
=====================

splitting by sample metadata
----------------------------

splitting by taxonomy
---------------------

merging
-------

Miscellaneous 
=============

summarizing
-----------

adding metadata
---------------

sorting samples
---------------


Metadata description language
=============================

