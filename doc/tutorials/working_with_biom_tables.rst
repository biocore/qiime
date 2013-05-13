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

The `single_rarefaction.py <../scripts/single_rarefaction.html>`_ script allows the user to randomly subsample each of their samples to an equal count. This is essential for many diversity analyses, which assume that the count of observations for each sample is equal. `multiple_rarefactions.py <../scripts/multiple_rarefactions.html>`_ is a variant on this script which allows the user to generate rarified OTU tables at a range of different sampling depths (it is the equivalent of running `single_rarefaction.py <../scripts/single_rarefaction.html>`_ multiple times at different depths). `multiple_rarefactions_even_depth.py <../scripts/multiple_rarefactions_even_depth.html>`_ is another variant which allows the user to create many replicate subsampled OTU tables at a single depth.

filtering observations (otus)
-----------------------------

To filter observations (usually OTUs in QIIME) from a BIOM table based on their abundance, the number of samples they appear in, or by providing a list of OTUs you want to remove (e.g., chimeric OTUs) or retain, you can use `filter_otus_from_otu_table.py <../scripts/filter_otus_from_otu_table.html>`_.

filtering samples
-----------------

To filter samples from a BIOM table based on their total observation count, metadata associated with the samples, or by providing a list of samples you want to retain you can use `filter_samples_from_otu_table.py <../scripts/filter_samples_from_otu_table.html>`_. See the **Metadata description language** section for information on how to specify which samples to retain or remove based on metadata.

filtering taxa
--------------

To filter observations (usually OTUs in QIIME) from a BIOM table based on their taxonomic assignments, you can use `filter_taxa_from_otu_table.py <../scripts/filter_taxa_from_otu_table.html>`_. 

Splitting and merging
=====================

splitting by sample metadata
----------------------------

To split a single BIOM table into multiple BIOM tables based on their value for some sample metadata field (i.e., in the mapping file) you can use `split_otu_table.py <../scripts/split_otu_table.html>`_. For example, if you have a field in your mapping file called ``BodySite`` with values ``Gut``, ``Skin``, and ``Tongue``, calling ``split_otu_table.py -f BodySite`` will result in three BIOM tables, one each containing the samples at each body site. 

splitting by taxonomy
---------------------

To split a single BIOM table into multiple BIOM tables based on taxonomy of the OTUs you can use `split_otu_table_by_taxonomy.py <../scripts/split_otu_table_by_taxonomy.html>`_. Here you will specify some taxonomic level, and generate an OTU table for each different taxa at that level. This is useful, for example, to create per-phylum OTU tables.

merging
-------

To combine multiple BIOM tables into a single BIOM table, you can use `merge_otu_tables.py <../scripts/merge_otu_tables.html>`_. The main thing that you need to watch out for here is that the OTU ids and sample ids are compatible. If they are overlapping (e.g., you have ``OTU1`` in more than one table), their counts will be merged. 

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

