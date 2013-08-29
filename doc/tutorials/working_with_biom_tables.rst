.. _working_with_biom_tables:

=================================
Working with BIOM tables in QIIME
=================================

The Biological Observation Matrix (or BIOM, canonically pronounced *biome*) table is the core data type for downstream analyses in QIIME. It is a matrix of counts of observations on a per-sample basis. Most commonly, the observations are OTUs or taxa, and the samples are the unit of sampling in a study (e.g., a microbiome sample from the skin of one individual at one time point). These tables are often referred to as OTU tables in QIIME (but really an OTU table is one type of a BIOM table). To learn more about the BIOM format, you should review the `documentation for the BIOM format <http://biom-format.org/>`_. 

This document covers some of the tools commonly used in QIIME for manipulating or summarizing BIOM tables. Some of this code is in QIIME, and some is in the BIOM project, but since the BIOM project is a core dependency of QIIME users won't necessarily know which is which. This document is intended to point users toward information on a number of scripts in QIIME. See the ``script help pages <../scripts/index.html>`_ or call a script with the ``-h`` parameter to get information about how to work with it.

Filtering
=========

Even sampling (rarefaction)
---------------------------

The `single_rarefaction.py <../scripts/single_rarefaction.html>`_ script allows the user to randomly subsample each of their samples to an equal count. This is essential for many diversity analyses, which assume that the count of observations for each sample is equal. `multiple_rarefactions.py <../scripts/multiple_rarefactions.html>`_ is a variant on this script which allows the user to generate rarified OTU tables at a range of different sampling depths (it is the equivalent of running `single_rarefaction.py <../scripts/single_rarefaction.html>`_ multiple times at different depths). `multiple_rarefactions_even_depth.py <../scripts/multiple_rarefactions_even_depth.html>`_ is another variant which allows the user to create many replicate subsampled OTU tables at a single depth.

Filtering observations/OTUs
-----------------------------

To filter observations (usually OTUs in QIIME) from a BIOM table based on their abundance, the number of samples they appear in, or by providing a list of OTUs that you want to remove (e.g., chimeric OTUs) or retain, you can use `filter_otus_from_otu_table.py <../scripts/filter_otus_from_otu_table.html>`_. You can find additional discussion of filtering observations in :ref:`filtering_contamination_otus`.

Filtering samples
-----------------

To filter samples from a BIOM table based on their total observation count, metadata associated with the samples, or by providing a list of samples you want to retain you can use `filter_samples_from_otu_table.py <../scripts/filter_samples_from_otu_table.html>`_. See :ref:`metadata_description` for information on how to specify which samples to retain or remove based on metadata.

Filtering taxa
--------------

To filter observations (usually OTUs in QIIME) from a BIOM table based on their taxonomic assignments, you can use `filter_taxa_from_otu_table.py <../scripts/filter_taxa_from_otu_table.html>`_. 

Splitting and merging
=====================

Splitting by sample metadata
----------------------------

To split a single BIOM table into multiple BIOM tables, each containing data on a subset of the samples based on their value for some sample metadata field (i.e., in the mapping file), you can use `split_otu_table.py <../scripts/split_otu_table.html>`_. For example, if you have a field in your mapping file called ``BodySite`` with values ``Gut``, ``Skin``, and ``Tongue``, calling ``split_otu_table.py -f BodySite`` will result in three BIOM tables, one each containing the samples at each body site. 

Splitting by taxonomy
---------------------

To split a single BIOM table into multiple BIOM tables, each containing data on a subset of the observations/OTUs based on their taxonomy, you can use `split_otu_table_by_taxonomy.py <../scripts/split_otu_table_by_taxonomy.html>`_. Here you will specify some taxonomic level, and an OTU table will be generated for each different taxonomic group at that level. This is useful, for example, to create per-domain or per-phylum OTU tables.

Merging
-------

To combine multiple BIOM tables into a single BIOM table, you can use `merge_otu_tables.py <../scripts/merge_otu_tables.html>`_. The main thing that you need to watch out for here is that the OTU ids and sample ids are compatible in each of the tables. If they are overlapping (e.g., you have ``OTU1`` in more than one table), their counts will be summed.

Miscellaneous 
=============

Summarizing
-----------

To summarize the information in a BIOM file, you can call `biom summarize-table <http://biom-format.org/documentation/summarizing_biom_tables.html>`_. This will provide you with information including the number of samples in the BIOM table, the number of observations/OTUs in the BIOM table, and a summary of the total counts observed in each sample.

Adding metadata
---------------

To add either sample or observation metadata to a BIOM file, you can call `add_metadata.py <http://biom-format.org/documentation/adding_metadata.html>`_. Sample metadata would be the information that is stored in your mapping file. It's useful to write this information to the BIOM file so all data and metadata can be stored in the same location for archiving or publication. In QIIME, observation metadata will often be the taxonomy associated with each OTU. 

Sorting samples
---------------

To sort the samples in the BIOM table based on some metadata about the samples (e.g, pH or time) or based on their order in a file, you can call `sort_otu_table.py <../scripts/sort_otu_table.html>`_. This can be useful if you're using the BIOM table with some code downstream that will plot some information about the samples in the order that they appear in the BIOM table.

