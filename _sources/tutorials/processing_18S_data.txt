.. _processing_18S_data:

====================
Analysis of 18S data
====================

Most of the steps for analysis of 18S, or mixed 16S/18S, are identical to the standard 16S pipeline described in the `QIIME Illumina Overview Tutorial <./illumina_overview_tutorial.html>`_ or the `QIIME 454 Overview Tutorial <./tutorial.html>`_, with the main difference being the use of a non-default reference database. Alternative QIIME-compatible reference databases can be found on the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_. Please refer to the `Fungal ITS Analysis Tutorial <./fungal_its_analysis.html>`_ for an example of how to perform a QIIME analysis with a non-default reference database.

Separating OTU Tables According to Domain
=========================================

It may be desirable to split the OTU table according to domain for mixed 16S/18S datasets.  To do this, you can use ``split_otu_table_by_taxonomy.py`` to split at the domain level (2)::

	split_otu_table_by_taxonomy.py -i otu_table.biom -L 2 -o split_otu_tables/

The output directory, ``split_otu_tables``, will contain an OTU table for archaea, bacteria, and eukaryotes, which can be utilized in downstream QIIME analyses.
