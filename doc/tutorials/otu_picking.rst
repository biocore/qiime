.. _otu_picking:

===============================
OTU picking strategies in QIIME
===============================

QIIME provides three high-level protocols for OTU picking. These can be described as de novo, closed-reference, and open-reference OTU picking, and are accessible through `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_, `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_, and `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_. Each of these protocols are described in this document, and commands are provided which illustrate how to run each of these with uclust and usearch 6.1 (i.e, usearch61).

Description of OTU picking processes
====================================

TODO: Discussion of protocols and pros and cons of each

Running the OTU picking workflows
=================================

The same workflow commands are used for running OTU picking with usearch61 and uclust. To run the methods with usearch, you will need to either pass in a parameters file or specify -m usearch61 on the command line, depending on what workflow you are using. See :ref:`qiime_parameter_files` for information on parameter files.

To run the methods with usearch, you will need to either pass in a parameters file or specify -m usearch61 on the command line, depending on what workflow you are using.

Conventions used in these examples
----------------------------------

It's a good idea, particularly for when running these workflows in parallel, to specify absolute paths for your input and output files. That is indicated here with ``$PWD``, but in practice it will often looks something like ``/home/ubuntu/my-analysis/seqs.fna``.

The reference-based OTU picking workflows require that the user provide reference files. Here we define some environment variables to point to those locations. These paths will likely be different on your system. You can download QIIME-compatible reference files from the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_. In this example we're working with the Greengenes 12_10 reference OTU collection. You can set environment variables to point to these as follows::

	export reference_seqs=/home/ubuntu/qiime_software/gg_otus-12_10-release/rep_set/97_otus.fasta
	export reference_tree=/home/ubuntu/qiime_software/gg_otus-12_10-release/trees/97_otus.tree
	export reference_tax=/home/ubuntu/qiime_software/gg_otus-12_10-release/taxonomy/97_otu_taxonomy.txt

De novo OTU picking
-------------------

With uclust::

	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/dn_uc/

With usearch61::
	
	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/dn_us/ -p $PWD/usearch_params.txt

where the following information is in ``usearch_params.txt``::
	
	pick_otus:otu_picking_method usearch61

The key output files are ``otu_table.biom``, the OTU table, and ``rep_set.tre``, the phylogenetic tree relating the OTUs in the OTU table.

You can find an additional example using de novo OTU picking in :ref:`tutorial`.

Closed-reference OTU picking
----------------------------

With uclust::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/cr_uc/ -r $reference_seqs -t $reference_tax

With usearch61::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/cr_us/ -r $reference_seqs -t $reference_tax -p $PWD/usearch_ref_params.txt

where the following information is in ``usearch_ref_params.txt``::
	
	pick_otus:otu_picking_method usearch61_ref

The key output file is ``otu_table.biom``, the OTU table. Note that there is no phylogenetic tree generated in this protocol - as all OTUs are defined by reference sequences, it is assumed that a tree already exists (which would likely be better than the one generated here).

You can find an additional example using closed-reference OTU picking in :ref:`illumina_overview_tutorial`.

Open-reference OTU picking
--------------------------

With uclust::

	pick_open_reference_otus.py -i seqs.fna -o or_uc/ -r $reference_seqs

With usearch61::

	pick_open_reference_otus.py -i seqs.fna -o or_us/ -r $reference_seqs -m usearch61

The key output files are ``otu_table.biom``, the OTU table, and ``rep_set.tre``, the phylogenetic tree relating the OTUs in the OTU table.

You can find an additional example using open-reference OTU picking in :ref:`open_reference_illumina`.

Alternative processing parameters
=================================

De-replication of sequences
--------------------------

If you're interested only in dereplicated sequences as your OTU picking process, that is a special case of de novo clustering where the similarity threshold is 100%. To achieve that you can do the following.

With uclust::
	
	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/derep_uc/ -p $PWD/uclust_dereplication_params.txt

where the following is in $PWD/uclust_dereplication_params.txt::
	
	pick_otus:similarity 1.0

With usearch61::
	
	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/derep_us/ -p $PWD/usearch_dereplication_params.txt

where the following information is in ``usearch_dereplication_params.txt``::
	
	pick_otus:otu_picking_method usearch61
	pick_otus:similarity 1.0

Running usearch in size-order mode
----------------------------------

If you're interested in running the usearch OTU pickers in size-order mode (meaning that accepts are prioritized by the size of the cluster rather than the percent identity), add the following lines to a parameters file::

	pick_otus:sizeorder True 
	pick_otus:maxaccepts 16
	pick_otus:maxrejects 64

For example, in de novo mode::

	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/dn_us_sizeorder/ -p $PWD/dn_sizeorder_params.txt

where the following information is in ``dn_sizeorder_params.txt``::
	
	pick_otus:otu_picking_method usearch61
	pick_otus:sizeorder True 
	pick_otus:max_accepts 16
	pick_otus:max_rejects 64

In closed-reference mode::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -o $PWD/cr_us_sizeorder/ -r $reference_seqs -t $reference_tax -p $PWD/cr_sizeorder_params.txt

where the following information is in ``cr_sizeorder_params.txt``::
	
	pick_otus:otu_picking_method usearch61_ref
	pick_otus:sizeorder True 
	pick_otus:max_accepts 16
	pick_otus:max_rejects 64

In open-reference mode::

	pick_open_reference_otus.py -i seqs.fna -o or_us_sizeorder/ -r $reference_seqs -m usearch61 -p $PWD/or_sizeorder_params.txt

where the following information is in ``or_sizeorder_params.txt``::
	
	pick_otus:sizeorder True 
	pick_otus:max_accepts 16
	pick_otus:max_rejects 64



Citing these tools
==================

If using these tools you should cite both QIIME and usearch or uclust. 