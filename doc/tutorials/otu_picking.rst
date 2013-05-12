.. _otu_picking:

===============================
OTU picking strategies in QIIME
===============================

QIIME provides three high-level protocols for OTU picking. These can be described as de novo, closed-reference, and open-reference OTU picking, and are accessible through `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_, `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_, and `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_. Each of these protocols are described in this document, and commands are provided which illustrate how to run each of these with uclust and usearch 6.1 (i.e, usearch61).

Description of QIIME's OTU picking protocols
============================================

De novo OTU picking
-------------------

In a de novo OTU picking process, reads are clustered against one another without any external reference sequence collection. ``pick_de_novo_otus.py`` is the primary interface for de novo OTU picking in QIIME, and includes taxonomy assignment, sequence alignment, and tree-building steps. A benefit of de novo OTU picking is that all reads are clustered. A drawback is that there is no existing support for running this in parallel in QIIME, so it can be too slow to apply to large datasets (e.g., more than 10 million reads). 

You **must** use de novo OTU picking if:

*  You do not have a reference sequence collection to cluster against, for example because you're working with an infrequently used marker gene.

You **cannot** use de novo OTU picking if:

*  You are comparing non-overlapping amplicons, such as the V2 and the V4 regions of the 16S rRNA.
*  You working with very large data sets, like a full HiSeq 2000 run. (Technically, you can use de novo OTU picking here, but you literally might wait a month for ``pick_de_novo_otus.py`` to run.)

Pros:

*  All reads are clustered

Cons:

*  Speed. Does not run in parallel.

Closed-reference OTU picking
----------------------------

In a closed-reference OTU picking process, reads are clustered against a reference sequence collection and any reads which do not hit a sequence in the reference sequence collection are excluded from downstream analyses. ``pick_closed_reference_otus.py`` is the primary interface for closed-reference OTU picking in QIIME. If the user provides taxonomic assignments for sequences in the reference database, those are assigned to OTUs.

You **must** use closed-reference OTU picking if:

*  You are comparing non-overlapping amplicons, such as the V2 and the V4 regions of the 16S rRNA. Your reference sequences must span both of the regions being sequenced.

You **cannot** use de novo OTU picking if:

*  You do not have a reference sequence collection to cluster against, for example because you're working with an infrequently used marker gene.

Pros:

*  Speed. Closed-reference OTU picking is fully parallelizable, so is useful for extremely large data sets.
*  Better trees and taxonomy. Because all OTUs are already defined in your reference sequence collection, you may already have a tree and a taxonomy that you trust for those OTUs. You have the option of using those, or building a tree and taxonomy from your sequence data.

Cons:

*  Inability to detect novel diversity with respect to your reference sequence collection. Because reads that don't hit the reference sequence collection are discarded, your analyses only focus on the diversity that you "already know about". Also, depending on how well-characterized the environment that you're working in is, you may end up throwing away a small fraction of your reads (e.g., discarding 1-10% of the reads is common for 16S-based human microbiome studies, where databases like Greengenes cover most of the organisms that are typically present) or a large fraction of your reads (e.g, discarding 50-80% of the reads has been observed for "unusual" environments like the Guerrero Negro microbial mats). 

Open-reference OTU picking
--------------------------

In an open-reference OTU picking process, reads are clustered against a reference sequence collection and any reads which do not hit the reference sequence collection are subsequently clustered de novo. ``pick_open_reference_otus.py`` is the primary interface for open-reference OTU picking in QIIME, and includes taxonomy assignment, sequence alignment, and tree-building steps.

**Open-reference OTU picking with** ``pick_open_reference_otus.py`` **is the preferred strategy for OTU picking among the QIIME developers.**

You **cannot** use de novo OTU picking if:

*  You are comparing non-overlapping amplicons, such as the V2 and the V4 regions of the 16S rRNA.
*  You do not have a reference sequence collection to cluster against, for example because you're working with an infrequently used marker gene.

Pros:

*  All reads are clustered.
*  Speed. Open-reference OTU picking is partially run in parallel. In particular, the *subsampled open reference OTU picking* process implemented in ``pick_open_reference_otus.py`` is much faster than ``pick_de_novo_otus.py`` as some strategies are applied to run several pieces of the workflow in parallel.

Cons:

*  Speed. Some steps of this workflow do still run serially. For data sets with a lot of novel diversity with respect to the reference sequence collection, this can still take days to run.

Running the OTU picking workflows
=================================

The same workflow commands are used for running OTU picking with usearch61 and uclust. To run the methods with usearch, you will need to either pass in a parameters file or specify ``-m usearch61`` on the command line, depending on what workflow you are using. See :ref:`qiime_parameter_files` for information on parameter files.

To run the methods with usearch, you will need to either pass in a parameters file or specify -m usearch61 on the command line, depending on what workflow you are using.

Conventions used in these examples
----------------------------------

It's a good idea, particularly for when running these workflows in parallel, to specify absolute paths for your input and output files. That is indicated here with ``$PWD``, but in practice it will often looks something like ``$HOME/my-analysis/seqs.fna``.

The reference-based OTU picking workflows require that the user provide reference files (the reference sequence collection). Here we define some environment variables to point to those locations. These paths will likely be different on your system. You can download QIIME-compatible reference files from the `QIIME resources page <http://qiime.org/home_static/dataFiles.html>`_. In this example we're working with the Greengenes 12_10 reference OTU collection. You can set environment variables to point to these as follows::

	export QIIME_DIR=$HOME/qiime_software
	export reference_seqs=$QIIME_DIR/gg_otus-12_10-release/rep_set/97_otus.fasta
	export reference_tree=$QIIME_DIR/gg_otus-12_10-release/trees/97_otus.tree
	export reference_tax=$QIIME_DIR/gg_otus-12_10-release/taxonomy/97_otu_taxonomy.txt

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

Open-reference OTU picking
--------------------------

With uclust::

	pick_open_reference_otus.py -i seqs.fna -o or_uc/ -r $reference_seqs

With usearch61::

	pick_open_reference_otus.py -i seqs.fna -o or_us/ -r $reference_seqs -m usearch61

The key output files are ``otu_table.biom``, the OTU table, and ``rep_set.tre``, the phylogenetic tree relating the OTUs in the OTU table.

You can find an additional example using open-reference OTU picking in :ref:`illumina_overview_tutorial`.

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
