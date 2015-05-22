.. _otu_picking:

===============================
OTU picking strategies in QIIME
===============================

QIIME provides three high-level protocols for OTU picking. These can be described as de novo, closed-reference, and open-reference OTU picking, and are accessible through `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_, `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_, and `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_. Each of these protocols are briefly described in this document; for a more detailed discussion of these OTU picking protocols, please see `Rideout et al. (2014) <https://peerj.com/articles/545/>`_.

**Open-reference OTU picking with** ``pick_open_reference_otus.py`` **is the preferred strategy for OTU picking among the QIIME developers.**

.. note:: QIIME does not actually implement OTU picking algorithms, but rather wraps external OTU clustering tools. For this reason, it is important to cite the OTU clustering tools that you used directly, in addition to citing QIIME. There are a number of OTU clustering tools available through QIIME's workflows, including open source (e.g., `SortMeRNA <http://www.ncbi.nlm.nih.gov/pubmed/23071270>`_, `SUMACLUST <http://metabarcoding.org/>`_, and `swarm <https://peerj.com/articles/593/>`_) and closed source tools (e.g., `uclust and usearch <http://www.ncbi.nlm.nih.gov/pubmed/20709691>`_). ``uclust`` is the default OTU clustering tool used in QIIME's workflows. We are currently evaluating changing the default OTU clustering tool to one of the open source alternatives for future versions of QIIME.

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

You **cannot** use closed-reference OTU picking if:

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

You **cannot** use open-reference OTU picking if:

*  You are comparing non-overlapping amplicons, such as the V2 and the V4 regions of the 16S rRNA.
*  You do not have a reference sequence collection to cluster against, for example because you're working with an infrequently used marker gene.

Pros:

*  All reads are clustered.
*  Speed. Open-reference OTU picking is partially run in parallel. In particular, the *subsampled open reference OTU picking* process implemented in ``pick_open_reference_otus.py`` is much faster than ``pick_de_novo_otus.py`` as some strategies are applied to run several pieces of the workflow in parallel.

Cons:

*  Speed. Some steps of this workflow do still run serially. For data sets with a lot of novel diversity with respect to the reference sequence collection, this can still take days to run.

Running the OTU picking workflows
=================================

Please refer to the script usage examples in `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_, `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_, and `pick_open_reference_otus.py <../scripts/pick_open_reference_otus.html>`_, and the `QIIME Illumina Overview Tutorial <./illumina_overview_tutorial.html>`_ and the `QIIME 454 Overview Tutorial <./tutorial.html>`_ for examples of how to use QIIME's OTU picking workflows.

Alternative processing parameters
=================================

Dereplication of sequences
---------------------------

If you're interested only in dereplicating sequences as your OTU picking process, that is a special case of de novo clustering where the similarity threshold is 100%. To achieve that you can do the following::

	pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/derep_uc/ -p $PWD/dereplication_params.txt

where the following is in $PWD/dereplication_params.txt::

	pick_otus:similarity 1.0

Running usearch in size-order mode
----------------------------------

If you're interested in running the usearch OTU pickers in size-order mode (meaning that accepts are prioritized by the size of the cluster rather than the percent identity), add the following lines to a parameters file::

	pick_otus:otu_picking_method usearch61
	pick_otus:sizeorder True
	pick_otus:maxaccepts 16
	pick_otus:maxrejects 64

Pass this parameters file via ``-p`` to any of the three OTU picking workflows in QIIME.
