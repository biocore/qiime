.. _chimera_checking:

=====================================
Chimera checking sequences with QIIME
=====================================

ChimeraSlayer 
-------------

This document describes how to chimera check sequences with QIIME via ChimeraSlayer, and subsequently exclude chimeric sequences when building your tree and OTU table so they will not be included in downstream analyses. (Note that because we recommend applying chimera checking after OTU and representative sequence picking, for the purposes of this document `chimeric sequence` should be considered synonymous to `chimeric OTU`.)

ChimeraSlayer requires aligned sequences, you should apply chimera checking after you've aligned your sequences with PyNAST, but before you apply the lanemask filtering to your alignment.

Applying ChimeraSlayer is done as follows::

	identify_chimeric_seqs.py -m ChimeraSlayer -i rep_set_aligned.fasta -a reference_set_aligned.fasta -o chimeric_seqs.txt

Before building your phylogenetic tree, you should now remove chimeric sequences from your alignment using your chimeric sequence list (``chimeric_seqs.txt`` in this example) and the alignment to as follows::

	filter_fasta.py -f rep_set_aligned.fasta -o non_chimeric_rep_set_aligned.fasta -s chimeric_seqs.txt -n
	
Don't forget to pass the ``-n`` parameter to ``filter_fasta.py`` -- this tells ``filter_fasta.py`` to discard the sequences you've passed via ``-s``, rather than to keep only those sequences. Then pass the resulting subalignment (``non_chimeric_rep_set_aligned.fasta``) in the downstream steps ``filter_alignment.py`` prior to tree-building.

You'll also need to exclude the chimeric sequences at the ``make_otu_table.py`` step. To do that, you pass the chimeric sequence list to ``make_otu_table.py`` via the ``-e`` parameter::

	make_otu_table.py -i otu_map.txt -o otu_table.biom -e chimeric_seqs.txt -t taxonomy.txt

This will cause the chimeric sequences to be excluded from all downstream analysis that you'll be performing, as the OTU table is always the input provided (for example to beta diversity, alpha rarefaction, etc.).

USEARCH 6.1
-----------

USEARCH 6.1 performs reference based chimera detection, like ChimeraSlayer, but also can perform de novo chimera detection based upon abundances of input sequences.

Important notes about USEARCH 6.1 usage:

.. note::

   * 1.  Input sequences should be demultiplexed sequences (e.g. output of ``split_libraries.py``) that are not already clustered.
   * 2.  The reference database should not be aligned.
   * 3.  The reference sequences need to be in the same orientation as the query sequences.  Use ``adjust_seq_orientation.py`` to reverse complement your reads if needed.
   * 4.  Chimera checking should be done first, followed by filtering chimeras out of the input reads, and these filtered sequences can then be clustered with ``pick_otus.py``.
   
An example step by step process for removing chimeras with USEARCH 6.1 starts as follows, using the seqs.fna file (output of ``split_libraries.py``) as the input sequence file: ::

    identify_chimeric_seqs.py -i seqs.fna -m usearch61 -o usearch_checked_chimeras/ -r gg_97_otus_4feb2011.fasta
    
Next, filter the input seqs.fna file by passing the chimeras.txt file created in the previous step and specifying that these should be removed with the -n option: ::
    
    filter_fasta.py -f seqs.fna -o seqs_chimeras_filtered.fna -s usearch_checked_chimeras/chimeras.txt -n
    
Finally, use the ``pick_otus.py`` to cluster these sequences using usearch61: ::

    pick_otus.py -m usearch61 -i seqs_chimeras_filtered.fna -o usearch61_picked_otus/
    
The remainder of the processing is the same as described in the main `tutorial <./tutorial.html#step-2-pick-representative-sequences-for-each-otu>`_.