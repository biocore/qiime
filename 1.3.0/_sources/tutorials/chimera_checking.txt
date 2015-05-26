.. _chimera_checking:

=====================================
Chimera checking sequences with QIIME
=====================================

This document describes how to chimera check sequences with QIIME, and subsequently exclude chimeric sequences when building your tree and OTU table so they will not be included in downstream analyses. (Note that because we recommend applying chimera checking after OTU and representative sequence picking, for the purposes of this document `chimeric sequence` should be considered synonymous to `chimeric OTU`.)

We recommend using QIIME's ChimeraSlayer wrapper to chimera check your sequences. As ChimeraSlayer requires aligned sequences, you should apply chimera checking after you've aligned your sequences with PyNAST, but before you apply the lanemask filtering to your alignment.

Applying ChimeraSlayer is done as follows::

	identify_chimeric_seqs.py -m ChimeraSlayer -i rep_set_aligned.fasta -a reference_set_aligned.fasta -o chimeric_seqs.txt

Before building your phylogenetic tree, you should now remove chimeric sequences from your alignment using your chimeric sequence list (``chimeric_seqs.txt`` in this example) and the alignment to as follows::

	filter_fasta.py -f rep_set_aligned.fasta -o non_chimeric_rep_set_aligned.fasta -s chimeric_seqs.txt -n
	
Don't forget to pass the ``-n`` parameter to ``filter_fasta.py`` -- this tells ``filter_fasta.py`` to discard the sequences you've passed via ``-s``, rather than to keep only those sequences. Then pass the resulting subalignment (``non_chimeric_rep_set_aligned.fasta``) in the downstream steps ``filter_alignment.py`` prior to tree-building.

You'll also need to exclude the chimeric sequences at the ``make_otu_table.py`` step. To do that, you pass the chimeric sequence list to ``make_otu_table.py`` via the ``-e`` parameter::

	make_otu_table.py -i otu_map.txt -o otu_table.txt -e chimeric_seqs.txt -t taxonomy.txt

This will cause the chimeric sequences to be excluded from all downstream analysis that you'll be performing, as the OTU table is always the input provided (for example to beta diversity, alpha rarefaction, etc.).
