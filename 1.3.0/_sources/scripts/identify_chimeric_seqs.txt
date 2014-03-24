.. _identify_chimeric_seqs:

.. index:: identify_chimeric_seqs.py

*identify_chimeric_seqs.py* -- Identify chimeric sequences in input FASTA file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

A FASTA file of sequences, can be screened to remove chimeras (sequences generated due to the PCR amplification of multiple templates or parent sequences). QIIME currently includes a taxonomy-assignment-based approach, blast_fragments, for identifying sequences as chimeric and the ChimeraSlayer algorithm. 

1. Blast_fragments approach: 

The reference sequences (-r) and id-to-taxonomy map (-t) provided are the same format as those provided to `assign_taxonomy.py <./assign_taxonomy.html>`_. The reference sequences are in fasta format, and the id-to-taxonomy map contains tab-separated lines where the first field is a sequence identifier, and the second field is the taxonomy separated by semi-colons (e.g., Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium). The reference collection should be derived from a chimera-checked database (such as the full greengenes database), and filtered to contain only sequences at, for example, a maximum of 97% sequence identity.

2. ChimeraSlayer:

ChimeraSlayer uses BLAST to identify potential chimera parents and computes the optimal branching alignment of the query against two parents.
We suggest to use the pynast aligned representative sequences as input.



**Usage:** :file:`identify_chimeric_seqs.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-t, `-`-id_to_taxonomy_fp
		Path to tab-delimited file mapping sequences to assigned taxonomy. Each assigned taxonomy is provided as a comma-separated list. [default: None; REQUIRED when method is blast_fragments]
	-r, `-`-reference_seqs_fp
		Path to reference sequences (used to build a blast db when method blast_fragments). [default: None; REQUIRED when method blast_fragments if no blast_db is provided;]
	-a, `-`-aligned_reference_seqs_fp
		Path to (Py)Nast aligned reference sequences. REQUIRED when method ChimeraSlayer [default: None]
	-b, `-`-blast_db
		Database to blast against. Must provide either --blast_db or --reference_seqs_fp when method is blast_fragments [default: None]
	-m, `-`-chimera_detection_method
		Chimera detection method. Choices: blast_fragments or ChimeraSlayer. [default:ChimeraSlayer]
	-n, `-`-num_fragments
		Number of fragments to split sequences into (i.e., number of expected breakpoints + 1) [default: 3]
	-d, `-`-taxonomy_depth
		Number of taxonomic divisions to consider when comparing taxonomy assignments [default: 4]
	-e, `-`-max_e_value
		Max e-value to assign taxonomy [default: 1e-30]
	-R, `-`-min_div_ratio
		Min divergence ratio (passed to ChimeraSlayer). If set to None uses ChimeraSlayer default value.  [default: None]
	-k, `-`-keep_intermediates
		Keep intermediate files, useful for debugging  [default: False]
	-o, `-`-output_fp
		Path to store output [default: derived from input_seqs_fp]


**Output:**

The result of `identify_chimeric_seqs.py <./identify_chimeric_seqs.html>`_ is a text file that identifies which sequences are chimeric.


**blast_fragments example:**

For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in `assign_taxonomy.py <./assign_taxonomy.html>`_) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to the directory "identify_chimeras/" and using default parameters (e.g. chimera detection method ("-m blast_fragments"), number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")):

::

	identify_chimeric_seqs.py -i repr_set_seqs.fasta -t taxonomy_assignment.txt -r ref_seq_set.fna -o chimeric_seqs.txt

**ChimeraSlayer Example:**

Identify chimeric sequences using the ChimeraSlayer algorithm against a user provided reference data base. The input sequences need to be provided in aligned (Py)Nast format. The reference data base needs to be provided as aligned FASTA (-a). Note that the reference database needs to be the same that was used to build the alignment of the input sequences!

::

	identify_chimeric_seqs.py -m ChimeraSlayer -i repr_set_seqs_aligned.fasta -a ref_seq_set_aligned.fasta -o chimeric_seqs.txt


