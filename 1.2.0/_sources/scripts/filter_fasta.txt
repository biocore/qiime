.. _filter_fasta:

.. index:: filter_fasta.py

*filter_fasta.py* -- This script can be applied to remove sequences from a fasta file based on input criteria.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`filter_fasta.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-input_fasta_fp
		Path to the input fasta file
	-o, `-`-output_fasta_fp
		The output fasta filepath
	
	**[OPTIONAL]**
		
	-m, `-`-otu_map
		An OTU map where sequences ids are those which should be retained
	-s, `-`-seq_id_fp
		A list of sequence identifiers (or tab-delimited lines with a seq identifier in the first field) which should be retained
	-a, `-`-subject_fasta_fp
		A fasta file where the seq ids should be retained.
	-p, `-`-seq_id_prefix
		Keep seqs where seq_id starts with this prefix
	-n, `-`-negate
		Discard passed seq ids rather than keep passed seq ids [default: False]


**Output:**




**Keep all sequences that show up in an OTU map:**

::

	filter_fasta.py -f inseqs.fasta -o filtered_seqs.fasta -m uclust_ref_otus.txt

**Discard all sequences that show up in chimera checking output. NOTE: It is very important to pass -n here as this tells the script to negate the request, or discard all sequences that are listed via -s. This is necessary to remove the identified chimeras from inseqs.fasta:**

::

	filter_fasta.py -f inseqs.fasta -o non_chimeric_seqs.fasta -s chimeric_seqs.txt -n

**Keep all sequences listed in a text file:**

::

	filter_fasta.py -f inseqs.fasta -o filtered_seqs.fasta -s seqs_to_keep.txt


