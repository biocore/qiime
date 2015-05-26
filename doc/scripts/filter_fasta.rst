.. _filter_fasta:

.. index:: filter_fasta.py

*filter_fasta.py* -- This script can be applied to remove sequences from a fasta or fastq file based on input criteria.
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
		An OTU map where sequences ids are those which should be retained.
	-s, `-`-seq_id_fp
		A list of sequence identifiers (or tab-delimited lines with a seq identifier in the first field) which should be retained.
	-b, `-`-biom_fp
		A biom file where otu identifiers should be retained.
	-a, `-`-subject_fasta_fp
		A fasta file where the seq ids should be retained.
	-p, `-`-seq_id_prefix
		Keep seqs where seq_id starts with this prefix.
	`-`-sample_id_fp
		Keep seqs where seq_id starts with a sample id listed in this file. Must be newline delimited and may not contain a header.
	-n, `-`-negate
		Discard passed seq ids rather than keep passed seq ids. [default: False]
	`-`-mapping_fp
		Mapping file path (for use with --valid_states). [default: None]
	`-`-valid_states
		Description of sample ids to retain (for use with --mapping_fp). [default: None]


**Output:**




**OTU map-based filtering:**

Keep all sequences that show up in an OTU map.

::

	filter_fasta.py -f inseqs.fasta -o otu_map_filtered_seqs.fasta -m otu_map.txt

**Chimeric sequence filtering:**

Discard all sequences that show up in chimera checking output. NOTE: It is very important to pass -n here as this tells the script to negate the request, or discard all sequences that are listed via -s. This is necessary to remove the identified chimeras from inseqs.fasta.

::

	filter_fasta.py -f inseqs.fasta -o non_chimeric_seqs.fasta -s chimeric_seqs.txt -n

**Sequence list filtering:**

Keep all sequences from as fasta file that are listed in a text file.

::

	filter_fasta.py -f inseqs.fasta -o list_filtered_seqs.fasta -s seqs_to_keep.txt

**biom-based filtering:**

Keep all sequences that are listed as observations in a biom file.

::

	filter_fasta.py -f inseqs.fastq -o biom_filtered_seqs.fastq -b otu_table.biom

**fastq filtering:**

Keep all sequences from a fastq file that are listed in a text file (note: file name must end with .fastq to support fastq filtering).

::

	filter_fasta.py -f inseqs.fastq -o list_filtered_seqs.fastq -s seqs_to_keep.txt

**sample id list filtering:**

Keep all sequences from a fasta file where the sample id portion of the sequence identifier is listed in a text file (sequence identifiers in fasta file must be in post-split libraries format: sampleID_seqID).

::

	filter_fasta.py -f sl_inseqs.fasta -o sample_id_list_filtered_seqs.fasta --sample_id_fp map.txt


