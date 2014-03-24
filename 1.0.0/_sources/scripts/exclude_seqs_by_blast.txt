.. _exclude_seqs_by_blast:

.. index:: exclude_seqs_by_blast.py

*exclude_seqs_by_blast.py* -- Exclude contaminated sequences using BLAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**



This code is designed to allow users of the QIIME workflow to conveniently exclude unwanted sequences from their data. This is mostly useful for excluding human sequences from runs to comply with Internal Review Board (IRB) requirements, but may also have other uses (e.g. perhaps excluding a major bacterial contaminant). Sequences from a run are searched against a user-specified subject database, where BLAST hits are screened by e-value and the percentage of the query that aligns to the sequence.

For human screening THINK CAREFULLY about the data set that you screen against. Are you excluding human non-coding sequences? What about mitochondrial sequences? This point is CRITICAL because submitting human sequences that are not IRB-approved is BAD.

(e.g. you would NOT want to just screen against just the coding sequences of the human genome as found in the KEGG .nuc files, for example)

One valid approach is to screen all putative 16S rRNA sequences against greengenes to ensure they are bacterial rather than human.

WARNING: You cannot use this script if there are spaces in the path to the database of fasta files because formatdb cannot handle these paths (this is a limitation of NCBI's tools and we have no control over it).



**Usage:** :file:`exclude_seqs_by_blast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-querydb
		The path to a FASTA file containing query sequences
	-d, `-`-subjectdb
		The path to a FASTA file to BLAST against
	-o, `-`-outputfilename
		 The base path/filename to save results. Sequences passing the screen, failing the screen, raw BLAST results and the log will be saved to your filename + '.screened', '.excluded', '.raw_blast_results', and '.sequence_exclusion_log' respectively.
	
	**[OPTIONAL]**
		
	-e, `-`-e_value
		The e-value cutoff for blast queries [DEFAULT: 1e-10]
	-p, `-`-percent_aligned
		The %% alignment cutoff for blast queries [DEFAULT: 0.97]
	`-`-blastmatroot
		Path to a folder containing blast matrices. [DEFAULT: None]
	`-`-working_dir
		Working dir for BLAST [DEFAULT: /tmp]
	-M, `-`-max_hits
		Max hits parameter for BLAST. CAUTION: Because filtering on alignment percentage occurs after BLAST, a max hits value of 1 in combination with an alignment percent filter could miss valid contaminants. [DEFAULT: 100]
	-W, `-`-word_size
		Word size to use for BLAST search [DEFAULT: 28]


**Output:**

Four output files are generated based on the supplied outputpath + unique suffixes:

1. "filename_prefix".screened: A FASTA file of sequences that did pass the screen (i.e. matched the database and passed all filters).

2. "filename_prefix".excluded: A FASTA file of sequences that did not pass the screen.

3. "filename_prefix".raw_blast_results: Contains the raw BLAST results from the screening.

4. "filename_prefix".sequence_exclusion_log: A log file summarizing the options used and results obtained.



**Examples:**

The following is a simple example, where the user can take a given FASTA file (i.e. resulting FASTA file from `pick_rep_set.py <./pick_rep_set.html>`_) and blast those sequences against a reference FASTA file containing the set of sequences which are considered contaminated:

::

	exclude_seqs_by_blast.py -i repr_set_seqs.fasta -d ref_seq_set.fna -o exclude_seqs/

Alternatively, if the user would like to change the percent of aligned sequence coverage ("-p") or the maximum E-value ("-e"), they can use the following command:

::

	exclude_seqs_by_blast.py -i repr_set_seqs.fasta -d ref_seq_set.fna -o exclude_seqs/ -p 0.95 -e 1e-10


