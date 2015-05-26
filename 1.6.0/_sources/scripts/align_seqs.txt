.. _align_seqs:

.. index:: align_seqs.py

*align_seqs.py* -- Align sequences using a variety of alignment methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script aligns the sequences in a FASTA file to each other or to a template sequence alignment, depending on the method chosen. Currently, there are three methods which can be used by the user:

1. PyNAST (Caporaso et al., 2009) - The default alignment method is PyNAST, a python implementation of the NAST alignment algorithm.  The NAST algorithm aligns each provided sequence (the "candidate" sequence) to the best-matching sequence in a pre-aligned database of sequences (the "template" sequence).  Candidate sequences are not permitted to introduce new gap characters into the template database, so the algorithm introduces local mis-alignments to preserve the existing template sequence.

2. MUSCLE (Edgar, 2004) - MUSCLE is an alignment method which stands for MUltiple Sequence Comparison by Log-Expectation.

3. INFERNAL (Nawrocki, Kolbe, & Eddy, 2009) - Infernal ("INFERence of RNA ALignment") is for an alignment method for using RNA structure and sequence similarities.



**Usage:** :file:`align_seqs.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-t, `-`-template_fp
		Filepath for template against [default: /Users/caporaso/data/greengenes_core_sets/core_set_aligned_imputed.fasta_11_8_07.no_dots]
	-m, `-`-alignment_method
		Method for aligning sequences. Valid choices are: pynast, infernal, clustalw, muscle, infernal, mafft [default: pynast]
	-a, `-`-pairwise_alignment_method
		Method for performing pairwise alignment in PyNAST. Valid choices are muscle, pair_hmm, clustal, blast, uclust, mafft [default: uclust]
	-d, `-`-blast_db
		Database to blast against when -m pynast [default: created on-the-fly from template_alignment]
	`-`-muscle_max_memory
		Maximum memory allocation for the muscle alignment method (MB) [default: 80% of available memory, as detected by MUSCLE]
	-o, `-`-output_dir
		Path to store result file [default: <ALIGNMENT_METHOD>_aligned]
	-e, `-`-min_length
		Minimum sequence length to include in alignment [default: 75% of the median input sequence length]
	-p, `-`-min_percent_id
		Minimum percent sequence identity to closest blast hit to include sequence in alignment [default: 0.75]


**Output:**

All aligners will output a fasta file containing the alignment and log file in the directory specified by `-`-output_dir (default <alignment_method>_aligned). PyNAST additionally outputs a failures file, containing the sequences which failed to align. So the result of %prog will be up to three files, where the prefix of each file depends on the user supplied FASTA file:

1. "..._aligned.fasta" - This is a FASTA file containing all aligned sequences.

2. "..._failures.fasta" - This is a FASTA file containing all sequences which did not meet all the criteria specified. (PyNAST only)

3. "..._log.txt" - This is a log file containing information pertaining to the results obtained from a particular method (e.g. BLAST percent identity, etc.).


**Alignment with MUSCLE:**

One could also use the MUSCLE algorithm. The following command can be used to align sequences (i.e. the resulting FASTA file from `pick_rep_set.py <./pick_rep_set.html>`_), where the output is written to the directory "muscle_alignment/":

::

	align_seqs.py -i $PWD/unaligned.fna -m muscle -o $PWD/muscle_alignment/

**Alignment with PyNAST:**

The default alignment method is PyNAST, a python implementation of the NAST alignment algorithm. The NAST algorithm aligns each provided sequence (the "candidate" sequence) to the best-matching sequence in a pre-aligned database of sequences (the "template" sequence). Candidate sequences are not permitted to introduce new gap characters into the template database, so the algorithm introduces local mis-alignments to preserve the existing template sequence. The quality thresholds are the minimum requirements for matching between a candidate sequence and a template sequence. The set of matching template sequences will be searched for a match that meets these requirements, with preference given to the sequence length. By default, the minimum sequence length is 150 and the minimum percent id is 75%. The minimum sequence length is much too long for typical pyrosequencing reads, but was chosen for compatibility with the original NAST tool.

The following command can be used for aligning sequences using the PyNAST method, where we supply the program with a FASTA file of unaligned sequences (i.e. resulting FASTA file from `pick_rep_set.py <./pick_rep_set.html>`_, a FASTA file of pre-aligned sequences (this is the template file, which is typically the Greengenes core set - available from http://greengenes.lbl.gov/), and the results will be written to the directory "pynast_aligned/":

::

	align_seqs.py -i $PWD/unaligned.fna -t $PWD/core_set_aligned.fasta.imputed -o $PWD/pynast_aligned_defaults/

Alternatively, one could change the minimum sequence length ("-e") requirement and minimum sequence identity ("-p"), using the following command:

::

	align_seqs.py -i $PWD/unaligned.fna -t core_set_aligned.fasta.imputed -o $PWD/pynast_aligned/ -e 500 -p 95.0

**Alignment with Infernal:**

An alternative alignment method is to use Infernal. Infernal is similar to the PyNAST method, in that you supply a template alignment, although Infernal has several distinct differences. Infernal takes a multiple sequence alignment with a corresponding secondary structure annotation. This input file must be in Stockholm alignment format. There is a fairly good description of the Stockholm format rules at: http://en.wikipedia.org/wiki/Stockholm_format. Infernal will use the sequence and secondary structural information to align the candidate sequences to the full reference alignment. Similar to PyNAST, Infernal will not allow for gaps to be inserted into the reference alignment. Using Infernal is slower than other methods, and therefore is best used with sequences that do not align well using PyNAST.

The following command can be used for aligning sequences using the Infernal method, where we supply the program with a FASTA file of unaligned sequences, a STOCKHOLM file of pre-aligned sequences and secondary structure (this is the template file - an example file can be obtained from: http://bmf.colorado.edu/QIIME/seed.16s.reference_model.sto.zip), and the results will be written to the directory "infernal_aligned/":

::

	align_seqs.py -m infernal -i $PWD/unaligned.fna -t $PWD/seed.16s.reference_model.sto -o $PWD/infernal_aligned/


