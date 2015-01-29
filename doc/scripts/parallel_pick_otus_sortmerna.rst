.. _parallel_pick_otus_sortmerna:

.. index:: parallel_pick_otus_sortmerna.py

*parallel_pick_otus_sortmerna.py* -- Parallel pick otus using SortMeRNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script works like the `pick_otus.py <./pick_otus.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_pick_otus_sortmerna.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to input fasta file.
	-o, `-`-output_dir
		Directory where output should be written.
	-r, `-`-refseqs_fp
		Path to reference fasta file.
	
	**[OPTIONAL]**
		
	-s, `-`-similarity
		Sequence similarity threshold [default: 0.97]
	`-`-sortmerna_db
		Pre-existing database to search against when using -m sortmerna [default: None]
	`-`-sortmerna_e_value
		Maximum E-value when clustering [default = 1]
	`-`-sortmerna_coverage
		Mininum percent query coverage (of an alignment) to consider a hit, expressed as a fraction between 0 and 1 [default: 0.97]
	`-`-sortmerna_tabular
		Output alignments in the Blast-like tabular format with two additional columns including the CIGAR string and the percent query coverage [default: False]
	`-`-sortmerna_best_N_alignments
		Must be set together with --sortmerna_tabular. This option specifies how many alignments per read will be written [default: 1]
	`-`-sortmerna_max_pos
		The maximum number of positions per seed to store  in the indexed database [default: 10000]
	`-`-threads
		Specify the number of threads to use per job. Use --jobs_to_start to specify the number of jobs.[default: 1]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 1]
	-R, `-`-retain_temp_files
		Retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		Path to cluster jobs script (defined in qiime_config)  [default: `start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 1]


**Output:**

The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an OTU identifier which is the reference sequence identifier, and the remaining fields correspond to sequence identifiers assigned to that OTU. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file.


**Example:**

Pick OTUs by searching $PWD/inseqs.fasta against $PWD/refseqs.fasta with reference-based sortmerna and write the output to the $PWD/smr_otus/ directory. This is a closed-reference OTU picking process. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_pick_otus_sortmerna.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/smr_otus/


