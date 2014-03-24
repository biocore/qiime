.. _parallel_pick_otus_blast:

.. index:: parallel_pick_otus_blast.py

*parallel_pick_otus_blast.py* -- Parallel pick otus using BLAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `pick_otus.py <./pick_otus.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_pick_otus_blast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to input_fasta_fp
	-o, `-`-output_dir
		Path to store output files
	
	**[OPTIONAL]**
		
	-e, `-`-max_e_value
		Max E-value [default: 1e-10]
	-s, `-`-similarity
		Sequence similarity threshold [default: 0.97]
	-r, `-`-refseqs_fp
		Full path to template alignment [default: None]
	-b, `-`-blast_db
		Database to blast against [default: None]
	`-`-min_aligned_percent
		Minimum percent of query sequence that can be aligned to consider a hit (BLAST OTU picker only) [default: 0.5]
	-N, `-`-pick_otus_fp
		Full path to scripts/`pick_otus.py <./pick_otus.html>`_ [default: /Users/caporaso/code/Qiime/scripts/`pick_otus.py <./pick_otus.html>`_]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]
	-P, `-`-poller_fp
		Full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /Users/caporaso/code/Qiime/scripts/`poller.py <./poller.html>`_]
	-R, `-`-retain_temp_files
		Retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		Path to cluster jobs script (defined in qiime_config)  [default: /Users/caporaso/code/Qiime/scripts/`start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		Full path to python executable [default: python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 6]


**Output:**

The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file.


**Example:**

Pick OTUs by blasting /home/qiime_user/inseqs.fasta against /home/qiime_user/refseqs.fasta and write the output to the /home/qiime_user/out/ directory.

::

	parallel_pick_otus_blast.py -i /home/qiime_user/inseqs.fasta -r /home/qiime_user/refseqs.fasta -o /home/qiime_user/out/


