.. _parallel_align_seqs_pynast:

.. index:: parallel_align_seqs_pynast.py

*parallel_align_seqs_pynast.py* -- Parallel sequence alignment using PyNAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

A wrapper for the `align_seqs.py <./align_seqs.html>`_ PyNAST option, intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_align_seqs_pynast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-o, `-`-output_dir
		Path to the output directory
	
	**[OPTIONAL]**
		
	-a, `-`-pairwise_alignment_method
		Method to use for pairwise alignments [default: uclust]
	-d, `-`-blast_db
		Database to blast against [default: created on-the-fly from template_alignment]
	-e, `-`-min_length
		Minimum sequence length to include in alignment [default: 75% of the median input sequence length]
	-p, `-`-min_percent_id
		Minimum percent sequence identity to closest blast hit to include sequence in alignment [default: 75.0]
	-N, `-`-align_seqs_fp
		Full path to Qiime/scripts/`align_seqs.py <./align_seqs.html>`_ [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`align_seqs.py <./align_seqs.html>`_]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]
	-P, `-`-poller_fp
		Full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`poller.py <./poller.html>`_]
	-R, `-`-retain_temp_files
		Retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		Path to cluster jobs script (defined in qiime_config)  [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		Full path to python executable [default: /usr/local/bin/python2.7]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 60]
	-t, `-`-template_fp
		Filepath for template against [default: /Users/jistombaugh/python_software/core_set_aligned.fasta.imputed]


**Output:**

This results in a multiple sequence alignment (FASTA-formatted).


**Example:**

Align the input file (-i) against using PyNAST and write the output (-o) to $PWD/pynast_aligned_seqs/. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_align_seqs_pynast.py -i $PWD/inseqs.fasta -o $PWD/pynast_aligned_seqs/


