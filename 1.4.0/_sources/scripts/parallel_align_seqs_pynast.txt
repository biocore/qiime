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
		Minimum sequence length to include in alignment [default: 150]
	-p, `-`-min_percent_id
		Minimum percent sequence identity to closest blast hit to include sequence in alignment [default: 75.0]
	-N, `-`-align_seqs_fp
		Full path to Qiime/scripts/`align_seqs.py <./align_seqs.html>`_ [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`align_seqs.py <./align_seqs.html>`_]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 1]
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
		Full path to python executable [default: /Library/Frameworks/Python.framework/Versions/2.7/bin/python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 60]
	-t, `-`-template_fp
		Filepath for template against [default: /data/greengenes_core_sets/core_set_aligned.fasta.imputed]


**Output:**

This results in a multiple sequence alignment (FASTA-formatted).


**Example:**

Align the input file (-i) against /home/qiime_user/pynast_test_template.fasta (-t) via 5 (-O) independent jobs and write the output (-o) to /home/qiime_user/out/:

::

	parallel_align_seqs_pynast.py -i /home/qiime_user/10_seq.fasta -O 5 -t /home/qiime_user/pynast_test_template.fasta -o /home/qiime_user/out/


