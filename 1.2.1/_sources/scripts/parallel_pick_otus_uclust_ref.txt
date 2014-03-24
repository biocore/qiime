.. _parallel_pick_otus_uclust_ref:

.. index:: parallel_pick_otus_uclust_ref.py

*parallel_pick_otus_uclust_ref.py* -- Parallel pick otus using uclust_ref
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script works like the `pick_otus.py <./pick_otus.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_pick_otus_uclust_ref.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to input_fasta_fp
	-o, `-`-output_dir
		Path to store output files
	-r, `-`-refseqs_fp
		Full path to reference collection
	
	**[OPTIONAL]**
		
	-s, `-`-similarity
		Sequence similarity threshold [default: 0.97]
	-z, `-`-enable_rev_strand_match
		Enable reverse strand matching for uclust otu picking, will double the amount of memory used. [default: False]
	-A, `-`-optimal_uclust
		Pass the --optimal flag to uclust for uclust otu picking. [default: False]
	-E, `-`-exact_uclust
		Pass the --exact flag to uclust for uclust otu picking. [default: False]
	`-`-max_accepts
		Max_accepts value to uclust and uclust_ref [default: 20]
	`-`-max_rejects
		Max_rejects value to uclust and uclust_ref [default: 500]
	`-`-stepwords
		Stepwords value to uclust and uclust_ref [default: 20]
	`-`-word_length
		W value to uclust and uclust_ref [default: 12]
	`-`-uclust_stable_sort
		Deprecated: stable sort enabled by default, pass --uclust_suppress_stable_sort to disable [default: True]
	`-`-suppress_uclust_stable_sort
		Don't pass --stable-sort to uclust [default: False]
	-d, `-`-save_uc_files
		Enable preservation of intermediate uclust (.uc) files that are used to generate clusters via uclust. [default: True]
	-N, `-`-pick_otus_fp
		Full path to scripts/`pick_otus.py <./pick_otus.html>`_ [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`pick_otus.py <./pick_otus.html>`_]
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
		Path to `cluster_jobs.py <./cluster_jobs.html>`_ script  [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`start_parallel_jobs.py <./start_parallel_jobs.html>`_]
	-W, `-`-suppress_polling
		Suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		Full path to python executable [default: /opt/local/bin/python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 60]


**Output:**

The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an OTU identifier which is the reference sequence identifier, and the remaining fields correspond to sequence identifiers assigned to that OTU. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file.


**Example:**

Pick OTUs with uclust_ref by searching /home/qiime/inseqs.fasta against /home/qiime/refseqs.fasta and write the output to the /home/qiime/out/ directory.

::

	parallel_pick_otus_uclust_ref.py -i /home/qiime/inseqs.fasta -r /home/qiime/refseqs.fasta -o /home/qiime/out/


