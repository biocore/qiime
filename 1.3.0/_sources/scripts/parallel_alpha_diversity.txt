.. _parallel_alpha_diversity:

.. index:: parallel_alpha_diversity.py

*parallel_alpha_diversity.py* -- Parallel alpha diversity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `alpha_diversity.py <./alpha_diversity.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_alpha_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path, must be directory [REQUIRED]
	-o, `-`-output_path
		Output path, must be directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-t, `-`-tree_path
		Path to newick tree file, required for phylogenetic metrics [default: None]
	-N, `-`-alpha_diversity_fp
		Full path to scripts/`alpha_diversity.py <./alpha_diversity.html>`_ [default: /Users/caporaso/code/Qiime/scripts/`alpha_diversity.py <./alpha_diversity.html>`_]
	-m, `-`-metrics
		Metrics to use, comma delimited
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
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]


**Output:**

The resulting output will be the same number of files as supplied by the user. The resulting files are tab-delimited text files, where the columns correspond to alpha diversity metrics and the rows correspond to samples and their calculated diversity measurements. 


**Example:**

Apply the observed_species, chao1, PD_whole_tree metrics (-m) to all otu tables in /home/qiime_user/rare/ (-i) and write the resulting output files to /home/qiime_user/out/ (-o, will be created if it doesn't exist). Use the tree file rep_set.tre (-t) when necessary.

::

	parallel_alpha_diversity.py -i /home/qiime_user/rare/ -o /home/qiime_user/out -m observed_species,chao1,PD_whole_tree -t /home/qiime_user/rep_set.tre


