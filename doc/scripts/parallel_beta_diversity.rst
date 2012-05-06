.. _parallel_beta_diversity:

.. index:: parallel_beta_diversity.py

*parallel_beta_diversity.py* -- Parallel beta diversity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `beta_diversity.py <./beta_diversity.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_beta_diversity.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path, must be directory [REQUIRED]
	-o, `-`-output_path
		Output path, must be directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-m, `-`-metrics
		Beta-diversity metric(s) to use. A comma-separated list should be provided when multiple metrics are specified. [default: unweighted_unifrac,weighted_unifrac]
	-t, `-`-tree_path
		Path to newick tree file, required for phylogenetic metrics [default: None]
	-N, `-`-beta_diversity_fp
		Full path to scripts/`beta_diversity.py <./beta_diversity.html>`_ [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`beta_diversity.py <./beta_diversity.html>`_]
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
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]
	-f, `-`-full_tree
		By default, each job removes calls _fast_unifrac_setup to remove unused parts of the tree. pass -f if you already have a minimal tree, and this script will run faster


**Output:**

The output of %prog is a folder containing text files, each a distance matrix between samples.


**Apply beta_diversity.py in parallel:**

Apply the unweighted_unifrac and weighted_unifrac metrics (modify with -m) to all otu tables in rarefied_otu_tables (-i) and write the resulting output files to bdiv/ (-o, will be created if it doesn't exist). Use the rep_set.tre (-t) to compute phylogenetic diversity metrics. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_beta_diversity.py -i $PWD/rarefied_otu_tables/ -o $PWD/bdiv/ -t $PWD/rep_set.tre


