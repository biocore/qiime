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
	-m, `-`-metrics
		Metrics to use, comma delimited
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
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]


**Output:**

The resulting output will be the same number of files as supplied by the user. The resulting files are tab-delimited text files, where the columns correspond to alpha diversity metrics and the rows correspond to samples and their calculated diversity measurements. 


**Example:**

Apply the observed_species, chao1, PD_whole_tree metrics (-m) to all otu tables in rarefied_otu_tables/ (-i) and write the resulting output files to adiv/ (-o, will be created if it doesn't exist). Use the rep_set.tre (-t) to compute phylogenetic diversity metrics. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_alpha_diversity.py -i $PWD/rarefied_otu_tables -o $PWD/adiv -m observed_species,chao1,PD_whole_tree -t $PWD/rep_set.tre


