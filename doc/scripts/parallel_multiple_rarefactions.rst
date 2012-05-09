.. _parallel_multiple_rarefactions:

.. index:: parallel_multiple_rarefactions.py

*parallel_multiple_rarefactions.py* -- Parallel multiple file rarefaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `multiple_rarefactions.py <./multiple_rarefactions.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_multiple_rarefactions.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input filepath, (the otu table) [REQUIRED]
	-o, `-`-output_path
		Write output rarefied otu tables here makes dir if it doesn't exist [REQUIRED]
	-m, `-`-min
		Min seqs/sample [REQUIRED]
	-x, `-`-max
		Max seqs/sample (inclusive) [REQUIRED]
	
	**[OPTIONAL]**
		
	-n, `-`-num-reps
		Num iterations at each seqs/sample level [default: 10]
	`-`-lineages_included
		Deprecated: lineages are now included by default. Pass --supress_lineages_included to prevent output OTU tables from including taxonomic (lineage) information for each OTU. Note: this will only work if lineage information is in the input OTU table.
	`-`-suppress_lineages_included
		Exclude taxonomic (lineage) information for each OTU.
	-N, `-`-single_rarefaction_fp
		Full path to scripts/`single_rarefaction.py <./single_rarefaction.html>`_ [default: /Users/jistombaugh/Dropbox/Qiime_work/scripts/`single_rarefaction.py <./single_rarefaction.html>`_]
	-s, `-`-step
		Levels: min, min+step... for level <= max [default: 1]
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


**Output:**

The result of `parallel_multiple_rarefactions.py <./parallel_multiple_rarefactions.html>`_ consists of a number of files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.biom, and are named in the following way: rarefaction_100_0.txt, where "100" corresponds to the sequences per sample and "0" for the iteration.


**OTU tables of different depths:**

Build rarefied otu tables containing 10 (-m) to 140 (-x) sequences in steps of 10 (-s) with 2 (-n) repetions per number of sequences, from otu_table.biom (-i). Write the output files to the rarefied_otu_tables directory (-o, will be created if it doesn't exist). The name of the output files will be of the form rarefaction_<num_seqs>_<reptition_number>.biom. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_multiple_rarefactions.py -o $PWD/rarefied_otu_tables/ -m 10 -x 140 -s 10 -n 2 -i $PWD/otu_table.biom

**OTU tables of the same depth:**

Build 8 rarefied otu tables each containing exactly 100 sequences per sample (even depth rarefaction). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_multiple_rarefactions.py -o $PWD/even_otu_tables/ -m 100 -x 100 -n 8 -i $PWD/otu_table.biom


