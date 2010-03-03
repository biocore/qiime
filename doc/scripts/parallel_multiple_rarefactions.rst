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
		input filepath, (the otu table) [REQUIRED]
	-o, `-`-output_path
		write output rarefied otu tables here makes dir if it doesn't exist [REQUIRED]
	-m, `-`-min
		min seqs/sample [REQUIRED]
	-x, `-`-max
		max seqs/sample (inclusive) [REQUIRED]
	-s, `-`-step
		levels: min, min+step... for level <= max [REQUIRED]
	
	**[OPTIONAL]**
		
	-n, `-`-num-reps
		num iterations at each seqs/sample level [default: 1]
	`-`-small_included
		samples containing fewer seqs than the rarefaction level are included in the output but not rarefied [default: False]
	`-`-lineages_included
		output rarefied otu tables will include taxonomic (lineage) information for each otu, if present in input otu table [default: False]
	-N, `-`-single_rarefaction_fp
		full path to scripts/`single_rarefaction.py <./single_rarefaction.html>`_ [default: /Users/Jesse/Qiime/scripts/`single_rarefaction.py <./single_rarefaction.html>`_]
	-P, `-`-poller_fp
		full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /Users/Jesse/Qiime/qiime/parallel/`poller.py <./poller.html>`_]
	-R, `-`-retain_temp_files
		retain temporary files after runs complete (useful for debugging) [default: False]
	-S, `-`-suppress_submit_jobs
		Only split input and write commands file - don't submit jobs [default: False]
	-T, `-`-poll_directly
		Poll directly for job completion rather than running poller as a separate job. If -T is specified this script will not return until all jobs have completed. [default: False]
	-U, `-`-cluster_jobs_fp
		path to `cluster_jobs.py <./cluster_jobs.html>`_ script  [default: /software/scripts/`cluster_jobs.py <./cluster_jobs.html>`_]
	-W, `-`-suppress_polling
		suppress polling of jobs and merging of results upon completion [default: False]
	-X, `-`-job_prefix
		job prefix [default: descriptive prefix + random chars]
	-Y, `-`-python_exe_fp
		full path to python executable [default: /usr/local/bin/python]
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 60]


**Output:**

The result of `parallel_multiple_rarefactions.py <./parallel_multiple_rarefactions.html>`_ consists of a number of files, which depend on the minimum/maximum number of sequences per samples, steps and iterations. The files have the same otu table format as the input otu_table.txt, and are named in the following way: rarefaction_100_0.txt, where "100" corresponds to the sequences per sample and "0" for the iteration.


**Example**

Build rarefied otu tables containing 100 (-m) to 2000 (-x) sequences in steps of 100 (-s) with 5 (-n) repetions per number of sequences, from otu_table.txt (-i). Write the output files to the rare directory (-o, will be created if it doesn't exist). The name of the output files will be of the form rare/rarefaction_<num_seqs>_<reptition_number>.txt

::

	parallel_multiple_rarefactions.py -o rare -m 100 -x 2000 -s 100 -n 5 -i otu_table.txt


