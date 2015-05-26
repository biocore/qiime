.. _parallel_pick_otus_trie:

.. index:: parallel_pick_otus_trie.py

*parallel_pick_otus_trie.py* -- Parallel pick otus using a trie
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script performs like the `pick_otus.py <./pick_otus.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel. The script uses the first p bases of each read to sort all reads into separate buckets and then each buckets is processed separately. Note that in cases of amplicon sequencing we do not expect the buckets to be even sized, but rather a few buckets make up the majority of reads. Thus, not all combination of prefix length p and number of CPUS -O make sense. Good combinations for a small desktop multicore system would be -p 5 (default) and -O 4. For larger clusters, we suggest -p 10 and -O 20. Increasing -p to a value much larger than 10 will lead to lots of temporary files and many small jobs, so likely will not speed up the OTU picking. On the other hand, the max speed-up is bounded by the size of the largest buckets, so adding more cores will not always increase efficiency.


**Usage:** :file:`parallel_pick_otus_trie.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Full path to input_fasta_fp
	-o, `-`-output_dir
		Path to store output files
	
	**[OPTIONAL]**
		
	-p, `-`-prefix_length
		Prefix length used to split the input. Must be smaller than the shortest seq in input! [default: 5]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 4]
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

The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file.


**Example:**

Pick OTUs by building a trie out of $PWD/inseqs.fasta and write the output to the $PWD/trie_otus/ directory. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_pick_otus_trie.py -i $PWD/seqs.fna -o $PWD/trie_otus/

**Example:**

Pick OTUs by building a trie out of $PWD/inseqs.fasta and write the output to the $PWD/trie_otus/ directory. Split the input according to the first 10 bases of each read and process each set independently.

::

	parallel_pick_otus_trie.py -i $PWD/seqs.fna -o $PWD/trie_otus/ -p 10


