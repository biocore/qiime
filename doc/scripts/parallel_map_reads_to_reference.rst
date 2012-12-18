.. _parallel_map_reads_to_reference:

.. index:: parallel_map_reads_to_reference.py

*parallel_map_reads_to_reference.py* -- 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**




**Usage:** :file:`parallel_map_reads_to_reference.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_seqs_filepath
		Path to input sequences file
	-o, `-`-output_dir
		Directory to store results
	-r, `-`-refseqs_fp
		Path to reference sequences
	
	**[OPTIONAL]**
		
	-t, `-`-observation_metadata_fp
		Path to observation metadata (e.g., taxonomy, EC, etc) [default: None]
	-m, `-`-assignment_method
		Method for picking OTUs.  Valid choices are: usearch blat bwa-short. [default: usearch]
	-e, `-`-evalue
		Max e-value to consider a match [default: 1e-10]
	-s, `-`-min_percent_id
		Min percent id to consider a match [default: 0.75]
	`-`-max_diff
		MaxDiff to consider a match (applicable for -m bwa) -- see the aln section of "man bwa" for details [default (defined by bwa): 0.04]
	`-`-queryalnfract
		Min percent of the query seq that must match to consider a match (usearch only) [default: 0.35]
	`-`-targetalnfract
		Min percent of the target/reference seq that must match to consider a match (usearch only) [default: 0.0]
	`-`-max_accepts
		Max_accepts value (usearch only) [default: 1]
	`-`-max_rejects
		Max_rejects value to (usearch only) [default: 32]
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




::

	parallel_map_reads_to_reference.py -i $PWD/query_nt.fasta -r $PWD/refseqs_pr.fasta -o $PWD/usearch_mapped

::

	parallel_map_reads_to_reference.py -i $PWD/query_nt.fasta -r $PWD/refseqs_pr.fasta -o $PWD/blat_mapped -m blat

::

	parallel_map_reads_to_reference.py -i $PWD/query_nt.fasta -r $PWD/refseqs_nt.fasta -o $PWD/bwa-short_mapped -m bwa-short


