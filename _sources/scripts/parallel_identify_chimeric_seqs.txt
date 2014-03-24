.. _parallel_identify_chimeric_seqs:

.. index:: parallel_identify_chimeric_seqs.py

*parallel_identify_chimeric_seqs.py* -- Parallel chimera detection
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script works like the `identify_chimeric_seqs.py <./identify_chimeric_seqs.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_identify_chimeric_seqs.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-a, `-`-aligned_reference_seqs_fp
		Path to (Py)Nast aligned reference sequences. REQUIRED when method ChimeraSlayer [default: /Users/caporaso/data/greengenes_core_sets/core_set_aligned_imputed.fasta_11_8_07.no_dots]
	-t, `-`-id_to_taxonomy_fp
		Path to tab-delimited file mapping sequences to assigned taxonomy. Each assigned taxonomy is provided as a comma-separated list. [default: None; REQUIRED when method is blast_fragments]
	-r, `-`-reference_seqs_fp
		Path to reference sequences (used to build a blast db when method blast_fragments). [default: None; REQUIRED when method blast_fragments if no blast_db is provided;]
	-b, `-`-blast_db
		Database to blast against. Must provide either --blast_db or --reference_seqs_fp when method is blast_fragments [default: None]
	-m, `-`-chimera_detection_method
		Chimera detection method. Choices: blast_fragments or ChimeraSlayer. [default:ChimeraSlayer]
	-n, `-`-num_fragments
		Number of fragments to split sequences into (i.e., number of expected breakpoints + 1) [default: 3]
	-d, `-`-taxonomy_depth
		Number of taxonomic divisions to consider when comparing taxonomy assignments [default: 4]
	-e, `-`-max_e_value
		Max e-value to assign taxonomy [default: 1e-30]
	`-`-min_div_ratio
		Min divergence ratio (passed to ChimeraSlayer). If set to None uses ChimeraSlayer default value.  [default: None]
	-o, `-`-output_fp
		Path to store output [default: derived from input_seqs_fp]
	-O, `-`-jobs_to_start
		Number of jobs to start [default: 2]
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

The result of `parallel_identify_chimeric_seqs.py <./parallel_identify_chimeric_seqs.html>`_ is a text file that identifies which sequences are chimeric.


**blast_fragments example:**

For each sequence provided as input, the blast_fragments method splits the input sequence into n roughly-equal-sized, non-overlapping fragments, and assigns taxonomy to each fragment against a reference database. The BlastTaxonAssigner (implemented in `assign_taxonomy.py <./assign_taxonomy.html>`_) is used for this. The taxonomies of the fragments are compared with one another (at a default depth of 4), and if contradictory assignments are returned the sequence is identified as chimeric. For example, if an input sequence was split into 3 fragments, and the following taxon assignments were returned:

==========  ==========================================================
fragment1:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
fragment2:  Archaea;Euryarchaeota;Halobacteriales;uncultured
fragment3:  Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium
==========  ==========================================================

The sequence would be considered chimeric at a depth of 3 (Methanobacteriales vs. Halobacteriales), but non-chimeric at a depth of 2 (all Euryarchaeota).

blast_fragments begins with the assumption that a sequence is non-chimeric, and looks for evidence to the contrary. This is important when, for example, no taxonomy assignment can be made because no blast result is returned. If a sequence is split into three fragments, and only one returns a blast hit, that sequence would be considered non-chimeric. This is because there is no evidence (i.e., contradictory blast assignments) for the sequence being chimeric. This script can be run by the following command, where the resulting data is written to $PWD/blast_fragments_chimeric_seqs.txt and using default parameters (i.e., number of fragments ("-n 3"), taxonomy depth ("-d 4") and maximum E-value ("-e 1e-30")). ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_identify_chimeric_seqs.py -i $PWD/inseqs.fasta -t $PWD/id_to_tax.txt -r $PWD/refseqs.fasta -o $PWD/blast_fragments_chimeric_seqs.txt -m blast_fragments

**ChimeraSlayer Example:**

Identify chimeric sequences using the ChimeraSlayer algorithm against a user provided reference database. The input sequences need to be provided in aligned (Py)Nast format and the reference database needs to be provided as aligned FASTA (-a). Note that the reference database needs to be the same that was used to build the alignment of the input sequences! ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	parallel_identify_chimeric_seqs.py -i $PWD/inseqs_aligned.fasta -o $PWD/chimera_slayer_chimeric_seqs.txt


