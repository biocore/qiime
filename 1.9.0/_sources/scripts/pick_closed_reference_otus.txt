.. _pick_closed_reference_otus:

.. index:: pick_closed_reference_otus.py

*pick_closed_reference_otus.py* -- Closed-reference OTU picking/Shotgun UniFrac workflow.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script picks OTUs using a closed reference and constructs an OTU table.
Taxonomy is assigned using a pre-defined taxonomy map of reference sequence OTU
to taxonomy. If full-length genomes are provided as the reference sequences,
this script applies the Shotgun UniFrac method.

**Note:** If most or all of your sequences are failing to hit the reference,
your sequences may be in the reverse orientation with respect to your reference
database. To address this, you should add the following line to your parameters
file (creating one, if necessary) and pass this file as -p:

pick_otus:enable_rev_strand_match True

Be aware that this doubles the amount of memory used.




**Usage:** :file:`pick_closed_reference_otus.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input sequences
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-r, `-`-reference_fp
		The reference sequences [default: /Users/jairideout/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta]. NOTE: If you do not pass -r to this script, you will be using QIIME's default reference sequences. In this case, QIIME will copy the corresponding reference tree to the output directory. This is the tree that should be used to perform phylogenetic diversity analyses (e.g., with `core_diversity_analyses.py <./core_diversity_analyses.html>`_).
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	-t, `-`-taxonomy_fp
		The taxonomy map [default: /Users/jairideout/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt]
	-s, `-`-assign_taxonomy
		Assign taxonomy to each sequence using `assign_taxonomy.py <./assign_taxonomy.html>`_ (this will override --taxonomy_fp, if provided) [default: False]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]
	`-`-suppress_taxonomy_assignment
		Skip the taxonomy assignment step, resulting in an OTU table without taxonomy (this will override --taxonomy_fp and --assign_taxonomy, if provided) [default: False]


**Output:**




Pick OTUs, assign taxonomy, and create an OTU table against a reference set of OTUs. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_w_tax/ -t $PWD/taxa.txt

Pick OTUs and create an OTU table against a reference set of OTUs without adding taxonomy assignments. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus/

Pick OTUs, assign taxonomy, and create an OTU table against a reference set of OTUs using usearch_ref. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_usearch/ -p $PWD/usearch_params.txt -t $PWD/taxa.txt

Pick OTUs using usearch_ref, assign taxonomy, and create an OTU table against a reference set of OTUs using usearch_ref. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).

::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_usearch_ref/ -p $PWD/usearch5.2_params.txt -t $PWD/taxa.txt

Pick OTUs, assign taxonomy, and create an OTU table against a reference set of OTUs using sortmerna. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented here as $PWD, but will generally look something like /home/ubuntu/my_analysis/). 

::

	pick_closed_reference_otus.py -i $PWD/seqs.fna -r $PWD/refseqs.fna -o $PWD/otus_sortmerna/ -p $PWD/sortmerna_params.txt -t $PWD/taxa.txt


