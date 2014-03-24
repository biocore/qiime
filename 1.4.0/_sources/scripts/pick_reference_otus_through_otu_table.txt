.. _pick_reference_otus_through_otu_table:

.. index:: pick_reference_otus_through_otu_table.py

*pick_reference_otus_through_otu_table.py* -- Reference OTU picking/Shotgun UniFrac workflow.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script picks OTUs using a reference-based method and constructs an OTU table. Taxonomy is assigned using a pre-defined taxonomy map of reference sequence OTU to taxonomy. If full-length genomes are provided as the reference sequences, this script applies the Shotgun UniFrac method.


**Usage:** :file:`pick_reference_otus_through_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fp
		The input sequences
	-r, `-`-reference_fp
		The reference sequences
	-o, `-`-output_dir
		The output directory
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters . [if omitted, default values will be used]
	-t, `-`-taxonomy_fp
		The taxonomy map [default: None]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-a, `-`-parallel
		Run in parallel where available [default: False]
	-O, `-`-jobs_to_start
		Number of jobs to start. NOTE: you must also pass -a to run in parallel, this defines the number of jobs to be started if and only if -a is passed [default: 1]


**Output:**




Pick OTUs, assign taxonomy, and create an OTU table against a reference set of OTUs.

::

	pick_reference_otus_through_otu_table.py -i inseqs.fasta -r refseqs.fasta -o out -p qiime_parameters.txt -t taxa.txt

Pick OTUs and create an OTU table against a reference set of OTUs without adding taxonomy assignments.

::

	pick_reference_otus_through_otu_table.py -i inseqs.fasta -r refseqs.fasta -o out -p qiime_parameters.txt


