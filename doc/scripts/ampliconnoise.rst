.. _ampliconnoise:

.. index:: ampliconnoise.py

*ampliconnoise.py* -- Run AmpliconNoise
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


The steps performed by this script are:

1. Split input sff.txt file into one file per sample

2. Run scripts required for PyroNoise

3. Run scripts required for SeqNoise

4. Run scripts requred for Perseus (chimera removal)

5. Merge output files into one file similar to the output of `split_libraries.py <./split_libraries.html>`_

This script produces a denoised fasta sequence file such as:
>PC.355_41
CATGCTGCCTC...
...
>PC.636_23
CATGCTGCCTC...
...

Additionally, the intermediate results of the ampliconnoise pipeline are
written to an output directory.

Ampliconnoise must be installed and correctly configured, and parallelized
steps will be called with mpirun, not qiime's `start_parallel_jobs_torque.py <./start_parallel_jobs_torque.html>`_ script.




**Usage:** :file:`ampliconnoise.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		The mapping filepath
	-i, `-`-sff_filepath
		Sff.txt filepath
	-o, `-`-output_filepath
		The output file
	
	**[OPTIONAL]**
		
	-n, `-`-np
		Number of processes to use for mpi steps. Default: 2
	`-`-chimera_alpha
		Alpha value to Class.pl used for chimera removal  Default: -3.8228
	`-`-chimera_beta
		Beta value to Class.pl used for chimera removal  Default: 0.62
	`-`-seqnoise_resolution
		-s parameter passed to seqnoise. Default is 25.0 for titanium, 30.0 for flx
	-d, `-`-output_dir
		Directory for ampliconnoise intermediate results. Default is output_filepath_dir
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters. [if omitted, default values will be used]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: False]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	`-`-suppress_perseus
		Omit perseus from ampliconnoise workflow
	`-`-platform
		Sequencing technology, options are 'titanium','flx'. [default: flx]
	`-`-truncate_len
		Specify a truncation length for ampliconnoise.  Note that is this is not specified, the truncate length is chosen by the --platform option (220 for FLX, 400 for Titanium) [default: None]


**Output:**

a fasta file of sequences, with labels as:'>sample1_0' , '>sample1_1' ...


Run ampliconnoise, write output to anoise_out.fna, compatible with output of `split_libraries.py <./split_libraries.html>`_

::

	ampliconnoise.py -i Fasting_Example.sff.txt -m Fasting_Map.txt -o anoise_out.fna


