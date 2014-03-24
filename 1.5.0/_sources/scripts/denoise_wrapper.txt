.. _denoise_wrapper:

.. index:: denoise_wrapper.py

*denoise_wrapper.py* -- Denoise a flowgram file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script will denoise a flowgram file in .sff.txt format, which is the output of sffinfo.


**Usage:** :file:`denoise_wrapper.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_file
		Path to flowgram files (.sff.txt), comma separated
	-f, `-`-fasta_file
		Path to fasta file from `split_libraries.py <./split_libraries.html>`_
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to output directory [default: denoised_seqs/]
	-n, `-`-num_cpus
		Number of CPUs [default: 1]
	`-`-force_overwrite
		Overwrite files in output directory [default: False]
	-m, `-`-map_fname
		Name of mapping file, Has to contain field LinkerPrimerSequence. [REQUIRED unless --primer specified]
	-p, `-`-primer
		Primer sequence [REQUIRED unless --map_fname specified]
	`-`-titanium
		Select Titanium defaults for denoiser, otherwise use FLX defaults [default: False]


**Output:**

This script results in a OTU like mapping file along with a sequence file of denoised (FASTA-format). Note that the sequences coming from denoising are no real OTUs, and have to be sent to `pick_otus.py <./pick_otus.html>`_ if the users wishes to have a defined similarity threshold.


**Example:**

Denoise flowgrams in file 454Reads.sff.txt, discard flowgrams not in seqs.fna, and extract primer from map.txt:

::

	denoise_wrapper.py -i 454Reads.sff.txt -f seqs.fna -m map.txt

**Multi-core Example:**

Denoise flowgrams in file 454Reads.sff.txt using 2 cores on your machine in parallel:

::

	denoise_wrapper.py -n 2 -i 454Reads.sff.txt -f seqs.fna -m map.txt


