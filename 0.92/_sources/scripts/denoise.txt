.. _denoise:

.. index:: denoise.py

*denoise.py* -- Denoise a flowgram file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script will denoise a flowgram file in  .sff.txt format, which is the output of sffinfo.


**Usage:** :file:`denoise.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_file
		Path to flowgram file (.sff.txt)
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to output directory [default: denoised_seqs/]
	-k, `-`-keep_intermediates
		Do not delete intermediate files -- useful for debugging [default: False]
	-c, `-`-cut-off
		Cut-off value (passed to pyroNoise) [default: 0.05]
	-s, `-`-precision
		Precision (passed to pyroNoise)[default: 15.0]
	-n, `-`-num_cpus
		Number of CPUs [default: 1]


**Output:**

This script results in a OTU mapping file along with a sequence file of denoised (FASTA-format). Note that the sequences coming from denoising are no real OTUs, and have to be sent to `pick_otus.py <./pick_otus.html>`_ if the users wishes to have a defined similarity threshold. 


**Example:**

Denoise flowgrams in file 454Reads.sff.txt:

::

	denoise.py -i 454Reads.sff.txt

**Multi-core Example:**

Denoise flowgrams in file 454Reads.sff.txt using 2 cores on your machine in parallel (requires mpirun):

::

	denoise.py -n 2 -i 454Reads.sff.txt


