.. _denoising_454_data:

Denoising of 454 Data Sets 
--------------------------

**QIIME script:** :file:`denoise.py`

The pyrosequencing technology employed by 454 sequencing machines produce characteristic sequencing errors, mostly imprecise signals for longer homopolymers runs. Most of the sequences contain no or only a few errors, but a few sequences contain enough errors to be classified as an additional rare OTU. The goal for the denoising procedure is to reduce the amount of erroneous OTUs and thus increasing the accuracy of the whole QIIME pipeline.

Currently, QIIME supports denoising using Chris Quince's PyroNoise (Quince et al., 2009), which needs to be installed separately.

The input to the denoising script is a textual representation of 454's .sff files, produced by 454's own tool sffinfo from the initial .sff file::

	sffinfo 454Reads.sff > 454Reads.sff.txt

**Input Arguments:**

.. note::

	-i SFF_FP, `-`-input_file=SFF_FP [REQUIRED]

		This is the path to the flowgram file (.sff.txt). 

	-o OUTPUT_DIR, `-`-output_dir=OUTPUT_DIR [Default: pyronoise_picked_otus/]

		This is the location where the resulting output should be written.

	-n NUM_CPUS, `-`-num_cpus=NUM_CPUS [Default: 1]

		This is the number of CPUs that should be used. 

	-s PRECISION, `-`-precision=PRECISION [Default: 15.0]

		This is the precision that should be used (passed to pyroNoise). 

	-c CUT_OFF, `-`-cut-off=CUT_OFF [Default: 0.05]

		This is the cut-off that should be used (passed to pyroNoise).

	-k, `-`-keep_intermediates [Default: False]

		If this parameter is passed, then the script will not delete intermediate PyroNoise files - which is useful for debugging.

**Output:**

The output of this script produces two files 1) a denoised FASTA set of cluster centroids and 2) an OTU mapping of flowgram identifiers to centroids.

**Examples:**

To denoise the flowgram sequences in :file:`454Reads.sff.txt`, you can use the 
following command::

	denoise.py -i 454Reads.sff.txt -o 454Reads_out/

which produces these two output files:

	* :file:`454Reads.fasta`: A denoised set of cluster centroids.
	* :file:`454Reads_otu.txt`: A mapping of flowgram identifiers to centroids

On a multi-processor machine pyroNoise can be run in parallel using mpirun, where the number of processors is passed to the script via -n, as shown by the following command::

	denoise.py -i 454Reads.sff.txt -o 454Reads_out/ -n 4

Since PyroNoise's steep computational requirement, you should limit the application to small data sets. Barcodes and primers are not taken into account here, and barcoded samples should be denoised in separate steps. See Chris's PyroNoise web site for details or use a combination of `split_libraries.py <./scripts/split_libraries.html>`_ and :file:`sfffile` (from the 454 software package) to separate the sequences into different sets.

References
^^^^^^^^^^

Quince, C., Lanzen, A., Curtis, T. P., Davenport, R. J., Hall, N., Head, I. M., et al. (2009). Accurate determination of microbial diversity from 454 pyrosequencing data. Nat Methods, 6(9), 639-641.