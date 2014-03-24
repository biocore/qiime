.. _denoiser_preprocess:

.. index:: denoiser_preprocess.py

*denoiser_preprocess.py* -- Run phase of denoiser algorithm: prefix clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The script `denoiser_preprocess.py <./denoiser_preprocess.html>`_ runs the first clustering phase
which groups reads based on common prefixes.


**Usage:** :file:`denoiser_preprocess.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_file
		Path to flowgram file [REQUIRED]
	
	**[OPTIONAL]**
		
	-f, `-`-fasta_file
		Path to fasta input file [default: None]
	-s, `-`-squeeze
		Use run-length encoding for prefix filtering [default: False]
	-l, `-`-log_file
		Path to log file [default: preprocess.log]
	-p, `-`-primer
		Primer sequence used for the amplification [default: CATGCTGCCTCCCGTAGGAGT]
	-o, `-`-output_dir
		Path to output directory [default: /tmp/]


**Output:**


prefix_dereplicated.sff.txt: human readable sff file containing the flowgram of the
                             cluster representative of each cluster.

prefix_dereplicated.fasta: Fasta file containing the cluster representative of each cluster.

prefix_mapping.txt: This file contains the actual clusters. The cluster centroid is given first,
                    the cluster members follw after the ':'.   



Run program on flowgrams in 454Reads.sff. Remove reads which are not in split_lib_filtered_seqs.fasta. 
Remove primer CATGCTGCCTCCCGTAGGAGT from reads before running phase I

::

	denoiser_preprocess.py -i Fasting_Example.sff.txt -f seqs.fna -p CATGCTGCCTCCCGTAGGAGT 


