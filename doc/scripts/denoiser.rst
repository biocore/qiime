.. _denoiser:

.. index:: denoiser.py

*denoiser.py* -- Remove noise from  454 sequencing data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The denoiser removes sequencing noise characteristic to pyrosequencing by flowgram clustering. For a detailed explanation of the underlying algorithm see (Reeder and Knight, Nature Methods 7(9), 2010).


**Usage:** :file:`denoiser.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_files
		Path to flowgram files (.sff.txt), comma separated
	
	**[OPTIONAL]**
		
	-f, `-`-fasta_fp
		Path to fasta input file. Reads not in the fasta file are filtered out before denoising. File format is as produced by `split_libraries.py <./split_libraries.html>`_ [default: None]
	-o, `-`-output_dir
		Path to output directory [default: random dir in ./]
	-c, `-`-cluster
		Use cluster/multiple CPUs for flowgram alignments [default: False]
	-p, `-`-preprocess_fp
		Do not do preprocessing (phase I),instead use already preprocessed data in PREPROCESS_FP
	`-`-checkpoint_fp
		Resume denoising from checkpoint. Be careful when changing parameters for a resumed run. Requires -p option.  [default: None]
	-s, `-`-squeeze
		Use run-length encoding for prefix filtering in phase I [default: False]
	-S, `-`-split
		Split input into per library sets and denoise separately [default: False]
	`-`-force
		Force overwrite of existing directory [default: False]
	`-`-primer
		Primer sequence [default: CATGCTGCCTCCCGTAGGAGT]
	-n, `-`-num_cpus
		Number of cpus, requires -c [default: 1]
	-m, `-`-max_num_iterations
		Maximal number of iterations in phase II. None means unlimited iterations [default: None]
	-b, `-`-bail_out
		Stop clustering in phase II with clusters smaller or equal than BAILde [default: 1]
	`-`-percent_id
		Sequence similarity clustering threshold, expressed as a fraction between 0 and 1 [default: 0.97]
	`-`-low_cut_off
		Low clustering threshold for phase II [default: 3.75]
	`-`-high_cut_off
		High clustering threshold for phase III [default: 4.5]
	`-`-low_memory
		Use slower, low memory method [default: False]
	-e, `-`-error_profile
		Path to error profile [default= <qiime-install-path>/qiime/support_files/denoiser/Data/FLX_error_profile.dat]
	`-`-titanium
		Shortcut for -e <qiime-install-path>/qiime/support_files/denoiser/Data//Titanium_error_profile.dat --low_cut_off=4 --high_cut_off=5 . Warning: overwrites all previous cut-off values [DEFAULT: False]


**Output:**



centroids.fasta: The cluster representatives of each cluster

singletons.fasta: contains all unclustered reads

denoiser_mapping.txt: This file contains the actual clusters. The cluster centroid is given first,
                    the cluster members follow after the ':'.

checkpoints/ : directory with checkpoints

Note that the centroids and singleton files are disjoint. For most downstream analyses one wants to cat the two files.



Run denoiser on flowgrams in 454Reads.sff.txt with read-to-barcode mapping in seqs.fna,
put results into Outdir, log progress in Outdir/denoiser.log

::

	denoiser.py -i 454Reads.sff.txt -f seqs.fna -v -o Outdir

**Multiple sff.txt files:**

Run denoiser on two flowgram files in 454Reads_1.sff.txt and 454Reads_2.sff.txt
with read-to-barcode mapping in seqs.fna, put results into Outdir,
log progress in Outdir/denoiser.log

::

	denoiser.py -i 454Reads_1.sff.txt,454Reads_2.sff.txt -f seqs.fna -v -o Outdir

**Denoise multiple library separately:**

Run denoiser on flowgrams in 454Reads.sff.txt with read-to-barcode mapping in seqs.fna,
split input files into libraries and process each library separately,
put results into Outdir, log progress in Outdir/denoiser.log

::

	denoiser.py -S -i 454Reads.sff.txt -f seqs.fna -v -o Outdir

**Resuming a failed run:**

Resume a previous denoiser run from breakpoint stored in Outdir_from_failed_run/checkpoints/checkpoint100.pickle.
The checkpoint option requires the -p or --preprocess option, which usually can be set to the output dir of the failed run.
All other arguments must be identical to the failed run.

::

	denoiser.py -i 454Reads.sff.txt -f seqs.fna -v -o Outdir_resumed -p Outdir_from_failed_run --checkpoint Outdir_from_failed_run/checkpoints/checkpoint100.pickle


