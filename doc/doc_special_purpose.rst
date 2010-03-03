.. _doc_special_purpose:

.. toctree::
   :maxdepth: 1

=========================
Special-Purpose Tutorials
=========================

Parallel Runs 
-------------

QIIME supports running several of its slower steps in parallel in a cluster (or other multiple processor/core) environment. Currently, these include:

	* Assignment of taxonomy with BLAST, via :file:`Qiime/qiime/parallel/assign_taxonomy_blast.py`
	* Assignment of taxonomy with RDP, via :file:`Qiime/qiime/parallel/assign_taxonomy_rdp.py`
	* Sequence alignment with PyNAST, via :file:`Qiime/qiime/parallel/align_seqs_pynast.py`

Writing a cluster_jobs Script Specific to your Cluster Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To make QIIME parallelization useful in different computing environments users are required to provide a script which can start jobs on their system, referred to here as a 'cluster_jobs' script. The cluster_jobs script takes as its two parameterts:

	1. A single file which lists the commands to be run (referred to as a 'jobs_list' file), with one command per line
	2. A string to use as a prefix when constructing unique job identifiers.

The lines in an example jobs_list file might be:

.. note::

	* python pick_otus.py -i inseqs_file1.fasta 
	* python pick_otus.py -i inseqs_file2.fasta 
	* python pick_otus.py -i inseqs_file3.fasta 

If passed to your cluster_jobs script, this should start three separate jobs corresponding to each of the commands.

The call to the cluster_jobs script in QIIME's parallel scripts looks like the following::

	CLUSTER_JOBS_FP -ms job_list.txt JOB_ID

where CLUSTER_JOBS_FP is the path to your cluster_jobs script and is passed to the parallel scripts via the -U parameter. JOB_ID is intended to be used as a prefix by the cluster_jobs script when creating a unique identifier for each job (and will be passed to the parallel scripts via -X). The same JOB_ID is also used by the QIIME parallel scripts when creating names for temporary files and directories. The -ms indicates that the job files should be made (-m) and submitted (-s).

Once you have written a cluster_jobs script for your specific environment that can be called via the above interface, running QIIME jobs in parallel should be straight-forward. The parallel variants of the scripts use the same parameters as the serial versions of the scripts, with some additional options in the parallel scripts. Options -N through -Z (capital N through capital Z) are reserved in QIIME for parallel scripts, and in most cases the defaults can be defined in your :file:`qiime_config` file.

Example Run of PyNAST in Parallel 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command will start a parallel run of PyNAST, which uses the same interface as the `align_seqs.py <./scripts/align_seqs.html>`_ script, where the results are written the an output directory "parallel_aligned_seqs/"::

	parallel_align_seqs_pynast.py -i repr_set_seqs.fasta -t /ref_set_seqs.fasta -o /home/caporaso/out 

The important thing to note is that this command is that same that would be used to call serial (single processor) PyNAST, except that instead of calling `parallel_align_seqs_pynast.py <./scripts/parallel_align_seqs_pynast.html>`_, you would call `align_seqs.py <./scripts/align_seqs.html>`_ to start the run on a single processor. The output from this parallel run is the same as the output would be from the serial run. 

Details of the Parallelization 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section provides some information on details of the parallelization which are hidden from the user, but provided for users who are interested in what is happening behind-the-scenes.

The parallelization works as follows. First, the input file (-i) is split into JOBS_TO_START (-O) different roughly equal-sized files. The serial version of the script -- `align_seqs.py <./scripts/align_seqs.html>`_ -- is then called on each of these split files as a separate job. Each of these jobs therefore writes its own output files (alignment, log, and failure files). One additional job, the poller, is started to monitor each of the jobs via their output files. When all expected output files exist, the poller will merge the individual output files and clean up any temporary files including the output files created by each of the individual runs. Cleaning up temporary files can be suppressed by passing -R, which is useful for debugging. Bypassing the polling system all-together can be achieved by passing -W.

The qiime_config File 
---------------------

First things first: you should not edit or remove :file:`Qiime/qiime_config`. 

Some QIIME scripts, at this stage primarily the parallel scripts, read default values from a :file:`qiime_config` file. The default location of this file is in your top-level QIIME directory (:file:`Qiime/qiime_config`). QIIME scripts pull default values from this file which are system-specific, such as paths to executable files, and these can be overwritten for convenience. The recommended procedure for overwriting these defaults is to copy the :file:`qiime_config` file to either :file:`~/.qiime_config` or a location specified by the environment variable $QIIME_CONFIG_FP.

The Qiime configuration values should only be modified in these copies of the :file:`qiime_config` file, as changes to the :file:`Qiime/qiime_config` version may be overwritten in future QIIME updates.

When defaults are loaded, all three locations are checked in order of precedence. Lowest precedence is given to the :file:`Qiime/qiime_config` file, as these are defaults defined by the QIIME development team and are likely not relevant to other users' environments. Higher precedence is given to the file specified by $QIIME_CONFIG_FP, and this is envisioned to be used for defining system-wide defaults. Finally, highest precedence is given to :file:`~/.qiime_config`, so users have the ability to overwrite defaults defined elsewhere to have maximum control over their environment (e.g., if testing an experimental version of their :file:`cluster_jobs` script). Note that these values are defaults: the scripts typically allow overwriting of these values via their command line interfaces.

Note that users can have up to three separate :file:`qiime_config` files, and one is provided by default with QIIME. At least one :file:`qiime_config` file must be present in one of the three locations, or scripts that rely on :file:`qiime_config` file will raise an error. Not all values need to be defined in all :file:`qiime_config` files, but all values must be defined at least once. This is one more reason why you should not edit or remove :file:`Qiime/qiime_config`: when new values are added in the future they will be defined in Qiime's default copy, but not in your local copies.

There is a script that prints the current :file:`qiime_config` settings in the scripts folder::

	print_qiime_config.py

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
----------

Quince, C., Lanzen, A., Curtis, T. P., Davenport, R. J., Hall, N., Head, I. M., et al. (2009). Accurate determination of microbial diversity from 454 pyrosequencing data. Nat Methods, 6(9), 639-641.