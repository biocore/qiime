.. _parallel_qiime:

Using parallel QIIME
---------------------

QIIME supports running several of its slower steps in parallel in a cluster (or other multiple processor/core) environment. Currently, these include:

	* Assignment of taxonomy with BLAST, via :file:`Qiime/scripts/parallel_assign_taxonomy_blast.py`
	* Assignment of taxonomy with RDP, via :file:`Qiime/scripts/parallel_assign_taxonomy_rdp.py`
	* Sequence alignment with PyNAST, via :file:`Qiime/scripts/parallel_align_seqs_pynast.py`

QIIME achieves support of parallelization in different environments by requiring users to define a script which is responsible for making and starting the jobs when provided with a list of commands. This script is referred to in QIIME as the cluster_jobs script. An example cluster_jobs script which can be used for parallel runs in multicore/multiprocessor environments is packaged with QIIME as :file:`Qiime/scripts/start_parallel_jobs.py`.

Enabling parallel runs in QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To enable parallel runs in QIIME you will first need to determine if the :file:`Qiime/scripts/start_parallel_jobs.py` script will work for your purposes. If you're running in a multi-processor or multi-core environment with no queueing system, then it should work for you. If you are running in a more complex environment (e.g, a cluster) you'll need to write a custom cluster jobs script. This is discussed below.

If :file:`Qiime/scripts/start_parallel_jobs.py` will work for you, or you've defined your own cluster jobs script, you'll next need to define the number of jobs that you would like QIIME to start by default. This is done by editing the ``jobs_to_start`` value in your :file:`qiime_config` file. The default value is 1, corresponding to no parallelization. Follow the instructions on creating a custom :file:`qiime_config` (i.e., don't modify :file:`Qiime/qiime/support_files/qiime_config`, but instead copy that file to :file:`$HOME/.qiime_config` and edit that version). Then modify the jobs_to_start value to one that makes sense for your environment. For example, if you are running on a dual-core laptop, you probably want 2. (Note that this will likely prevent you from doing anything else with your laptop while QIIME is running in parallel.) If you're running on an 8 processor desktop machine, you'd want to set jobs_to_start to a maximum of 8 -- lower might be better if you'd like to reserve one or more processors for other work while running parallel QIIME. Note that setting jobs_to_start (e.g., 5 on a dual core system) to a value that is too high will reduce the performance of parallel QIIME. You can overwrite the jobs_to_start value via the command line interface of the parallel scripts -- you are just setting the default value here. If you installed QIIME using Qiime/setup.py, you will also need to set the qiime_scripts_dir value in your qiime_config file to the directory containing the QIIME scripts. By default, this will likely be /usr/local/bin. If you specified a different location by passing --install-scripts= to setup.py, you should set qiime_scripts_dir to this value.

Writing a cluster_jobs Script Specific to your Cluster Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To make QIIME parallelization useful in different computing environments users can provide a script which can start jobs on their system, referred to here as a 'cluster_jobs' script. The cluster_jobs script takes as its two parameters:

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

To avoid passing -U CLUSTER_JOBS_FP to each call to a parallel script, you should define the cluster_jobs_fp value in your :file:`qiime_config`.

Example Run of PyNAST in Parallel 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following command will start a parallel run of PyNAST, which uses the same interface as the `align_seqs.py <../scripts/align_seqs.html>`_ script, where the results are written the an output directory "parallel_aligned_seqs/"::

	parallel_align_seqs_pynast.py -i repr_set_seqs.fasta -t /ref_set_seqs.fasta -o /home/caporaso/out 

The important thing to note is that this command is that same that would be used to call serial (single processor) PyNAST, except that instead of calling `parallel_align_seqs_pynast.py <../scripts/parallel_align_seqs_pynast.html>`_, you would call `align_seqs.py <../scripts/align_seqs.html>`_ to start the run on a single processor. The output from this parallel run is the same as the output would be from the serial run. 

Details of the Parallelization 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section provides some information on details of the parallelization which are hidden from the user, but provided for users who are interested in what is happening behind-the-scenes.

The parallelization works as follows. First, the input file (-i) is split into JOBS_TO_START (-O) different roughly equal-sized files. The serial version of the script -- `align_seqs.py <../scripts/align_seqs.html>`_ -- is then called on each of these split files as a separate job. Each of these jobs therefore writes its own output files (alignment, log, and failure files). One additional job, the poller, is started to monitor each of the jobs via their output files. When all expected output files exist, the poller will merge the individual output files and clean up any temporary files including the output files created by each of the individual runs. Cleaning up temporary files can be suppressed by passing -R, which is useful for debugging. Bypassing the polling system all-together can be achieved by passing -W.
