.. _parallel_qiime:

Using parallel QIIME
====================

Some steps in QIIME can be run in parallel, either directly or though workflow scripts that wrap them. Some of the individual steps that can be run in parallel are assignment of taxonomy with BLAST (:doc:`../scripts/parallel_assign_taxonomy_blast`) or the RDP classifier (:doc:`../scripts/parallel_assign_taxonomy_rdp`), sequence alignment with PyNAST (:doc:`../scripts/parallel_align_seqs_pynast`), computation of alpha diversity (:doc:`../scripts/parallel_alpha_diversity`) and beta diversity (:doc:`../scripts/parallel_beta_diversity`), and reference-based OTU picking (:doc:`../scripts/parallel_pick_otus_uclust_ref`). The workflow scripts that wrap these steps allow you to run them in parallel by passing the ``-a`` parameter.

How parallel QIIME works
-------------------------

The parallel steps in QIIME are designed to reduce the wall time of computationally expensive steps and make use of as many as possible of the types of parallel systems that biologists are frequently working with. This includes multi-core/multi-processor systems without queueing systems (e.g., laptop or desktop machines, or single AWS cloud instances) or with queueing systems (e.g., local compute clusters or AWS clusters built with `StarCluster <http://star.mit.edu/cluster/>`_). 

QIIME achieves support of parallelization in these different environments by working on the assumption that a list of commands can be generated which are then submitted (somehow) to run on the system. The `somehow` could, for example, be simply by issuing one system call per command on a system with no job queueing or by creating a PBS job file per command and submitting those via  ``qsub`` on a cluster running torque. A parallel QIIME script thus builds a list of commands that should be run, and passes that list to the user's `cluster jobs` script. The `cluster jobs` script is either one of several that ship with QIIME, or one that the user writes to customize how jobs are submitted in their environment. 

The `cluster jobs` scripts that ship with QIIME are:

 * :doc:`../scripts/start_parallel_jobs`
 * :doc:`../scripts/start_parallel_jobs_torque`
 * :doc:`../scripts/start_parallel_jobs_sc`

Enabling parallel runs in QIIME
-------------------------------

By default, QIIME is configured to run parallel jobs on systems without a queueing system (e.g., your laptop, desktop, or single AWS instance), where you tell a parallel QIIME script how many jobs should be submitted. If you're working on one of these types of systems, you'll call one of  QIIME's parallel scripts or pass ``-a`` to one of QIIME's workflow scripts to run in parallel. You'll also pass ``-O`` with the number of jobs to run. For example, to run ``pick_reference_otus_through_otu_table.py`` via four parallel jobs, you can do the following::

	pick_reference_otus_through_otu_table.py -a -O 4 ...

where ``...`` is replaced with the remaining parameters to the script. 

You can specify a default number of jobs via the ``jobs_to_start`` parameter in ``qiime_config`` (see :doc:`../install/qiime_config` for help with setting up your ``qiime_config`` file) to avoid having to pass ``-O``. The default value for ``jobs_to_start`` is ``1``, corresponding to no parallelization (hence your having to pass ``-O`` to use parallel QIIME unless you've overridden the default). You should modify the ``jobs_to_start`` value to one that makes sense for your environment. For example, if you are running on a dual-core laptop, you probably want to specify ``2``. (Note that this will likely prevent you from doing anything else with your laptop while QIIME is running in parallel.) If you're running on an 8 processor machine, you'd want to set ``jobs_to_start`` to a maximum of ``8``, but ``7`` might be better if you'd like to reserve one processor for other work while running parallel QIIME. Note that setting ``jobs_to_start`` to a value that is too high (e.g., ``5`` on a quad-core system) will cause your job to take longer to complete than if you specify a value that makes sense for your environment. Finally, you can always override your default ``jobs_to_start`` value by passing ``-O`` to a parallel script: you are just setting that default value in ``qiime_config``.

If you are running in a more complex environment (e.g, a cluster) you'll need to determine if one of the QIIME `cluster jobs` scripts will work for you, or whether you'll need to write a custom `cluster jobs` script (discussed below). In either of these cases, you'll overwrite the ``cluster_jobs_fp`` value in your ``qiime_config`` file to be the full path to the `cluster jobs` script that QIIME should use, or just the name of the script if it is in a directory in your ``$PATH`` environment variable.

.. warning:: 
	
	Before starting parallel jobs with QIIME, you should run ``print_qiime_config.py -t`` to confirm that the changes you've made in ``qiime_config`` have been recognized by QIIME. This is very important as it allows you to ensure that the correct ``cluster_jobs_fp`` is being used in your environment (and therefore that you're not about to issue 100 ``system`` calls on the head node of your cluster, which would likely make your system administrator very angry - you've been warned!). 

Writing a cluster jobs script specific to your parallel environment
-------------------------------------------------------------------

To make QIIME parallelization useful in different computing environments users can provide a script which can start jobs on their system, referred to here as a `cluster jobs` script. The `cluster jobs` script takes exactly two parameters:

	1. A single file which lists the commands to be run (referred to as a `jobs list` file), with one command per line.
	2. A string to use as a prefix when constructing unique job identifiers.

The lines in an example `jobs list` file might be::

	pick_otus.py -i inseqs_file1.fasta 
	pick_otus.py -i inseqs_file2.fasta 
	pick_otus.py -i inseqs_file3.fasta 

If passed to your `cluster jobs` script, this should start three separate jobs corresponding to each of the commands.

The call to the `cluster jobs` script from QIIME's parallel scripts looks like the following (so your script must adhere to this interface)::

	CLUSTER_JOBS_FP -ms job_list.txt JOB_ID

where ``CLUSTER_JOBS_FP`` is the path to your `cluster jobs` script and is passed to the parallel scripts via the ``-U`` parameter (or you can define it with the ``cluster_jobs_fp`` variable in your ``qiime_config``). ``JOB_ID`` is intended to be used as a prefix by the `cluster jobs` script when creating a unique identifier for each job. The same ``JOB_ID`` is also used by the QIIME parallel scripts when creating names for temporary files and directories, but your script does not necessarily need to do anything with this information if it's not useful to you. The ``-ms`` indicates that the `job files` should be made (``-m``) and submitted (``-s``).

Once you have written a `cluster jobs` script for your specific environment that can be called via the above interface, running QIIME jobs in parallel should be straight-forward. The parallel variants of the scripts use the same parameters as the serial versions of the scripts, with some additional options in the parallel scripts.

