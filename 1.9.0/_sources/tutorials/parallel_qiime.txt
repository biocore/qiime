.. _parallel_qiime:

Using parallel QIIME
====================

Some scripts in QIIME can be run in parallel, either directly or through workflow scripts that wrap them. Some of the individual scripts that can be run in parallel are:

 - :doc:`../scripts/parallel_assign_taxonomy_blast`
 - :doc:`../scripts/parallel_assign_taxonomy_rdp`
 - :doc:`../scripts/parallel_align_seqs_pynast`
 - :doc:`../scripts/parallel_alpha_diversity`
 - :doc:`../scripts/parallel_beta_diversity`
 - :doc:`../scripts/parallel_pick_otus_uclust_ref`

To find a complete list of the parallel QIIME scripts, see the scripts beginning with the name ``parallel`` on :doc:`the script index <../scripts/index>`. The workflow scripts that wrap these scripts allow you to run them in parallel by passing the ``-a`` parameter.

How parallel QIIME works
------------------------

The parallel scripts in QIIME are designed to reduce the wall time of computationally expensive steps by making use of many types of parallel systems that biologists frequently work with. Parallel QIIME supports multi-core/multi-processor systems without queueing systems (e.g., laptop or desktop machines, or single AWS cloud instances) or systems with queueing systems (e.g., local compute clusters or AWS clusters built with `StarCluster <http://star.mit.edu/cluster/>`_).

QIIME supports parallelization in these different environments by assuming that a list of commands can be generated which are then submitted (somehow) to run on the system. The `somehow` could, for example, be simply issuing a single system call per command on a system with no job queueing, or by creating a PBS job file per command and submitting those via ``qsub`` on a cluster running torque. A parallel QIIME script thus builds a list of commands that should be run, and passes that list to the user's `cluster jobs` script. The `cluster jobs` script is either one of several that ship with QIIME, or one that the user writes to customize how jobs are submitted in their environment.

The `cluster jobs` scripts that ship with QIIME are:

 * :doc:`../scripts/start_parallel_jobs`
 * :doc:`../scripts/start_parallel_jobs_torque`
 * :doc:`../scripts/start_parallel_jobs_sc`
 * :doc:`../scripts/start_parallel_jobs_slurm`

Enabling parallel QIIME
-----------------------

By default, QIIME is configured to run parallel jobs on systems without a queueing system (e.g., your laptop, desktop, or single AWS instance), where you tell a parallel QIIME script how many jobs should be submitted. If you're working on one of these types of systems, you'll either run one of  QIIME's parallel scripts directly or pass ``-a`` to one of QIIME's workflow scripts to run in parallel. You'll also pass ``-O`` with the number of parallel jobs to run. For example, to run the `pick_closed_reference_otus.py <../scripts/pick_closed_reference_otus.html>`_ workflow with four parallel jobs::

	pick_closed_reference_otus.py -a -O 4 ...

where ``...`` is replaced with the remaining parameters to the script.

The following subsections describe how to set up parallel QIIME for your specific environment through your QIIME config file. See :doc:`../install/qiime_config` for help with setting up your QIIME config file.

Default number of parallel jobs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can specify a default number of parallel jobs via the ``jobs_to_start`` parameter in your QIIME config file to avoid having to pass ``-O`` each time you run a parallel script or workflow. The default value for ``jobs_to_start`` is ``1``, corresponding to no parallelization (hence your having to pass ``-O`` to use parallel QIIME unless you've overridden the default). You should modify the ``jobs_to_start`` value to one that makes sense for your environment. For example, if you are running on a dual-core laptop, you probably want to specify ``2``. (Note that this will likely prevent you from doing anything else with your laptop while QIIME is running in parallel.) If you're running on an 8 processor machine, you'd want to set ``jobs_to_start`` to a maximum of ``8``, but ``7`` might be better if you'd like to reserve one processor for other work while running parallel QIIME. Note that setting ``jobs_to_start`` to a value that is too high (e.g., ``5`` on a quad-core system) will cause your job to take longer to complete than if you specify a value that makes sense for your environment. Finally, you can always override your default ``jobs_to_start`` value by passing ``-O`` to a parallel script: you are just setting the *default* value in your QIIME config file.

Default cluster jobs script
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you are running in a more complex environment (e.g, a cluster), you'll need to determine if one of the QIIME `cluster jobs` scripts will work for you or whether you'll need to write a custom `cluster jobs` script (discussed below). In either of these cases, you'll set the ``cluster_jobs_fp`` value in your QIIME config file to be the absolute path to the `cluster jobs` script that QIIME should use, or just the name of the script if it is in a directory in your ``$PATH`` environment variable.

.. warning::

	Before starting parallel jobs with QIIME, you should run ``print_qiime_config.py -t`` to confirm that the changes you've made in your QIIME config have been recognized by QIIME. This is very important as it allows you to ensure that the correct ``cluster_jobs_fp`` is being used in your environment (and therefore that you're not about to issue 100 system calls on the head node of your cluster, which would likely make your system administrator very angry - you've been warned!).

.. warning::

	If you're using the QIIME workflow scripts in parallel mode (i.e., with the ``-a`` parameter), and submitting the workflow command as a job to the queueing system, that job must be able to submit other jobs to the queue. In other words, worker jobs on the cluster must have sufficient permission to submit jobs.

Specifying a temporary directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may specify a temporary directory in your QIIME config file by setting the ``temp_dir`` value. QIIME scripts use this temporary directory to store temporary files, which are cleaned up after the script completes. When setting up parallel QIIME to work in your environment, the temporary directory must meet the following requirements:

 - If you’re working in a cluster environment, the temporary directory must be shared across all of the worker nodes that QIIME jobs may be running on.
 - If you're working in a multi-user environment (e.g., a cluster or a server with more than one person running QIIME), the temporary directory must be unique to each user. There are several ways to accomplish this. For example, each user can create their own QIIME config file that points to their temporary directory. On a system with many users, it may be desirable to instead create a system-wide QIIME config file (specified via ``$QIIME_CONFIG_FP``) that uses ``$USER`` in the path to the temporary directory. For example, you could add the following line to the system-wide QIIME config file::

    temp_dir /tmp/$USER

Additional customization of parallel QIIME
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are additional values that may be specified in the QIIME config file to further customize your parallel QIIME environment. These values are described in :doc:`../install/qiime_config`.

Writing a cluster jobs script specific to your parallel environment
-------------------------------------------------------------------

To make parallel QIIME usable in different computing environments, users may provide a custom script which can start jobs on their system, referred to here as a `cluster jobs` script. The `cluster jobs` script takes exactly two parameters:

	1. A single file which lists the commands to be run (referred to as a `jobs list` file), with one command per line.
	2. A string to use as a prefix when constructing unique job identifiers.

The lines in an example `jobs list` file might be::

	pick_otus.py -i inseqs_file1.fasta
	pick_otus.py -i inseqs_file2.fasta
	pick_otus.py -i inseqs_file3.fasta

If passed to your `cluster jobs` script, this should start three separate jobs corresponding to each of the commands.

The call to the `cluster jobs` script from within QIIME's parallel scripts looks like the following (so your script must adhere to this interface)::

	CLUSTER_JOBS_FP -ms job_list.txt JOB_ID

where ``CLUSTER_JOBS_FP`` is the path to your `cluster jobs` script and is passed to the parallel scripts via the ``-U`` parameter (or you can define it with the ``cluster_jobs_fp`` variable in your QIIME config file, as described above). ``JOB_ID`` is intended to be used as a prefix by the `cluster jobs` script when creating a unique identifier for each job. The same ``JOB_ID`` is also used by the QIIME parallel scripts when creating names for temporary files and directories, but your script does not necessarily need to do anything with this information if it's not useful to you. The ``-ms`` indicates that the `job files` should be made (``-m``) and submitted (``-s``).

Once you have written a `cluster jobs` script for your specific environment that can be called via the above interface, running QIIME jobs in parallel should be straightforward. The parallel variants of the scripts use the same parameters as the serial versions of the scripts, with some additional options related to parallel execution.

The poller (and figuring out what the poller is waiting for)
------------------------------------------------------------

Most of the parallel QIIME scripts end by collating results to look like those generated by the non-parallel variant of the script. The `QIIME poller` is used to determine when all of the individual jobs are complete and at that time initiate the collation process. The poller determines when all jobs are completed by reading a ``check_run_complete_file`` that is generated by the parallel script, containing the paths to all of the expected output files. When all of the paths listed in the ``check_run_complete_file`` exist, the poller concludes that all jobs have finished and the collation process can begin.

Sometimes one of the parallel jobs will fail and its output files will not be written to the expected location. This will cause the poller to wait indefinitely. You can use the `identify_missing_files.py <../scripts/identify_missing_files.html>`_ script to identify which files the poller is still waiting on by running it with the path to the ``check_run_complete_file``. The ``check_run_complete_file`` will be called ``expected_out_files.txt`` and found in a temporary directory under the output directory for the parallel script. If you determine that the poller is still waiting on some files and you think that the job(s) that would generate those files are no longer running, you can identify the command that failed by looking for the missing output file name(s) in the `jobs list` file (also under the output directory for the parallel script, and having a filename ending with ``_jobs.txt``), and re-running those specific commands. In the future we hope to improve parallel job failure recovery in QIIME as we realize that this is fairly tedious.
