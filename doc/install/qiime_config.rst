.. _qiime_config:

Setting up your qiime config file
==================================

Your QIIME installation can be customized using a *qiime config* file.

You have two options of where your *qiime config* file can be stored, and you can use either or both of these options. First, you can place a file called ``.qiime_config`` in your home directory (i.e., ``$HOME/.qiime_config``). Second, you can put the *qiime config* file in a location of your choice, and define an environment variable called ``QIIME_CONFIG_FP`` whose value is the absolute path to the *qiime config* file. If you use both of these options, the settings that are defined in ``$HOME/.qiime_config`` will take precedence over the settings that are defined in ``$QIIME_CONFIG_FP``.

The values that you can set in your *qiime config* file are listed below. Defining any of these is optional. QIIME will use built-in defaults for any of these values that you don't set.

- ``assign_taxonomy_id_to_taxonomy_fp`` : id-to-taxonomy map to use with assign_taxonomy.py (and parallel versions), if you prefer to not use the default
- ``assign_taxonomy_reference_seqs_fp`` : reference database to use with assign_taxonomy.py (and parallel versions), if you prefer to not use the default
- ``blastmat_dir`` : directory where BLAST substitution matrices are stored
- ``blastall_fp`` : path to ``blastall`` executable
- ``cluster_jobs_fp`` : path to your *cluster jobs* file. This file is described in :doc:`../tutorials/parallel_qiime`.
- ``denoiser_min_per_core`` : minimum number of flowgrams to denoise per core in parallel denoiser runs
- ``jobs_to_start`` : default number of jobs to start when running QIIME in parallel.
- ``pick_otus_reference_seqs_fp`` : reference database to use with all OTU picking scripts and workflows, if you prefer to not use the default
- ``pynast_template_alignment_blastdb`` : template alignment to use with PyNAST as a pre-formatted BLAST database, if you prefer to not have a BLAST database constructed from the fasta filepath provided for ``pynast_template_alignment_fp``
- ``pynast_template_alignment_fp`` : template alignment to use with PyNAST as a fasta file, if you prefer to not use the default
- ``sc_queue`` : default queue to submit jobs to when running parallel QIIME on `StarCluster-based <http://star.mit.edu/cluster/>`_ Amazon Web Services clusters
- ``seconds_to_sleep`` : number of seconds to wait when checking whether parallel jobs have completed
- ``slurm_memory`` : amount of memory in megabytes to request per CPU, when using ``start_parallel_jobs_slurm.py``; will default to slurm's default
- ``slurm_queue`` : queue to submit jobs to, when using ``start_parallel_jobs_slurm.py``; will default to slurm's default
- ``temp_dir`` : directory for storing temporary files created by QIIME scripts. when a script completes successfully, any temporary files that it created are cleaned up. If you're working in a cluster environment, this directory must be shared across all of the worker nodes that QIIME jobs may be running on.
- ``topiaryexplorer_project_dir`` : directory where TopiaryExplorer is installed
- ``torque_queue`` : default queue to submit jobs to when using ``start_parallel_jobs_torque.py``

Viewing and testing your ``qiime_config`` settings
--------------------------------------------------

To see the qiime_config values as read by QIIME, run::

	print_qiime_config.py -t

We refer to the output of this command as the *print_qiime_config.py output*.

If you have added some of the above settings to a *qiime config* file, and you do not see those settings associated with the values in your print_qiime_config.py output, that means that QIIME did not find your *qiime config* file. Check that it is in one of the two allowed locations described above.

Default reference files
-----------------------

Several of the above settings will fall back to default files if their value is ``None`` in the print_qiime_config.py output. These are:

- ``assign_taxonomy_id_to_taxonomy_fp``
- ``assign_taxonomy_reference_seqs_fp``
- ``pick_otus_reference_seqs_fp``
- ``pynast_template_alignment_fp``

In the *QIIME default reference information* section of the print_qiime_config.py output, you'll be presented with a link to a page that describes what files are used as the defaults for each of these. These files come from the `qiime-default-reference project <http://github.com/biocore/qiime-default-reference>`_, one of QIIME's core dependencies. (We don't describe those here as the defaults may differ based on how recent your QIIME installation is.)
