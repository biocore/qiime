.. _qiime_config:

Setting up your qiime_config file 
==================================

First things first: you should not edit or remove :file:`Qiime/qiime/support_files/qiime_config`. Modifying or removing this file may cause QIIME to fail.

Some QIIME scripts read default values from a :file:`qiime_config` file. The default location of this file is (:file:`Qiime/qiime/support_files/qiime_config`). QIIME scripts pull default values from this file which are system-specific, such as paths to executable files, and these can be overwritten for convenience. The recommended procedure for overwriting these defaults is to copy the :file:`qiime_config` file to either :file:`$HOME/.qiime_config` or a location specified by the environment variable $QIIME_CONFIG_FP.

The QIIME configuration values should only be modified in these copies of the :file:`qiime_config` file, as changes to the :file:`Qiime/qiime/support_files/qiime_config` version may be overwritten in future QIIME updates.

When defaults are loaded, all three locations are checked in order of precedence. Lowest precedence is given to the :file:`Qiime/qiime/support_files/qiime_config` file, as these are defaults defined by the QIIME development team and are likely not relevant to other users' environments. Higher precedence is given to the file specified by $QIIME_CONFIG_FP, and this is envisioned to be used for defining system-wide defaults. Finally, highest precedence is given to :file:`~/.qiime_config`, so users have the ability to overwrite defaults defined elsewhere to have maximum control over their environment (e.g., if testing an experimental version of their :file:`cluster_jobs` script). Note that these values are defaults: the scripts typically allow overwriting of these values via their command line interfaces.

Note that users can have up to three separate :file:`qiime_config` files, and one is provided by default with QIIME. At least one :file:`qiime_config` file must be present in one of the three locations, or scripts that rely on :file:`qiime_config` file will raise an error. Not all values need to be defined in all :file:`qiime_config` files, but all values must be defined at least once. This is one more reason why you should not edit or remove :file:`Qiime/qiime/support_files/qiime_config`: when new values are added in the future they will be defined in QIIME's default copy, but not in your local copies.

Viewing and testing your ``qiime_config`` settings
----------------------------------------------------------------

To see the qiime_config values as read by QIIME, and to test your settings, you can call::

	print_qiime_config.py -t

Definition of values in qiime_config
------------------------------------

``cluster_jobs_fp`` : path to your *cluster jobs* file. This file is described in :doc:`../tutorials/parallel_qiime`.

``working_dir`` : a directory where work should be performed when running in parallel. You must be able to write to this directory. May be deprecated in favor of ``temp_dir`` in QIIME 1.8.0.

``blastmat_dir`` : directory where BLAST substitution matrices are stored

``blastall_fp`` : path to ``blastall`` executable

``pynast_template_alignment_fp`` : default template alignment to use with PyNAST as a fasta file

``pynast_template_alignment_blastdb`` : default template alignment to use with PyNAST as a pre-formatted BLAST database

``template_alignment_lanemask_fp`` : default alignment lanemask to use with ``filter_alignment.py``

``jobs_to_start`` : default number of jobs to start when running QIIME in parallel. don't make this more than the available cores/processors on your system

``seconds_to_sleep`` : number of seconds to wait when checking whether parallel jobs have completed

``temp_dir`` : directory for storing temporary files created by QIIME scripts. when a script completes successfully, any temporary files that it created are cleaned up (if you notice this isn't the case for some script, please let us know)

``denoiser_min_per_core`` : minimum number of flowgrams to denoise per core in parallel denoiser runs

``topiaryexplorer_project_dir`` : directory where TopiaryExplorer is installed

``torque_queue`` : default queue to submit jobs to when using parallel QIIME with torque

``assign_taxonomy_reference_seqs_fp`` : default reference database to use with assign_taxonomy.py (and parallel versions)

``assign_taxonomy_id_to_taxonomy_fp`` : default id-to-taxonomy map to use with assign_taxonomy.py (and parallel versions)

``sc_queue`` : default queue to submit jobs to when running parallel QIIME on StarCluster
