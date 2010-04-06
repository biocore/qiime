.. _qiime_config:

Setting up your qiime_config file 
==================================

First things first: you should not edit or remove :file:`Qiime/qiime/support_files/qiime_config`. 

Some QIIME scripts read default values from a :file:`qiime_config` file. The default location of this file is (:file:`Qiime/qiime/support_files/qiime_config`). QIIME scripts pull default values from this file which are system-specific, such as paths to executable files, and these can be overwritten for convenience. The recommended procedure for overwriting these defaults is to copy the :file:`qiime_config` file to either :file:`~/.qiime_config` or a location specified by the environment variable $QIIME_CONFIG_FP.

The Qiime configuration values should only be modified in these copies of the :file:`qiime_config` file, as changes to the :file:`Qiime/qiime/support_files/qiime_config` version may be overwritten in future QIIME updates.

When defaults are loaded, all three locations are checked in order of precedence. Lowest precedence is given to the :file:`Qiime/qiime/support_files/qiime_config` file, as these are defaults defined by the QIIME development team and are likely not relevant to other users' environments. Higher precedence is given to the file specified by $QIIME_CONFIG_FP, and this is envisioned to be used for defining system-wide defaults. Finally, highest precedence is given to :file:`~/.qiime_config`, so users have the ability to overwrite defaults defined elsewhere to have maximum control over their environment (e.g., if testing an experimental version of their :file:`cluster_jobs` script). Note that these values are defaults: the scripts typically allow overwriting of these values via their command line interfaces.

Note that users can have up to three separate :file:`qiime_config` files, and one is provided by default with QIIME. At least one :file:`qiime_config` file must be present in one of the three locations, or scripts that rely on :file:`qiime_config` file will raise an error. Not all values need to be defined in all :file:`qiime_config` files, but all values must be defined at least once. This is one more reason why you should not edit or remove :file:`Qiime/qiime_config`: when new values are added in the future they will be defined in Qiime's default copy, but not in your local copies.

To see the qiime_config values as read by QIIME, and to test your settings, you can call::

	Qiime/scripts/print_qiime_config.py -t

Setting qiime_scripts_dir
--------------------------
If you have installed QIIME using its ``setup.py`` script, you will need to set the ``qiime_scripts_dir`` value in your ``qiime_config`` file to the directory containing the QIIME scripts. By default, this will likely be ``/usr/local/bin/``. If you specified a different location by passing ``--install-scripts=`` to ``setup.py``, you should set ``qiime_scripts_dir`` to this value.
