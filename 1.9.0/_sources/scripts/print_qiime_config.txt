.. _print_qiime_config:

.. index:: print_qiime_config.py

*print_qiime_config.py* -- Print and optionally test QIIME configuration details
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Print QIIME configuration details and optionally perform tests of the QIIME base or full install.


**Usage:** :file:`print_qiime_config.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-t, `-`-test
		Test the QIIME install and configuration [default: False]
	-b, `-`-qiime_base_install
		SUPPRESSHELP
	-f, `-`-qiime_full_install
		If passed, report on dependencies required for the QIIME full install. To perform tests of the QIIME full install, you must also pass -t. [default: False]
	`-`-haiku
		SUPPRESSHELP


**Output:**

Prints QIIME configuration details to standard output.


**Example 1:**

Print basic QIIME configuration details:

::

	print_qiime_config.py

**Example 2:**

Print basic QIIME configuration details and test the base QIIME installation:

::

	print_qiime_config.py -t

**Example 3:**

Print basic QIIME configuration details and test the full QIIME installation:

::

	print_qiime_config.py -tf


