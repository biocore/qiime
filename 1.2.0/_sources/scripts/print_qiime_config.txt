.. _print_qiime_config:

.. index:: print_qiime_config.py

*print_qiime_config.py* -- Print out the qiime config settings.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

A simple scripts that prints out the qiime config settings and does some sanity checks.


**Usage:** :file:`print_qiime_config.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-t, `-`-test
		Test the qiime config for sanity [default: False]


**Output:**

This prints the qiime_config to stdout.


**Example 1:**

Print qiime config settings:

::

	print_qiime_config.py

**Example 2:**

Print and check qiime config settings for sanity:

::

	print_qiime_config.py -t


