.. _identify_missing_files:

.. index:: identify_missing_files.py

*identify_missing_files.py* -- This script checks for the existence expected file in parallel runs.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script checks for the existence expected file in parallel runs, and is useful for checking the status of a parallel run or for finding out what `poller.py <./poller.html>`_ is waiting on in a possibly failed run.


**Usage:** :file:`identify_missing_files.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-e, `-`-expected_out_fp
		The list of expected output files


**Output:**




Check for the existence of files listed in the expected_out_files.txt from a PyNAST alignment run, and print a warning for any that are missing.

::

	identify_missing_files.py -e ALIGN_BQ7_/expected_out_files.txt


