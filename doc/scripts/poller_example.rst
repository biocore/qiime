.. _poller_example:

.. index:: poller_example.py

*poller_example.py* -- Create python file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script is designed for use with parallel jobs to wait for their completion, and subsequently process the results and clean up. This script allows users to see it in action, and also to allow manual testing as this is a difficult process to unit test.
 
To test, call the example command below. The poller will begin running, at which time you can create the three polled files in POLLED_DIR. When all three are created, the poller will process the results, clean up, and exit.


**Usage:** :file:`poller_example.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-d, `-`-polled_dir
		Path to directory to poll
	
	**[OPTIONAL]**
		
	-P, `-`-poller_fp
		Full path to qiime/parallel/`poller.py <./poller.html>`_ [default: /Users/caporaso/code/Qiime/scripts/`poller.py <./poller.html>`_]
	-Y, `-`-python_exe_fp
		Full path to python executable [default: python]
	-c, `-`-suppress_custom_functions
		Use the default functions for checking run completion, processing results, and cleaning up (these are quiet) [default: False]


**Output:**

The poller waits for three files to be created:

 - <POLLED_DIR>/poller_test_0.txt
 - <POLLED_DIR>/poller_test_1.txt
 - <POLLED_DIR>/poller_test_2.txt
 - <POLLED_DIR> is defined via -d.

Existence of these three files is checked every 5 seconds with verbose_check_run_complete_f. When all three exist verbose_process_run_results_f 
is called, which cats all the files into a single file: <POLLED_DIR>/poller_test_completed.txt. Finally, verbose_clean_up_f is called which removes the original three files the poller was waiting on.


**Example usage:**

::

	poller_example.py -d /Users/caporaso/poller_test/

The actual call to the polling command is printed for reference just prior to calling it. This illustrates how to pass both functions and filepaths to the poller. For an example where the default (non-verbose) check_run_complete_f, process_run_results_f, and clean_up_f are used, pass -c. Again, the polling command will be printed just prior to calling:

::

	poller_example.py -d /Users/caporaso/poller_test/ -c


