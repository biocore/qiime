.. _poller:

.. index:: poller.py

*poller.py* -- Poller for parallel QIIME scripts.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Script for polling parallel runs to check completion. See Qiime/scripts/`poller_examples.py <./poller_examples.html>`_ for example usage.


**Usage:** :file:`poller.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-f, `-`-check_run_complete_file
		Path to file containing a list of files that must exist to declare a run complete [REQUIRED]
	
	**[OPTIONAL]**
		
	-r, `-`-check_run_complete_f
		Function which returns True when run is completed [default: qiime.parallel.poller.basic_check_run_complete_f]
	-p, `-`-process_run_results_f
		Function to be called when runs complete [default: qiime.parallel.poller.basic_process_run_results_f]
	-m, `-`-process_run_results_file
		Path to file containing a map of tmp filepaths which should be written to final output filepaths [default: None]
	-c, `-`-clean_up_f
		Function called after processing result [default: qiime.parallel.poller.basic_clean_up_f]
	-d, `-`-clean_up_file
		List of files and directories to remove after run [default: None]
	-t, `-`-time_to_sleep
		Time to wait between calls to status_callback_f (in seconds) [default: 60]


**Output:**

No output created.



