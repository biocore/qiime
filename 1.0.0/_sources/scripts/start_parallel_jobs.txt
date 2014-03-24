.. _start_parallel_jobs:

.. index:: start_parallel_jobs.py

*start_parallel_jobs.py* -- Starts multiple jobs in parallel on multicore or multiprocessor systems.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script is designed to start multiple jobs in parallel on systems with no queueing system, for example a multiple processor or multiple core laptop/desktop machine. This also serves as an example 'cluster_jobs' which users can use a template to define scripts to start parallel jobs in their environment.


**Usage:** :file:`start_parallel_jobs.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-m, `-`-make_jobs
		Make the job files [default: None]
	-s, `-`-submit_jobs
		Submit the job files [default: None]


**Output:**

No output is created.


**Example:**

Start each command listed in test_jobs.txt in parallel. The run id for these jobs will be RUNID. 

::

	start_parallel_jobs.py -ms test_jobs.txt RUNID


