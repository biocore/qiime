.. _start_parallel_jobs_torque:

.. index:: start_parallel_jobs_torque.py

*start_parallel_jobs_torque.py* -- Starts multiple jobs in parallel on torque/qsub based multiprocessor systems.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script is designed to start multiple jobs in parallel on cluster systems with a torque/qsub based scheduling system.


**Usage:** :file:`start_parallel_jobs_torque.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-m, `-`-make_jobs
		Make the job files [default: None]
	-s, `-`-submit_jobs
		Submit the job files [default: None]
	-q, `-`-queue
		Name of queue to submit to  [default: friendlyq]
	-j, `-`-job_dir
		Directory to store the jobs [default: jobs/]


**Output:**

No output is created.


**Example:**

Start each command listed in test_jobs.txt in parallel. The run id for these jobs will be RUNID. 

::

	start_parallel_jobs_torque.py -ms test_jobs.txt RUNID


