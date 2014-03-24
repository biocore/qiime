.. _start_parallel_jobs_sc:

.. index:: start_parallel_jobs_sc.py

*start_parallel_jobs_sc.py* -- Starts parallel jobs on Sun GridEngine queueing systems.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Starts multiple jobs in parallel on Sun GridEngine systems. This is designed to work with StarCluster EC2 instances, but may be applicable beyond there.


**Usage:** :file:`start_parallel_jobs_sc.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-m, `-`-make_jobs
		Make the job files [default: None]
	-s, `-`-submit_jobs
		Submit the job files [default: None]
	-q, `-`-queue_name
		The queue to submit jobs to [default: all.q]


**Output:**

No output is created.


**Job submission example:**

Start each command listed in test_jobs.txt in parallel. The run ID for these jobs will be RUNID.

::

	start_parallel_jobs_sc.py -ms test_jobs.txt RUNID

**Queue specification example:**

Submit the commands listed in test_jobs.txt to the specified queue.

::

	start_parallel_jobs_sc.py -ms test_jobs.txt -q all.q RUNID


