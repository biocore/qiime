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
		Name of queue to submit to [default: friendlyq]
	-j, `-`-job_dir
		Directory to store the jobs [default: jobs/]
	-w, `-`-max_walltime
		Maximum time in hours the job will run for [default: 72]
	-c, `-`-cpus
		Number of CPUs to use [default:1]
	-n, `-`-nodes
		Number of nodes to use [default:1]


**Output:**

No output is created.


**Job submission example:**

Start each command listed in test_jobs.txt in parallel. The run ID for these jobs will be RUNID.

::

	start_parallel_jobs_torque.py -ms test_jobs.txt RUNID

**Queue specification example:**

Submit the commands listed in test_jobs.txt to the specified queue.

::

	start_parallel_jobs_torque.py -ms test_jobs.txt -q friendlyq RUNID

**Jobs output directory specification example:**

Submit the commands listed in test_jobs.txt, with the jobs put under the my_jobs/ directory.

::

	start_parallel_jobs_torque.py -ms test_jobs.txt -j my_jobs/ RUNID


