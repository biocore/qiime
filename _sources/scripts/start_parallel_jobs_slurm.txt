.. _start_parallel_jobs_slurm:

.. index:: start_parallel_jobs_slurm.py

*start_parallel_jobs_slurm.py* -- Starts multiple jobs in parallel on slurm based multiprocessor systems.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script is designed to start multiple jobs in parallel on cluster systems with a slurm based scheduling system.


**Usage:** :file:`start_parallel_jobs_slurm.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-m, `-`-make_jobs
		Make the job files [default: False]
	-s, `-`-submit_jobs
		Submit the job files [default: False]
	-q, `-`-queue
		Name of queue to submit to [default: slurm's default]
	-K, `-`-mem_per_cpu
		Megabytes of memory to request per CPU [default: slurm's default]
	-j, `-`-job_dir
		Directory to store the jobs [default: jobs/]


**Output:**

No output is created.


**Job submission example:**

Start each command listed in test_jobs.txt in parallel. The run ID for these jobs will be RUNID.

::

	start_parallel_jobs_slurm.py -ms test_jobs.txt RUNID

**Queue specification example:**

Submit the commands listed in test_jobs.txt to the specified queue.

::

	start_parallel_jobs_slurm.py -ms test_jobs.txt -q himem RUNID

**Jobs output directory specification example:**

Submit the commands listed in test_jobs.txt, with the jobs put under the my_jobs/ directory.

::

	start_parallel_jobs_slurm.py -ms test_jobs.txt -j my_jobs/ RUNID


