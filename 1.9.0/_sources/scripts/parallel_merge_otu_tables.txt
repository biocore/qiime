.. _parallel_merge_otu_tables:

.. index:: parallel_merge_otu_tables.py

*parallel_merge_otu_tables.py* -- Parallel merge BIOM tables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script works like the `merge_otu_tables.py <./merge_otu_tables.html>`_ script, but is intended to make use of multicore/multiprocessor environments to perform analyses in parallel.


**Usage:** :file:`parallel_merge_otu_tables.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fps
		The otu tables in biom format (comma-separated)
	-o, `-`-output_dir
		The output otu table directory path
	
	**[OPTIONAL]**
		
	-C, `-`-cluster
		Submit to a torque cluster
	-Z, `-`-seconds_to_sleep
		Number of seconds to sleep between checks for run  completion when polling runs [default: 1]
	-X, `-`-job_prefix
		Job prefix [default: descriptive prefix + random chars]


**Output:**

The output consists of many files (i.e. merged_table.biom, merged_table.log and all intermediate merge tables). The .biom file contains the result of merging the individual BIOM tables. The resulting .log file contains a list of parameters passed to this script along with the output location of the resulting .txt file, the dependency hierarchy and runtime information for each individual merge.


**Example:**

Merge the OTU tables $PWD/t1.biom,$PWD/t2.biom,$PWD/t3.biom,$PWD/t4.biom and write the resulting output table to the $PWD/merged/ directory.

::

	parallel_merge_otu_tables.py -i $PWD/t1.biom,$PWD/t2.biom,$PWD/t3.biom,$PWD/t4.biom -o $PWD/merged/


