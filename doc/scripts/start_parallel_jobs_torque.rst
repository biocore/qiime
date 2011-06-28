.. _start_parallel_jobs_torque:

.. index:: start_parallel_jobs_torque.py

*start_parallel_jobs_torque.py* -- This script compares alpha 
    diversities based on a t_two_sample test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script compares the alpha 
 diversity of entries in a rarefaction file after they have been grouped 
 based on some category found in the mapping file based on a t_two_sample
 test.


**Usage:** :file:`start_parallel_jobs_torque.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-r, `-`-rarefaction_fp
		Path to rarefaction file [REQUIRED]
	-m, `-`-mapping_fp
		Path to the mapping file [REQUIRED]
	-c, `-`-category
		Category for comparison [REQUIRED]
	-d, `-`-depth
		Depth of rarefaction file to use [REQUIRED]
	-o, `-`-output_fp
		Output file path [REQUIRED]


**Output:**

Script generates an output nested dictionary which has as a first 
    key:value pair the category fed in, and a dictionary which gives the
    t_two_sample score for every possible combination of the values 
    under that category in the mapping file, saved as a text file into
    the directory specified by the output path.


**Explanation:    Inputs: mapping file lines (lines of a mapping file which associates
    to each OTU/sample a number of characteristics, given as file path),
    rarefaction file lines (lines of a rarefaction file that has scores 
    for each OTU/sample based on a certain rarefaction depth, given as a
    file path), depth (the depth score of the rarefaction file to use), 
    category (the category to compare OTU/samples on), output file path
    (a path to the output directory).:**

::

	Example: compare_alpha_diversity.py -r rarefaction_fp 
    -m mapping_fp -c 'Treatment' -d 10 -o /path/to/output.txt


