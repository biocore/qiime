.. _collate_alpha:

.. index:: collate_alpha.py

*collate_alpha.py* -- Collate alpha diversity results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

When performing batch analyses on the OTU table (e.g. rarefaction followed by alpha diversity), the result of `alpha_diversity.py <./alpha_diversity.html>`_ comprises many files, which need to be concatenated into a single file for generating rarefaction curves.  This script joins those files.
Input files are:
each file represents one (rarefied) otu table
each row in a file represents one sample
each column in a file represents one diversity metric

Output files are:
each file represents one diversity metric
each row in a file represents one (rarefied) otu table 
each column in a file represents one sample

The input directory should contain only otu tables. The output directory should be empty or nonexistant and the example file is optional.  

If you have a set of rarefied OTU tables, make sure the example file contains every sample present in the otu talbes. You should typically choose the file with the fewest sequences per sample, to avoid files with sparse samples omitted.



**Usage:** :file:`collate_alpha.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path (a directory)
	-o, `-`-output_path
		Output path (a directory).  will be created if needed
	
	**[OPTIONAL]**
		
	-e, `-`-example_path
		Example alpha_diversity analysis file, containing all samples and all metrics to be included in the collated result[Default: chosen automatically (see usage string)]


**Output:**

This script takes the resulting files from batch alpha diversity and collates them into (one file for each metric used).

This script transforms a series of files, named (e.g. alpha_rarefaction_20_0.txt, alpha_rarefaction_20_1.txt, etc.) into a (usually much smaller) set of files named (e.g. chao1.txt, PD_whole_tree.txt, etc.), where the columns correspond to samples and rows to the rarefaction files inputted, as shown by the following:

==========================  ====================    =========   ======  ======
\                           sequences per sample    iteration   PC.354  PC.355
==========================  ====================    =========   ======  ======
alpha_rarefaction_20_0.txt  20                      0           0.925   0.915 
alpha_rarefaction_20_1.txt  20                      1           0.9     0.89 
alpha_rarefaction_20_2.txt  20                      2           0.88    0.915 
alpha_rarefaction_20_3.txt  20                      3           0.91    0.93 
...                         ...                     ...         ...     ...
==========================  ====================    =========   ======  ====== 




**Example:**

The user inputs the results from batch alpha diversity (e.g. alpha_div/) and the location where the results should be written (e.g. collated_alpha/), as shown by the following command:

::

	collate_alpha.py -i alpha_div/ -o collated_alpha/


