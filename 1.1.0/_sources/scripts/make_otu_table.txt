.. _make_otu_table:

.. index:: make_otu_table.py

*make_otu_table.py* -- Make OTU table
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The script `make_otu_table.py <./make_otu_table.html>`_ tabulates the number of times an OTU is found in each sample, and adds the taxonomic predictions for each OTU in the last column if a taxonomy file is supplied.


**Usage:** :file:`make_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_map_fp
		Path to the input OTU map (i.e., the output from `pick_otus.py <./pick_otus.html>`_)
	
	**[OPTIONAL]**
		
	-t, `-`-taxonomy
		Path to taxonomy assignment, containing the assignments of \ taxons to sequences (i.e., resulting txt file from `assign_taxonomy.py <./assign_taxonomy.html>`_)  [default: None]
	-o, `-`-output_fp
		The output filepath


**Output:**

The output of `make_otu_table.py <./make_otu_table.html>`_ is a tab-delimited text file, where the columns correspond to Samples and rows correspond to OTUs and the number of times a sample appears in a particular OTU.


**Example:**

For this example the input is an OTU file containing sequence ids assigned to each OTU (i.e., resulting OTU file from `pick_otus.py <./pick_otus.html>`_) and a text file containing the taxonomy assignments (i.e., resulting text file from `assign_taxonomy.py <./assign_taxonomy.html>`_), where the output file is defined as otu_table.txt:

::

	make_otu_table.py -i seqs_otus.txt -t repr_set_tax_assignments.txt -o otu_table.txt


