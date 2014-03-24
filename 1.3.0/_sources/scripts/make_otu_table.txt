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
	-e, `-`-exclude_otus_fp
		A filepath listing OTU identifiers that should not be included in the OTU table (e.g., the output of `identify_chimeric_seqs.py <./identify_chimeric_seqs.html>`_)


**Output:**

The output of `make_otu_table.py <./make_otu_table.html>`_ is a tab-delimited text file, where the columns correspond to Samples and rows correspond to OTUs and the number of times a sample appears in a particular OTU.


Make an OTU table from an OTU map (i.e., result from `pick_otus.py <./pick_otus.html>`_) and a taxonomy assignment file (i.e., result from `assign_taxonomy.py <./assign_taxonomy.html>`_). Write the output file to otu_table.txt.

::

	make_otu_table.py -i otu_map.txt -t tax_assignments.txt -o otu_table.txt

Make an OTU table, excluding the sequences listed in chimeric_seqs.txt

::

	make_otu_table.py -i otu_map.txt -o otu_table.txt -e chimeric_seqs.txt


