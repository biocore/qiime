.. _filter_otu_table:

.. index:: filter_otu_table.py

*filter_otu_table.py* -- Filters OTU table by minimum OTU count and number of samples or by taxonomy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

After the OTU has been generated, the user may want to filter the table based on the number of samples within each OTU or by the number of sequences per OTU. This step is generally done to reduce the noise within the OTU table and can also reduce the overall size of the table, which is essential when performing analyses on large datasets. Along with filtering based on samples or sequences, the user can include and exclude specific taxon groups.


**Usage:** :file:`filter_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	
	**[OPTIONAL]**
		
	-c, `-`-min_count
		Retain OTUs with at least this many sequences [default=1]
	-s, `-`-min_samples
		Retain OTUs found in at least this many samples [default=2]
	-t, `-`-include_taxonomy
		List of taxonomy terms to include [default=]
	-e, `-`-exclude_taxonomy
		List of taxonomy terms to exclude [default=]
	-o, `-`-dir_path
		Directory prefix for all analyses [default=./]
	-p, `-`-seqs_per_sample
		Minimum sequences per sample to retain the sample. [default=None]


**Output:**

The result of `filter_otu_table.py <./filter_otu_table.html>`_ creates a new OTU table, where the filename uses the input OTU filename and appends "filtered.txt" to the end of the name.


**Examples:**

To filter the OTU table using the default parameters ("-c 1" and "-s 2"), then write the results to the current working directory, you can use the code as follows:

::

	filter_otu_table.py -i otu_table.txt

To filter by the number of samples (i.e., 5) within each OTU (keep only OTU's found in at least X number of samples), you can use the code as follows:

::

	filter_otu_table.py -i otu_table.txt -s 5

To filter by the number of sequences (i.e., 5) within each OTU (keep only OTU's with at least X sequences in the OTU), you can use the code as follows:

::

	filter_otu_table.py -i otu_table.txt -c 5

To include ("Bacteria") and exclude ("Proteobacteria") certain taxon groups (options -t / -e respectively), you can use the code as follows.  The include and exclude parameters must be used together:

::

	filter_otu_table.py -i otu_table.txt -t Bacteria -e Proteobacteria

**Filter samples by number of sequences:**

A user may want to remove samples that have low sequence coverage. NOTE: this feature is mutually exclusive from the other filtering options, so if you pass this, you will need to perform a subsequent filter to remove by the other options.

::

	filter_otu_table.py -i otu_table.txt -p 10


