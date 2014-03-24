.. _filter_by_metadata:

.. index:: filter_by_metadata.py

*filter_by_metadata.py* -- Filter OTU table by removal of specified metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This filter allows for the removal of sequences and OTUs that either do or don't match specified metadata, for instance, isolating samples from a specific set of studies or body sites. This script identifies samples matching the specified metadata criteria, and outputs a filtered mapping file and OTU table containing only the specified samples.


**Usage:** :file:`filter_by_metadata.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to the input OTU table (i.e., the output from `make_otu_table.py <./make_otu_table.html>`_)
	-m, `-`-map
		Path to the map file [REQUIRED]
	-s, `-`-states
		String containing valid states, e.g. 'STUDY_NAME:DOG'
	
	**[OPTIONAL]**
		
	-o, `-`-otu_outfile
		Name of otu output file, default is otu_filename.filtered.xls
	-p, `-`-map_outfile
		Name of map output file, default is map_filename.filtered.xls
	-n, `-`-num_seqs_per_otu
		Minimum counts across samples to keep OTU, default is only to keep OTUs that are present in the samples.


**Output:**

The result is a filtered OTU table and mapping file meeting the desired criteria.


**Examples:**

The following command can be used, where all options are passed (using the resulting OTU file from `make_otu_table.py <./make_otu_table.html>`_, the original Fasting_Map.txt, and keeping only the Control sequences in the Treatment field) with the resulting data being written to otu_table.txt.filtered.xls and Fasting_Map.txt.filtered.xls:

::

	filter_by_metadata.py -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control'

Some variations (not so useful on this dataset, but more useful on larger datasets) are: 

Keeping both Control and Fast in the Treatment field (i.e. keeping everything):

::

	filter_by_metadata.py -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control,Fast'

Excluding Fast in the Treatment field (same as the first example) - the syntax here is "*" to keep everything, then !Fast to eliminate the Fast group:

::

	filter_by_metadata.py -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:*,!Fast'

Keeping only samples with both Control in the Treatment field and 20061218 in the DOB field:

::

	        filter_by_metadata.py -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control;DOB:20061218'

Keeping only samples with Control in the Treatment field and OTUs with counts of at least 5 across samples:

::

	filter_by_metadata.py -i otu_table.txt -m Fasting_Map.txt -s 'Treatment:Control' -n 5

Note that the filtered mapping file will automatically exclude any columns that are the same for all the samples that are left, and will also exclude (except for SampleID) any columns that are different for all the samples that are left, making it more useful for downstream analyses with the coloring tools.


