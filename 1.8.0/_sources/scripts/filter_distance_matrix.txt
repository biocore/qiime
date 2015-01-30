.. _filter_distance_matrix:

.. index:: filter_distance_matrix.py

*filter_distance_matrix.py* -- Filter a distance matrix to contain only a specified set of samples.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Remove samples from a distance matrix based on a mapping file or an otu table or a list of sample ids.


**Usage:** :file:`filter_distance_matrix.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_distance_matrix
		The input distance matrix
	-o, `-`-output_distance_matrix
		Path to store the output distance matrix
	
	**[OPTIONAL]**
		
	`-`-sample_id_fp
		A list of sample identifiers (or tab-delimited lines with a sample identifier in the first field) which should be retained
	-t, `-`-otu_table_fp
		The otu table filepath
	-m, `-`-mapping_fp
		Path to the mapping file
	-s, `-`-valid_states
		String containing valid states, e.g. 'STUDY_NAME:DOB'
	`-`-negate
		Discard specified samples (instead of keeping them) [default: False]


**Output:**




Filter samples ids listed in sample_id_list.txt from dm.txt

::

	filter_distance_matrix.py -i dm.txt -o dm_out_sample_list.txt --sample_id_fp sample_id_list.txt

Filter samples ids in otu_table.biom from dm.txt

::

	filter_distance_matrix.py -i dm.txt -o dm_out_otu_table.txt -t otu_table.biom

Filter samples ids where DOB is 20061218 in Fasting_Map.txt. (Run "`filter_samples_from_otu_table.py <./filter_samples_from_otu_table.html>`_ -h" for additional information on how metadata filtering can be specified.)

::

	filter_distance_matrix.py -i dm.txt -o dm_out_mapping_file.txt -m Fasting_Map.txt -s "DOB:20061218"


