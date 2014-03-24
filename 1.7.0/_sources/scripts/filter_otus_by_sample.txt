.. _filter_otus_by_sample:

.. index:: filter_otus_by_sample.py

*filter_otus_by_sample.py* -- Filter OTU mapping file and sequences by SampleIDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This filter allows for the removal of sequences and OTUs containing user-specified Sample IDs, for instance, the removal of negative control samples. This script identifies OTUs containing the specified Sample IDs and removes its corresponding sequence from the sequence collection.


**Usage:** :file:`filter_otus_by_sample.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_map_fp
		Path to the input OTU map (i.e., the output from `pick_otus.py <./pick_otus.html>`_)
	-f, `-`-input_fasta_fp
		Path to the input fasta file
	-s, `-`-samples_to_extract
		This is a list of sample ids, which should be removed from the OTU file
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Path to the output directory


**Output:**

As a result a new OTU and sequence file is generated and written to a randomly generated folder where the name of the folder starts with "filter_by_otus" Also included in the folder, is another FASTA file containing the removed sequences, leaving the user with 3 files.


**Example:**

The following command can be used, where all options are passed (using the resulting OTU file from `pick_otus.py <./pick_otus.html>`_, FASTA file from `split_libraries.py <./split_libraries.html>`_ and removal of sample 'PC.636') with the resulting data being written to the output directory "filtered_otus/":

::

	filter_otus_by_sample.py -i seqs_otus.txt -f seqs.fna -s PC.636 -o filtered_otus/


