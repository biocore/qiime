.. _pool_by_metadata:

.. index:: pool_by_metadata.py

*pool_by_metadata.py* -- pool samples in OTU table and mapping file based on sample metadata from mapping file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

this script outputs a new otu table and mapping file with some samples removed and replaced with one pooled sample. the new pooled sample will have fields in the mapping file the same as its constituent samples, if all are idential. Else it will just say 'multipleValues'.


**Usage:** :file:`pool_by_metadata.py [options]`

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
		Name of otu output file, default is otu_filename.pooled.txt
	-p, `-`-map_outfile
		Name of map output file, default is map_filename.pooled.txt
	-l, `-`-pooled_sample_name
		New sample name used in new mapping file and new otu table


**Output:**

The result is a pooled OTU table and mapping file.


**Examples:**

The following command pools all the Control samples into one sample named 'pooledControl'. The resulting data is written to seqs_otus.txt.filtered.xls and Fasting_Map.txt.filtered.xls:

::

	pool_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control' -l pooledControl

Some variations are: 

Pooling all samples in both Control and Fast in the Treatment field (i.e. pooling everything):

::

	pool_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:Control,Fast' -l pooledSamples

Excluding Fast in the Treatment field - the syntax here is "*" to keep everything, then !Fast to eliminate the Fast group:

::

	pool_by_metadata.py -i seqs_otus.txt -m Fasting_Map.txt -s 'Treatment:*,!Fast' -l pooledNonFast


