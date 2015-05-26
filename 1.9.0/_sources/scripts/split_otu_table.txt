.. _split_otu_table:

.. index:: split_otu_table.py

*split_otu_table.py* -- Split a biom table into one table per combination of values found in the specified fields in the mapping file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script splits a biom table based on the cartesian product of the values
found in the mapping fields specified. It accepts any number of mapping fields
to split on. As an example assume the following was your mapping file data:

SampleID       Color       Habitat       Age
S1             Red         Stream        10
S2             Blue        Stream        20
S3             Blue        Lake          30
S4             Red         Stream        30

If we wanted to split a corresponding biom table by the 'Color' and 'Habitat'
fields simultanesouly, we would return 3 biom tables with the following samples
corresponding to the following groups:

(S1, S4): (Red, Stream)
(S2): (Blue, Stream)
(S3): (Blue, Lake)

Combinations which would result in no samples -- in our example (Red, Lake) -- 
are excluded and do not produce (empty) biom tables. The script optionally
produces split mapping files as well. 

The naming convention for split files is (assuming two fields):

input_table.biom -> input_table__field1_value1_field2_value2__.biom
input_mapping.txt -> input_mapping__field1_value1_field2_value2__.txt

So, from our example above:

input_table.biom -> (input_table__Color_Red_Habitat_Stream__.biom,
                     input_table__Color_Blue_Habitat_Stream__.biom,
                     input_table__Color_Blue_Habitat_Lake__.biom)



**Usage:** :file:`split_otu_table.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-biom_table_fp
		The input biom table file path.
	-m, `-`-mapping_fp
		The mapping file path.
	-f, `-`-fields
		Mapping columns to split biom table on, comma separated.
	-o, `-`-output_dir
		File path to the output directory to be created.
	
	**[OPTIONAL]**
		
	`-`-suppress_mapping_file_output
		Do not write out split mapping files.


**Output:**




Split otu_table.biom into per-study OTU tables, and store the results in ./per_study_otu_tables/

::

	split_otu_table.py -i otu_table.biom -m Fasting_Map.txt -f Treatment -o per_study_otu_tables

Split otu_table.biom into multiple biom tables based on the Treatment and Color of the samples

::

	split_otu_table.py -i otu_table.biom -m Fasting_Map.txt -f Treatment,Color -o ./per_study_otu_tables/


