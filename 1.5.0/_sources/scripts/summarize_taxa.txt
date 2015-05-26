.. _summarize_taxa:

.. index:: summarize_taxa.py

*summarize_taxa.py* -- Summarize taxa and store results in a new table or appended to an existing mapping file.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The `summarize_taxa.py <./summarize_taxa.html>`_ script provides summary information of the representation of taxonomic groups within each sample. It takes an OTU table that contains taxonomic information as input. The taxonomic level for which the summary information is provided is designated with the -L option. The meaning of this level will depend on the format of the taxon strings that are returned from the taxonomy assignment step. The taxonomy strings that are most useful are those that standardize the taxonomic level with the depth in the taxonomic strings. For instance, for the RDP classifier taxonomy, Level 2 = Domain (e.g. Bacteria), 3 = Phylum (e.g. Firmicutes), 4 = Class (e.g. Clostridia), 5 = Order (e.g. Clostridiales), 6 = Family (e.g. Clostridiaceae), and 7 = Genus (e.g. Clostridium). By default, the relative abundance of each taxonomic group will be reported, but the raw counts can be returned if -a is passed.


**Usage:** :file:`summarize_taxa.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Input OTU table filepath [REQUIRED]
	
	**[OPTIONAL]**
		
	-L, `-`-level
		Taxonomic level to summarize by. [default: 2,3,4,5,6]
	-m, `-`-mapping
		Input metadata mapping filepath. If supplied, then the taxon information will be added to this file. This option is  useful for coloring PCoA plots by taxon abundance or to  perform statistical tests of taxon/mapping associations.
	`-`-md_identifier
		The relevant observation metadat key [default: taxonomy]
	`-`-md_as_string
		Metadata is included as string [default: metadata is included as list]
	-d, `-`-delimiter
		Delimitor separating taxonomy levels. [default: ;]
	-r, `-`-relative_abundance
		DEPRECATED: please use -a/--absolute_abundance to disable relative abundance [default: ]
	-a, `-`-absolute_abundance
		If present, the absolute abundance of the lineage in  each sample is reported. By default, this script uses relative abundance [default: False]
	-l, `-`-lower_percentage
		If present, OTUs having higher absolute abundance are trimmed. To remove OTUs that make up more than 5% of the total dataset you would pass 0.05. [default: None]
	-u, `-`-upper_percentage
		If present, OTUs having lower absolute abundance are trimmed. To remove the OTUs that makes up less than 45% of the total dataset you would pass 0.45. [default: None]
	-t, `-`-transposed_output
		If present, the output will be written transposed from the regular output. This is helpful in cases when you want to use Site Painter to visualize your data [default: False]
	-o, `-`-output_dir
		Path to the output directory


**Output:**

There are two possible output formats depending on whether or not a mapping file is provided with the -m option. If a mapping file is not provided, a table is returned where the taxonomic groups are each in a row and there is a column for each sample. If a mapping file is provided, the summary information will be appended to this file. Specifically, a new column will be made for each taxonomic group to which the relative abundances or raw counts will be added to the existing rows for each sample. The addition of the taxonomic information to the mapping file allows for taxonomic coloration of Principal coordinates plots in the 3d viewer. As described in the `make_3d_plots.py <./make_3d_plots.html>`_ section, principal coordinates plots can be dynamically colored based on any of the metadata columns in the mapping file. Dynamic coloration of the plots by the relative abundances of each taxonomic group can help to distinguish which taxonomic groups are driving the clustering patterns.



**Examples:**

Summarize taxa based at taxonomic levels 2, 3, 4, 5, and 6, where the and write resulting taxa tables to the directory "./tax" 

::

	summarize_taxa.py -i otu_table.biom -o ./tax

**Examples:**

Summarize taxa based at taxonomic levels 2, 3, 4, 5, and 6, where the and write resulting mapping files to the directory "./tax" 

::

	summarize_taxa.py -i otu_table.biom -o tax_mapping/ -m Fasting_Map.txt


