.. _summarize_taxa:

.. index:: summarize_taxa.py

*summarize_taxa.py* -- Summarize Taxa
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The `summarize_taxa.py <./summarize_taxa.html>`_ script provides summary information of the representation of taxonomic groups within each sample. It takes an OTU table that contains taxonomic information as input. The taxonomic level for which the summary information is provided is designated with the -L option. The meaning of this level will depend on the format of the taxon strings that are returned from the taxonomy assignment step. The taxonomy strings that are most useful are those that standardize the taxonomic level with the depth in the taxonomic strings. For instance, for the RDP classifier taxonomy, Level 2 = Domain (e.g. Bacteria), 3 = Phylum (e.g. Firmicutes), 4 = Class (e.g. Clostridia), 5 = Order (e.g. Clostridiales), 6 = Family (e.g. Clostridiaceae), and 7 = Genus (e.g. Clostridium). By default, the relative abundance of each taxonomic group will be reported, but the raw counts can be returned if -a is passed.


**Usage:** :file:`summarize_taxa.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_file
		Path to read otu file [REQUIRED]
	
	**[OPTIONAL]**
		
	-o, `-`-output_file
		Path to write output file
	-L, `-`-level
		Level of taxonomy to use [default: 2]
	-m, `-`-mapping
		If supplied - the taxon information will be added to the mapping file. This mapping file can be used to color PCoA plots by taxon abundance or to perform statistical tests of taxon/mappingy associations.
	-d, `-`-delimiter
		Delimitor that separates taxonomy categories.[default: ;]
	-r, `-`-relative_abundance
		DEPRECATED: please use -a/--absolute_abundance to disable relative abundance [default: ]
	-a, `-`-absolute_abundance
		If present, reports the absolute abundance of the lineage in each sample. By default uses relative abundance [default: False]


**Output:**

There are two possible output formats depending on whether or not a mapping file is provided with the -m option. If a mapping file is not provided, a table is returned where the taxonomic groups are each in a row and there is a column for each sample. If a mapping file is provided, the summary information will be appended to this file. Specifically, a new column will be made for each taxonomic group to which the relative abundances or raw counts will be added to the existing rows for each sample. The addition of the taxonomic information to the mapping file allows for taxonomic coloration of Principal coordinates plots in the 3d viewer. As described in the `make_3d_plots.py <./make_3d_plots.html>`_ section, principal coordinates plots can be dynamically colored based on any of the metadata columns in the mapping file. Dynamic coloration of the plots by the relative abundances of each taxonomic group can help to distinguish which taxonomic groups are driving the clustering patterns.



**Examples:**

The following command can be used to summarize taxa based on the Class, where the default parameters are used (no mapping file, delimiter for RDP ("-d ;") and output relative abundance) and the results are written to the file "Class.txt":

::

	summarize_taxa.py -i otu_table.txt -L 4 -o Class.txt

Optionally the user can have the relative abundances added to the user-generated mapping file, by using the following command:

::

	summarize_taxa.py -i otu_table.txt -L 4 -m Mapping_file.txt

Alternatively, the user may want to output the raw counts of each lineage within a sample, which can be used in the next step for making pie charts, by using the following command:

::

	summarize_taxa.py -i otu_table.txt -L 4 -a


