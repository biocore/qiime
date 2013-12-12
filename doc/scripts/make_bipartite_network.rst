.. _make_bipartite_network:

.. index:: make_bipartite_network.py

*make_bipartite_network.py* -- This script makes a bipartite network connecting samples to observations. It is most suitable for visualization with cytoscape.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script was created to ease the process of making bipartite networks that have appeared in high profile publications including 10.1073/pnas.1217767110 and 10.1126/science.1198719. The script will take an biom table and a mapping file and produce an edge table which connects each sample in the biom table to the observations found in that sample. It is bipartite because there are two distinct node classes -- OTU and Sample. The 'OTU' node class does not have to be an operational taxonomic unit, it can be a KEGG category or metabolite etc. -- anything that is an observation. The edges are weighted by the abundance of the observation in the sample to which it is connected. The output files of this script are intended to be loaded into Cytoscape. The EdgeTable should be uploaded first, and then the NodeAttrTable file can be uploaded as node attributes to control coloring, sizing, and shaping as the user desires. The overall idea behind this script is to make bipartite network creation easier. To that end, the color, size, and shape options are used to provide fields in the NetworkViz tab of Cytoscape so that nodes can be appropriately presented. Those options are passed via comma separated strings (as in the example below). The most common visualization strategy is to color sample nodes by a metadata category like timepoint or pH, color OTU nodes by one of their taxonomic levels, and to scale OTU node size by abundance. This script makes this process easy (as well as a myriad of other visualiation strategies). Once the tables are created by this script they must be opened in cytoscape. This process is described in detail in the QIIME bipartite network tutorial available at: http://qiime.org/tutorials/making_cytoscape_networks.html All color, size, and shape options in this script default to 'NodeType'. OTU nodes have NodeType: otu, sample nodes have NodeType: sample. Thus, if you ran this script with defaults, you would only be able to change the shape, size, and color of the nodes depending on whether or not they were observations or samples. You would not be able to distinguish between two observations based on color, shape, or size. The script is flexible in that it allows you to pass any number of fields for the --{s,o}{shape,size,color}. This will allow you to distinguish between OTU and sample nodes in a huge number of different ways. The usage examples below show some of the common use cases and what options you would pass to emulate them. There are a couple of important considerations for using this script:

Note that the --md_fields option has a different meaning depending on the type of metadata in the biom table. Regardless of type, the md_fields will be the headers in the OTUNodeTable.txt. If the metadata is a dict or default dict, the md_fields will be used as keys to extract data from the biom file metadata. If the metadata is a list or a string, then the md_fields will be have no intrinsic relation to the columns they head. For example if md_fields=['k','p','c'] and the metadata contained in a given OTU was 'k__Bacteria;p__Actinobacter;c__Actino' the resulting OTUNodeTable would have k__Bacteria in the 'k' column, p__Actinobacter in the 'p' column, and c__Actino in the 'c' column. If one passed md_fields=['1.0','XYZ','Five'] then the OTUNodeTable would have columns headed by ['1.0','XYZ','Five'], but the metadata values in those columns would be the same (e.g. '1.0' column entry would be k__Bacteria etc.) If the number of elements in the metadata for a given OTU is not equal to the number of headers provided the script will adjust the OTU metadata. In the case where the metadata is too short, it will add 'Other' into the OTU metadata until the required length is reached.  In the case where the metadata is too long it will simply remove extra entries. This means you can end up with many observations which have the value of 'Other' if you have short taxonomic strings/lists for your observations.

The available fields for both sample and otu nodes are:
[NodeType, Abundance]

For observation nodes the additional fields available are:
any fields you passed for the md_fields

For sample nodes the additional fields available are
any fields found in the mapping file headers

If multiple fields are passed for a given option, they will be concatenated in the output with a '_' character. 



**Usage:** :file:`make_bipartite_network.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-biom_fp
		The input file path for biom table.
	-m, `-`-map_fp
		The input file path for mapping file.
	-o, `-`-output_dir
		Directory that will be created for storing the results.
	-k, `-`-observation_md_header_key
		Key to retrieve metadata (usually taxonomy) from the biom file.
	`-`-md_fields
		Metadata fields that will be the headers of the OTUNodeTable. If the biom table has metadata dictionaries, md_fields will be the keys extracted from the biom table metadata. Passed like "kingdom,phylum,class".
	
	**[OPTIONAL]**
		
	`-`-scolors
		Commas separated string specifying fields of interest for sample node coloring [default: NodeType].
	`-`-ocolors
		Commas separated string specifying fields of interest for observation node coloring [default: NodeType].
	`-`-sshapes
		Commas separated string specifying fields of interest for sample node shape [default: NodeType].
	`-`-oshapes
		Commas separated string specifying fields of interest for observation node shape [default: NodeType].
	`-`-ssizes
		Commas separated string specifying fields of interest for sample node size [default: NodeType].
	`-`-osizes
		Commas separated string specifying fields of interest for observation node size [default: NodeType].


**Output:**

The output of this script is four files:
    1. EdgeTable - table with connections between samples and observations. 
    2. OTUNodeTable - table with observations and their associated metadata. 
    3. SampleNodeTable - table with samples and their associated metadata.
    4. NodeAttrTable - table with the node attributes specified by the user with
       the given options.


Create an EdgeTable and NodeAttrTable that allow you to color sample nodes with one of their metadata categories (Treatment for our example), observation nodes (in this case OTUs) by their taxonomic level (class for our example), control observation node size by their abundance, and control node shape by whether its an observation or sample.

::

	make_bipartite_network.py -i otu_table.biom -m mapping_file.txt -k taxonomy --md_fields 'k,p,c,o,f' -o bipartite_network/ --scolors 'Treatment' --ocolors 'c' --osize 'Abundance'

Create an EdgeTable and NodeAttrTable that allow you to color sample nodes by a combination of their time point and diet, color observation nodes by their abundance and family, and node shape by whether the node is an observation or sample. Note that the names in the --md_fields are irrelevant as long as the field passed for --ocolors is available. The length is important however, since there are 5 levels in our OTU table. If fewer fewer than 5 fields were passed for --md_fields we would get an error.

::

	make_bipartite_network.py -i otu_table.biom -m mapping_file.txt -k taxonomy --md_fields 'a1,a2,a3,a4,f' -o bipartite_network_combo_colors/ --scolors 'TimePt,Diet' --ocolors 'f,Abundance'


