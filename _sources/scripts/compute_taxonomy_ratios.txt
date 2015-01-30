.. _compute_taxonomy_ratios:

.. index:: compute_taxonomy_ratios.py

*compute_taxonomy_ratios.py* -- Compute the log ratio abundance of specified taxonomic groups.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Compute the log ratio abundance of specified taxonomic groups. This method is based on the microbial dysbiosis index described in Gevers et al. 2014: http://www.ncbi.nlm.nih.gov/pubmed/24629344


**Usage:** :file:`compute_taxonomy_ratios.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input
		The input BIOM table [REQUIRED if not passing -s]
	-o, `-`-output
		Path to where the output will be written; this will be a new sample metadata mapping file [REQUIRED if not passing -s]
	`-`-increased
		Comma-separated list of taxa whose abundances are included in the numerator of the ratio [REQUIRED if not passing -s or -e]
	`-`-decreased
		Comma-separated list of taxa whose abundances are included in the denominator of the ratio [REQUIRED if not passing -s or -e]
	-e, `-`-index
		Apply an existing index. Options are: md [REQUIRED if not passing -s or --increased and --decreased]
	-n, `-`-name
		Column name for the index in the output file [default: 'index', or value passed as -e if provided]
	-m, `-`-mapping_file
		A mapping file containing data that should be included in the output file [default: no additional mapping file data is included in output]
	-k, `-`-key
		Metadata key to use for computing index [default: taxonomy]
	-s, `-`-show_indices
		List known indices and exit [default: False]


**Output:**

By default, a minimal QIIME mapping file is created containing two columns: SampleID and the index. If -m is provided, the information in that mapping file is merged into the default output mapping file.


**Example:**

Compute the microbial dysbiosis (MD) index

::

	compute_taxonomy_ratios.py -i table.biom.gz -e md -o md.txt

**Example:**

Compute the microbial dysbiosis (MD) index and add it to an existing mapping file

::

	compute_taxonomy_ratios.py -i table.biom.gz -e md -o map_w_md.txt -m map.txt

**Example:**

Compute the log of the abundance of p__Firmicutes plus p__Fusobacteria divided by the abundance of p__Bacteroidetes and write the results to custom_index.txt.

::

	compute_taxonomy_ratios.py -i table.biom.gz --increased p__Firmicutes,p__Fusobacteria --decreased p__Bacteroidetes -o custom_index.txt


