.. _summarize_taxa_through_plots:

.. index:: summarize_taxa_through_plots.py

*summarize_taxa_through_plots.py* -- A workflow script for performing taxonomy summaries and plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


The steps performed by this script are:

1. Summarize OTU by Category

2. Summarize Taxonomy

3. Plot Taxonomy Summary




**Usage:** :file:`summarize_taxa_through_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input otu table [REQUIRED]
	-m, `-`-mapping_fp
		Path to the mapping file [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters. [if omitted, default values will be used]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-c, `-`-mapping_category
		Summarize OTU table using this category. [default: None]


**Output:**

The results of this script is a folder ("wf_taxa_sum/") containing taxonomy summary files (at different levels) and a folder containing taxonomy summary plots. Additionally, if a mapping_catgory is supplied there will be a summarized OTU table.


**Examples:**

::

	summarize_taxa_through_plots.py -o wf_taxa_sum -i otu_table.txt -m inseqs1_mapping.txt -p custom_parameters.txt

Alternatively, the user can supply a mapping_category, where the OTU is summarized based on a mapping category:

::

	summarize_taxa_through_plots.py -o wf_taxa_sum -i otu_table.txt -m inseqs1_mapping.txt -p custom_parameters.txt -c Treatment


