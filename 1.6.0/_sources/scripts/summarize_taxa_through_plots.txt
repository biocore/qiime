.. _summarize_taxa_through_plots:

.. index:: summarize_taxa_through_plots.py

*summarize_taxa_through_plots.py* -- A workflow script for performing taxonomy summaries and plots
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


The steps performed by this script are: Summarize OTU by Category (optional, pass -c); Summarize Taxonomy; and Plot Taxonomy Summary


**Usage:** :file:`summarize_taxa_through_plots.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		The input otu table [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters. [if omitted, default values will be used]
	-m, `-`-mapping_fp
		Path to the mapping file [REQUIRED if passing -c]
	-f, `-`-force
		Force overwrite of existing output directory (note: existing files in output_dir will not be removed) [default: None]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]
	-c, `-`-mapping_category
		Summarize OTU table using this category. [default: None]
	-s, `-`-sort
		Sort the OTU Table [default: False]


**Output:**

The results of this script is a folder (specified by -o) containing taxonomy summary files (at different levels) and a folder containing taxonomy summary plots. Additionally, if a mapping_catgory is supplied there will be a summarized OTU table. The primary interface for this output are the OUTPUT_DIR/taxa_summary_plots/\*html files which are interactive plots that can be opened in a web browser (see the mouse-overs for interactivity).


**Plot taxa summaries for all samples:**

::

	summarize_taxa_through_plots.py -o taxa_summary -i otu_table.biom -m Fasting_Map.txt

**Plot taxa summaries on a categorical basis:**

Alternatively, the user can supply a mapping_category, where the OTU is summarized based on a sample metadata category:

::

	summarize_taxa_through_plots.py -o taxa_summary_by_treatment -i otu_table.biom -m Fasting_Map.txt -c Treatment


