.. _make_tep:

.. index:: make_tep.py

*make_tep.py* -- Makes TopiaryExplorer project file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script makes a TopiaryExplorer project file (.tep) and a jnlp file with the data location preloaded.

WARNING: The jnlp file relies on an absolute path, if you move the .tep file, the generated jnlp will no longer work. However, you can still open the .tep file from your normal TopiaryExplorer install.


**Usage:** :file:`make_tep.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-otu_table_fp
		Path to otu table in biom format
	-m, `-`-mapping_fp
		Path to mapping file
	-t, `-`-tree_fp
		Path to tree
	
	**[OPTIONAL]**
		
	-o, `-`-out_fp
		The output directory [default: None]
	-p, `-`-prefs_file_fp
		Path to prefs file
	-w, `-`-web_flag
		Web codebase jnlp flag [default: False]
	-u, `-`-url
		Url path


**Output:**

The result of this script is written to a .tep file and a .jnlp file, both with the name supplied by -o


**Example:**

Create .tep file and .jnlp file:

::

	make_tep.py -i otu_table.biom -m Fasting_Map.txt -t rep_set.tre


