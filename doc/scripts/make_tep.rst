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
		
	-i, `-`-otu_file
		Path to read otu table
	-m, `-`-mapping_file
		Path to read mapping file
	-t, `-`-tree_file
		Path to read tree
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		The output directory [default: None]
	-w, `-`-web
		Web codebase jnlp flag [default: False]
	-u, `-`-url_path
		Url path


**Output:**

The result of this script is written to a .tep file and a .jnlp file, both with the name supplied by -o


**Example:**

Create .tep file and .jnlp file:

::

	make_tep.py -i otu_table.txt -m mapping_file.txt -t otus.tre -o my_data


