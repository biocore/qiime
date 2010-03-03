.. _sra_spreadsheet_to_map_files:

.. index:: sra_spreadsheet_to_map_files.py

*sra_spreadsheet_to_map_files.py* -- Create mapping file from SRA submission spreadsheet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script reads an SRA submission spreadsheet and generates QIIME mapping files.


**Usage:** :file:`sra_spreadsheet_to_map_files.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_file
		The input SRA submission spreadsheet


**Output:**

Produces one map file per (STUDY, RUN_PREFIX) combination. Note that the output will include extra stuff not actually needed by QIIME. The intention is just to pull out the info needed for `split_libaries.py <./split_libaries.html>`_ and downstream analyses. Currently, this does not combine this with the data in the per-sample mapping file.


**Simple example:**

Take an SRA submission spreadsheet input_spreadsheet.txt and write out map files as a series of files input_spreadsheet_[STUDY].txt.map.

::

	sra_spreadsheet_to_map_files.py -i input_spreadsheet.txt


