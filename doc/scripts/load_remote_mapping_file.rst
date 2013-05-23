.. _load_remote_mapping_file:

.. index:: load_remote_mapping_file.py

*load_remote_mapping_file.py* -- Downloads and saves a remote mapping file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script exports, downloads, and saves a mapping file that is stored
remotely. Currently, the only type of remote mapping file that is supported is
a Google Spreadsheet, though other methods of remote storage may be supported
in the future.

For more information and examples pertaining to this script and remote mapping
files in general, please refer to the accompanying tutorial, which can be found
at http://qiime.org/tutorials/remote_mapping_files.html.



**Usage:** :file:`load_remote_mapping_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-k, `-`-spreadsheet_key
		The spreadsheet key that will be used to identify the Google Spreadsheet to load. This is the part of the Google Spreadsheet URL that comes after 'key='. You may instead provide the entire URL and the key will be extracted from it. If you provide the entire URL, you may need to enclose it in single quotes
	-o, `-`-output_fp
		The output filepath
	
	**[OPTIONAL]**
		
	-w, `-`-worksheet_name
		The name of the worksheet in the Google Spreadsheet that contains the mapping file. If the worksheet name contains spaces, please include quotes around the name. [default: the first worksheet in the Google Spreadsheet will be used]


**Output:**


The script outputs a single file, which is the metadata mapping file obtained
from the remote location (in QIIME-compatible format).



**Load mapping file from Google Spreadsheet:**

The following command exports and downloads a QIIME metadata mapping file from a Google Spreadsheet, using the data found in the first worksheet of the spreadsheet.

::

	load_remote_mapping_file.py -k 0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc -o example1_map.txt

**Load specific worksheet:**

The following command exports from a worksheet named 'Fasting_Map'.

::

	load_remote_mapping_file.py -k 0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc -w Fasting_Map -o example2_map.txt


