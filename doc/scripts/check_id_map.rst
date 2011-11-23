.. _check_id_map:

.. index:: check_id_map.py

*check_id_map.py* -- Checks user's metadata mapping file for required data, valid format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Specifically, we check that:

    - The filename does not contain spaces (warn + rewrite if it does)
    - There are headers for SampleID, LinkerPrimerSequence, and BarcodeSequence if barcodes are used (returns errors if these are absent or misspelled)
    - The BarcodeSequence and LinkerPrimerSequences fields have valid IUPAC DNA characters
    - There are not duplicate header fields (error)
    - There are not duplicate near-unique but not exactly unique values within each column (warning)
    - The headers do not contain invalid characters (alphanumeric and underscore only)
    - The data fields do not contain invalid characters (alphanumeric, underscore, space, and +-%./:,; characters)
    - SampleID fields are MIENS compliant (only alphanumeric and . characters)
    - There are no duplicates when the primer and barcodes are appended
    - If there is a field ReversePrimer for reverse primers (for removal with split_libraries), the characters are DNA IUPAC compliant and no fields are empty
    
    Errors and warnings are saved to a log file.  Errors are generally caused 
    by problems with the headers, and should be resolved before attempting to 
    correct any warnings.  Warnings can arise from invalid characters, 
    near-duplicate metadata, duplicate sample descriptions/barcodes, or missing
    data fields. Warnings will contain a reference to the cell (row,column) 
    that the warning arose from.
    
    In addition to the log file, a "corrected_mapping" file will be created.
    Invalid characters will be replaced by underscores in this corrected mapping
    file if there were any such characters in the input metadata mapping file.
    If there were no invalid characters to replace, the corrected mapping file 
    will contain comments saying as much.
    
    `check_id_map.py <./check_id_map.html>`_ should not raise exceptions itself under normal 
    circumstances, except for situations such as having a misformatted input 
    metadata mapping file.
    
    If pooled primers are used, separate with a comma.  For instance, a pooled
    set of three 27f primers (used to increase taxonomic coverage) could be
    specified in the LinkerPrimerSequence fields as such:
    AGGGTTCGATTCTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTGGTTCAG



**Usage:** :file:`check_id_map.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-map
		Metadata mapping file filepath
	-o, `-`-output_dir
		Required output directory for log file and corrected mapping file (by default, invalid characters will be converted to underscores)
	
	**[OPTIONAL]**
		
	-c, `-`-char_replace
		Changes the default character used to replace invalid characters found in the mapping file.  Must be a valid character (alphanumeric or underscore).  NOT IMPLEMENTED CURRENTLY [default: _]
	-b, `-`-not_barcoded
		Use -b if barcodes are not present. [default: False]
	-B, `-`-variable_len_barcodes
		Use -B if variable length barcodes are present to suppress warnings about barcodes of unequal length. [default: False]
	-p, `-`-disable_primer_check
		Use -p to disable checks for primers. [default: False]
	-v, `-`-verbose
		Turn on this flag to disable verbose output.  [default: True]
	-j, `-`-added_demultiplex_field
		Use -j to add a field to use in the mapping file as an additional demultiplexing option to the barcode.  All combinations of barcodes and the values in these fields must be unique. The fields must contain values that can be parsed from the fasta labels such as "plate=R_2008_12_09".  In this case, "plate" would be the column header and "R_2008_12_09" would be the field data (minus quotes) in the mapping file.  To use the run prefix from the fasta label, such as ">FLP3FBN01ELBSX", where "FLP3FBN01" is generated from the run ID, use "-j run_prefix" and set the run prefix to be used as the data under the column headerr "run_prefix".  [default: None]


**Output:**

A log file and corrected_mapping.txt file will be written to the mapping_info directory.


**Example:**

Check the test_mapping.txt mapping file for problems, supplying the required mapping file and output directory (in this case mapping_info)

::

	check_id_map.py -m test_mapping.txt -o mapping_info/


