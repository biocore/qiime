.. _validate_mapping_file:

.. index:: validate_mapping_file.py

*validate_mapping_file.py* -- Checks user's metadata mapping file for required data, valid format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Specifically, we check that:

    - The BarcodeSequence, LinkerPrimerSequences, and ReversePrimer fields 
       have valid IUPAC DNA characters, and BarcodeSequence characters
       are non-degenerate (error)
    - The SampleID, BarcodeSequence, LinkerPrimerSequence, and Description
       headers are present. (error)
    - There are not duplicate header fields (error)
    - There are not duplicate barcodes (error)
    - Barcodes are of the same length.  Suppressed when 
       variable_len_barcode flag is passed (warning)
    - The headers do not contain invalid characters (alphanumeric and 
       underscore only) (warning)
    - The data fields do not contain invalid characters (alphanumeric, 
       underscore, space, and +-%./:,; characters) (warning)
    - SampleID fields are MIENS compliant (only alphanumeric
       and . characters). (warning)
    - There are no duplicates when the primer and variable length 
       barcodes are appended (error)
    - There are no duplicates when barcodes and added demultiplex 
       fields (-j option) are combined (error)
    - Data fields are not found beyond the Description column (warning)
      
    Details about the metadata mapping file format can be found here:
    http://www.qiime.org/documentation/file_formats.html#metadata-mapping-files
    
    Errors and warnings are saved to a log file.  Errors can be caused by
    problems with the headers, invalid characters in barcodes or primers, or
    by duplications in SampleIDs or barcodes.
    
    Warnings can arise from invalid characters and variable length barcodes that
    are not specified with the --variable_len_barcode.
    Warnings will contain a reference to the cell (row,column) that the 
    warning arose from.
    
    In addition to the log file, a "corrected_mapping" file will be created.
    Any invalid characters will be replaced with '.' characters in
    the SampleID fields (to enforce MIENS compliance) and text in other data
    fields will be replaced with the character specified by the -c parameter,
    which is an underscore "_" by default.
    
    A html file will be created as well, which will show locations of 
    warnings and errors, highlighted in yellow and red respectively.  If no
    errors or warnings were present the file will display a message saying 
    such.  Header errors can mask other errors, so these should be corrected
    first.
    
    If pooled primers are used, separate with a comma.  For instance, a pooled
    set of three 27f primers (used to increase taxonomic coverage) could be
    specified in the LinkerPrimerSequence fields as such:
    AGGGTTCGATTCTGGCTCAG,AGAGTTTGATCCTGGCTTAG,AGAATTTGATCTTGGTTCAG



**Usage:** :file:`validate_mapping_file.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-m, `-`-mapping_fp
		Metadata mapping filepath
	
	**[OPTIONAL]**
		
	-o, `-`-output_dir
		Required output directory for log file, corrected mapping file, and html file. [default: ./]
	-v, `-`-verbose
		Enable printing information to standard out [default: True]
	-c, `-`-char_replace
		Changes the default character used to replace invalid characters found in the mapping file.  Must be a valid character (alphanumeric, period, or underscore).[default: _]
	-b, `-`-not_barcoded
		Use -b if barcodes are not present.  BarcodeSequence header still required.  [default: False]
	-B, `-`-variable_len_barcodes
		Use -B if variable length barcodes are present to suppress warnings about barcodes of unequal length. [default: False]
	-p, `-`-disable_primer_check
		Use -p to disable checks for primers.  LinkerPrimerSequence header still required. [default: False]
	-j, `-`-added_demultiplex_field
		Use -j to add a field to use in the mapping file as additional demultiplexing (can be used with or without barcodes).  All combinations of barcodes/primers and the these fields must be unique. The fields must contain values that can be parsed from the fasta labels such as "plate=R_2008_12_09".  In this case, "plate" would be the column header and "R_2008_12_09" would be the field data (minus quotes) in the mapping file.  To use the run prefix from the fasta label, such as ">FLP3FBN01ELBSX", where "FLP3FBN01" is generated from the run ID, use "-j run_prefix" and set the run prefix to be used as the data under the column header "run_prefix".  [default: None]
	-s, `-`-suppress_html
		Use -s to disable html file generation, can be useful for extremely large mapping files. [default: False]


**Output:**

A log file, html file, and corrected_mapping.txt file will be written to the current output directory.


**Example:**

Check the Fasting_Map.txt     mapping file for problems, supplying the required mapping file, and output     the results in the validate_mapping_file_output directory

::

	validate_mapping_file.py -m     Fasting_Map.txt -o validate_mapping_file_output


