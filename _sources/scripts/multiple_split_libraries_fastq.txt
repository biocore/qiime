.. _multiple_split_libraries_fastq:

.. index:: multiple_split_libraries_fastq.py

*multiple_split_libraries_fastq.py* -- Run split_libraries_fastq.py on multiple files.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**


This script runs `split_libraries_fastq.py <./split_libraries_fastq.html>`_ on data that are already demultiplexed
(split up according to sample, with one sample per file). The script supports
the following types of input:

- a directory containing many files, where each file is named on a per-sample
  basis (with different prefixes before the read number)
- a directory containing many directories, where each directory is named on a
  per-sample basis

This script assumes that the leading characters before the read indicator
(see --read_indicator) are matched between the read, barcode, and mapping files.
For example, sample1_L001_R1_001.fastq.gz, sample1_L001_I1_001.fastq.gz,
sample1_L001_mapping_001.txt would be matched up if "R1" is the read indicator,
"I1" is the barcode indicator, and "mapping" is the mapping file indicator.




**Usage:** :file:`multiple_split_libraries_fastq.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Input directory of directories or fastq files.
	-o, `-`-output_dir
		Output directory to write `split_libraries_fastq.py <./split_libraries_fastq.html>`_ results
	
	**[OPTIONAL]**
		
	-m, `-`-demultiplexing_method
		Method for demultiplexing. Can either be "sampleid_by_file" or "mapping_barcode_files". With the sampleid_by_file option, each fastq file (and/or directory name) will be used to generate the --sample_ids value passed to `split_libraries_fastq.py <./split_libraries_fastq.html>`_. The mapping_barcode_files option will search for barcodes and mapping files that match the input read files [default: sampleid_by_file]
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior of `split_libraries_fastq.py <./split_libraries_fastq.html>`_. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters [default: `split_libraries_fastq.py <./split_libraries_fastq.html>`_ defaults will be used]
	`-`-read_indicator
		Substring to search for to indicate read files [default: _R1_]
	`-`-barcode_indicator
		Substring to search for to indicate barcode files [default: _I1_]
	`-`-mapping_indicator
		Substring to search for to indicate mapping files [default: _mapping_]
	`-`-mapping_extensions
		Comma-separated list of file extensions used to identify mapping files. Only applies when --demultiplexing_method is "mapping_barcode_files" [default: txt,tsv]
	`-`-sampleid_indicator
		Text in fastq filename before this value will be used as output sample ids [default: _]
	`-`-include_input_dir_path
		Include the input directory name in the output sample id name. Useful in cases where the file names are repeated in input folders [default: False]
	`-`-remove_filepath_in_name
		Disable inclusion of the input filename in the output sample id names. Must use --include_input_dir_path if this option is enabled [default: False]
	`-`-leading_text
		Leading text to add to each `split_libraries_fastq.py <./split_libraries_fastq.html>`_ command [default: no leading text added]
	`-`-trailing_text
		Trailing text to add to each `split_libraries_fastq.py <./split_libraries_fastq.html>`_ command [default: no trailing text added]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]


**Output:**

The output of running `split_libraries_fastq.py <./split_libraries_fastq.html>`_ on many input files. See script description for more details.


**Example 1:**

Process an input folder of folders, with options specified to pair up reads, barcodes, and mapping files. A qiime_parameters.txt file is included to use the parameter split_libraries_fastq:barcode_type 12

::

	multiple_split_libraries_fastq.py -i input_folders -o output_folder --demultiplexing_method mapping_barcode_files --read_indicator reads --barcode_indicator barcode --mapping_indicator mapping -p qiime_parameters.txt

**Example 2:**

Process an input folder of files, with the option specified to generate sample ids using the filenames (default behavior is to use all text before the first underscore as the sample id)

::

	multiple_split_libraries_fastq.py -i input_files -o output_folder --demultiplexing_method sampleid_by_file

**Example 3:**

Process an input folder of folders, with an option specified to use the folder names as the sample ids. In this case, the fastq filenames themselves are not included, only the folder names are used.

::

	multiple_split_libraries_fastq.py -i input_folders_no_barcodes -o output_folder --demultiplexing_method sampleid_by_file --include_input_dir_path --remove_filepath_in_name

**Example 4:**

To see what commands would be executed by the script without actually running them, use the following command:

::

	multiple_split_libraries_fastq.py -i input_files -o output_folder -w


