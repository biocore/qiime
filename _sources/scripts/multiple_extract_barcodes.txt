.. _multiple_extract_barcodes:

.. index:: multiple_extract_barcodes.py

*multiple_extract_barcodes.py* -- Run extract_barcodes.py on multiple files.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script runs `extract_barcodes.py <./extract_barcodes.html>`_ on data that are already demultiplexed
(split up according to sample, with one sample per file). The script
supports the following types of input:

- a directory containing many files, where each file is named on a per-sample
  basis
- a directory containing many directories, where each directory is named on a
  per-sample basis
 
The script assumes that the leading/trailing characters before/after the read
number indicator (see --read1_indicator) are matched between forward and
reverse reads. For example:

- S0_L001_R1_001.fastq.gz and S0_L001_R2_001.fastq.gz would be matched up reads
- S0_L002_R1_00X.fastq.gz and S0_L002_R2_00X.fastq.gz would be matched up reads

The output directory used for each call to `extract_barcodes.py <./extract_barcodes.html>`_ uses the base
name of the input read 1 fastq file (a single directory would be problematic
since the output names for `extract_barcodes.py <./extract_barcodes.html>`_ can be the same for different
calls). Use the parameter --include_input_dir_path to also include the input
directory name in the output directory path, which may be preferable in the
case of an input folder of folders, and --remove_filepath_in_name can be used
in this case to prevent the input read 1 fastq file base name from being used
as part of the output directory name.




**Usage:** :file:`multiple_extract_barcodes.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Input directory of directories, or directory of paired fastq files.
	-o, `-`-output_dir
		Base output directory to write output folders
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior of `extract_barcodes.py <./extract_barcodes.html>`_. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters [default: `extract_barcodes.py <./extract_barcodes.html>`_ defaults will be used]
	`-`-paired_data
		Turn this option on if paired data are to be used. The type of paired data for `extract_barcodes.py <./extract_barcodes.html>`_ should be specified with -p. Forward and reverse reads will be searched for via the --read1_indicator and --read2_indicator parameters [default: False]
	`-`-read1_indicator
		Substring to search for to indicate read 1 [default: _R1_]
	`-`-read2_indicator
		Substring to search for to indicate read 2 [default: _R2_]
	`-`-leading_text
		Leading text to add to each `extract_barcodes.py <./extract_barcodes.html>`_ command [default: no leading text added]
	`-`-trailing_text
		Trailing text to add to each `extract_barcodes.py <./extract_barcodes.html>`_ command [default: no trailing text added]
	`-`-include_input_dir_path
		Include the input directory name in the output directory path. Useful in cases where the file names are repeated in input folders [default: False]
	`-`-remove_filepath_in_name
		Disable inclusion of the input filename in the output directory names. Must use --include_input_dir_path if this option is enabled [default: False]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]


**Output:**

The output of running `extract_barcodes.py <./extract_barcodes.html>`_ on many input files. See script description for more details.


**Example 1:**

Process an input folder of files, with default options used for `extract_barcodes.py <./extract_barcodes.html>`_:

::

	multiple_extract_barcodes.py -i input_files -o output_folder

**Example 2:**

Process an input folder of folders (with the filenames having _forward_ and _reverse_ containing the forward and reverse read filenames, respectively), using the `extract_barcodes.py <./extract_barcodes.html>`_ option for paired reads. The individual folder names are included in the output folder names, but not the filenames. Note: it is important to pass the --paired_data option if paired data are to be used in the `extract_barcodes.py <./extract_barcodes.html>`_ commands. Additionally, the paired fastq file type for `extract_barcodes.py <./extract_barcodes.html>`_ is specified with a qiime_parameters.txt file (with this value specified: extract_barcodes:input_type barcode_paired_end)

::

	multiple_extract_barcodes.py -i input_folders -o output_folder -p qiime_parameters.txt --paired_data --read1_indicator '_forward_' --read2_indicator '_reverse_' --include_input_dir_path --remove_filepath_in_name

**Example 3:**

To see what commands would be executed by the script without actually running them, use the following command:

::

	multiple_extract_barcodes.py -i input_files -o output_folder -w


