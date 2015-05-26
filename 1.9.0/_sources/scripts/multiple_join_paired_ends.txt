.. _multiple_join_paired_ends:

.. index:: multiple_join_paired_ends.py

*multiple_join_paired_ends.py* -- Run join_paired_ends.py on multiple files.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script runs `join_paired_ends.py <./join_paired_ends.html>`_ on data that are already demultiplexed
(split up according to sample, with one sample per pair of files). The script
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

If an optional --barcode_indicator file is used, it is searched for in the same
manner that the paired files are searched for, so if the default "_I1_" is
used, S0_L001_R1_001.fastq.gz and S0_L001_R2_001.fastq.gz would be matched up
with S0_L001_I1_001.fastq.gz as the barcode indicator file.

The output directory used for each call to `join_paired_ends.py <./join_paired_ends.html>`_ uses the base
name of the input read 1 fastq file (a single directory would be problematic
since the output names for `join_paired_ends.py <./join_paired_ends.html>`_ can be the same for different
calls). Use the parameter --include_input_dir_path to also include the input
directory name in the output directory path, which may be preferable in the
case of an input folder of folders, and --remove_filepath_in_name can be used
in this case to prevent the input read 1 fastq file base name from being used
as part of the output directory name.




**Usage:** :file:`multiple_join_paired_ends.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_dir
		Input directory of directories, or directory of paired fastq files.
	-o, `-`-output_dir
		Base output directory to write output folders
	
	**[OPTIONAL]**
		
	-p, `-`-parameter_fp
		Path to the parameter file, which specifies changes to the default behavior of `join_paired_ends.py <./join_paired_ends.html>`_. See http://www.qiime.org/documentation/file_formats.html#qiime-parameters [default: `join_paired_ends.py <./join_paired_ends.html>`_ defaults will be used]
	`-`-read1_indicator
		Substring to search for to indicate read 1 [default: _R1_]
	`-`-read2_indicator
		Substring to search for to indicate read 2 [default: _R2_]
	-b, `-`-match_barcodes
		Enable searching for matching barcodes [default: False]
	`-`-barcode_indicator
		Substring to search for to indicate barcode reads [default: _I1_]
	`-`-leading_text
		Leading text to add to each `join_paired_ends.py <./join_paired_ends.html>`_ command [default: no leading text added]
	`-`-trailing_text
		Trailing text to add to each `join_paired_ends.py <./join_paired_ends.html>`_ command [default: no trailing text added]
	`-`-include_input_dir_path
		Include the input directory name in the output directory path. Useful in cases where the file names are repeated in input folders [default: False]
	`-`-remove_filepath_in_name
		Disable inclusion of the input filename in the output directory names. Must use --include_input_dir_path if this option is enabled [default: False]
	-w, `-`-print_only
		Print the commands but don't call them -- useful for debugging [default: False]


**Output:**

The output of running `join_paired_ends.py <./join_paired_ends.html>`_ on many input files. See script description for more details.


**Example 1:**

Process an input folder of paired-up files (by filename, with the default _R1_ and _R2_ containing the forward and reverse reads filenames, respectively). An optional parameters file is passed with -p. This file can specify an optional parameter for `join_paired_ends.py <./join_paired_ends.html>`_, such as: join_paired_ends:pe_join_method SeqPrep

::

	multiple_join_paired_ends.py -i input_files -o output_folder -p qiime_parameters.txt

**Example 2:**

Process an input folder of folders (with the filenames having _forward_ and _reverse_ containing the forward and reverse read filenames, respectively). The individual folder names are included in the output folder names, but not the filenames. A matching barcode fastq file (indicated by _barcode_) is also included.

::

	multiple_join_paired_ends.py -i input_folders -o output_folder --read1_indicator '_forward_' --read2_indicator '_reverse_' --include_input_dir_path --remove_filepath_in_name -b --barcode_indicator '_barcode_'

**Example 3:**

To see what commands would be executed by the script without actually running them, use the following command:

::

	multiple_join_paired_ends.py -i input_files -o output_folder -w


