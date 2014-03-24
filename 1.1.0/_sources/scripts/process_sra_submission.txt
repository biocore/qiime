.. _process_sra_submission:

.. index:: process_sra_submission.py

*process_sra_submission.py* -- A workflow script for creating a second-stage SRA submission.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The steps performed by this script are:

1. Get fasta and qual from sff files.

2. Produce valid mapping file for library demultiplexing.

3. Demultiplex libraries.

4. Reduce sequence complexity by picking OTUs with cd-hit.

5. Pick a representative sequence for each OTU.

6. Blast the representative set sequences against 95% OTUs in greengenes to eliminate non-16S  sequences.

7. Make per-library files of "good" ids to pass to sfffile.

8. Use sfffile to make per-library sff files.

9. Use sfffile to quality-trim the barcodes, primers and linkers.

10. Move files around and make archive.

11. Finally, make the XML files for a second-stage submission.


**Usage:** :file:`process_sra_submission.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-u, `-`-input_submission_fp
		The tab-delimited text file with info about the submission [REQUIRED]
	-s, `-`-sff_dir
		Directory containing sff files [REQUIRED]
	-e, `-`-input_experiment_fp
		The tab-delimited text file with info about the experiment [REQUIRED]
	-r, `-`-reference_set_fp
		Path to reference set of 16S sequences [REQUIRED]
	-o, `-`-output_dir
		The output directory [REQUIRED]
	-p, `-`-parameter_fp
		Path to the parameter file [REQUIRED]
	
	**[OPTIONAL]**
		
	`-`-experiment_attribute_fp
		Three-column, tab-delimited file of experiment attributes [default: None]
	`-`-experiment_link_fp
		Three-column, tab-delimited file of experiment links [default: None]
	-w, `-`-print_only
		Print the commands but don't call them [default: False]
	`-`-remove_unassigned
		Comma-separated list of run prefixes for which to remove unassigned sequences [default: None]


**Output:**

Produces all the files necessary to make an SRA submission, including an archive of per-library SFF files, and XML files for the experiment, runs, and submission.


**Example:**

::

	process_sra_submission.py -s sff_dir -e experiment.txt -r greengenes_unaligned.fasta -u submission.txt -p sra_parameters.txt -o out_dir/


