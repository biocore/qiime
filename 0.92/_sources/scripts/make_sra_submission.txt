.. _make_sra_submission:

.. index:: make_sra_submission.py

*make_sra_submission.py* -- Makes SRA submission files (xml-format)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script makes the submission xml files for SRA (study, experiment, etc.).  This script assumes that there is a simple tab-delimited text input (allowing for examples and comments).


**Usage:** :file:`make_sra_submission.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-a, `-`-input_sample_fp
		The tab-delimited text file with info about samples [default: None]
	`-`-template_sample_fp
		The template file for samples [default: sample_template.xml]
	-t, `-`-input_study_fp
		The tab-delimited text file with info about the study [default: None]
	`-`-template_study_fp
		The template file for the study [default: study_template.xml]
	-u, `-`-input_submission_fp
		The tab-delimited text file with info about the submission [default: None]
	`-`-template_submission_fp
		The template file for the submission [default: submission_template.xml]
	-e, `-`-input_experiment_fp
		The tab-delimited text file with info about the experiment [default: None]
	-s, `-`-sff_dir
		The directory containing the demultiplexed sff files: 1 dir per run [default: None]


**Output:**

This script produces 3 xml-formatted files.


**Example:**

Read the sample data from sample.txt, the study data from study.txt, and the submission data from submission.txt, which writes out the corresponding XML files.

::

	make_sra_submission.py -a sample.txt -t study.txt -u submission.txt

Produces the files study.xml, submission.xml, sample.xml (based on the filenames of the .txt files) using the default xml templates.


