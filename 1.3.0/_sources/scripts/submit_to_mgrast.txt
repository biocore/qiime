.. _submit_to_mgrast:

.. index:: submit_to_mgrast.py

*submit_to_mgrast.py* -- This script submits a FASTA file to MG-RAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

This script takes a split-library FASTA file and generates individual FASTA files for each sample, then submits each sample FASTA file to MG-RAST, given the user provides an MG-RAST web-services authorization key and Project ID.  To get a web-services authorization key, the user should have an account on MG-RAST.  Once logged in, the user can go to their Account Management page and under Preferences they should click 'here', where they will see a Web Services section where they can click on the 'generate new key' if they have not already been provided one.


**Usage:** :file:`submit_to_mgrast.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	-w, `-`-web_key_auth
		The web services authorization key from MG-RAST
	-p, `-`-project_id
		The title to be used for the project
	-o, `-`-output_dir
		Path to the output directory


**Output:**

The resulting directory will contain all of the sample-separated FASTA files, along with a log html file, which informs the user of the jobs started on MG-RAST


**Example:**

The user can submit a post-split-library FASTA file, which will be loaded and processed into MG-RAST under the users account ('-w') and project ('-p'), as follows:

::

	submit_to_mgrast.py -i split_lib_seqs.fna -w user_mgrast_auth_key -p qiime_test_dataset -o ./output_dir


