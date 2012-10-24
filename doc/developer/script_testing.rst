
Automated script testing
^^^^^^^^^^^^^^^^^^^^^^^^

QIIME is well-tested: we make automated software testing a high priority, as you can see by reviewing the ``Qiime/tests/`` directory `here <https://github.com/qiime/qiime/tree/master/tests>`_. After installing QIIME, users can run the automated test suite by changing to the ``Qiime/tests`` directory and running::

	all_tests.py

These tests are additionally run on a nightly basis, and developers are alerted to test failures. 

One limitation of our testing strategy in QIIME however has been that it doesn't test the script interfaces, which are the primary point of interaction with users. So, for example, if script looks something like this:: 

	from qiime.filter import filter_samples_from_otu_table
	from qiime.filter import filter_otus_from_otu_table
	
	...
	
	if run_mode == 'filter_samples':
		filter_samples_from_otu_table(...)
	elif run_mode == 'filter_otus':
		filter_otus_from_otu_table(...)
	else:
		raise ValueError, "Unknown run_mode."

and a developer accidentally removes the ``from qiime.filter import filter_otus_from_otu_table`` line but only tests with ``run_mode == 'filter_samples'``, we may not catch the error before release of the software. These types of errors generally result in cryptic error messages, which frustrates users and increases our technical support load. As powerful as the underlying code may be, this simple, common mistake makes it useless to the end user.

To detect these types of errors, we've developed a script interface testing framework that makes use of the ``usage_examples`` in the ``script_info`` object that is implemented by all scripts. This framework relies on a set of example files that can be used with QIIME, as well as a script for running the script usage tests. This is all contained in the ``qiime_test_data`` repository. To obtain the test data, run the following command::

	git clone git://github.com/qiime-dev/qiime_test_data.git

The following sections first illustrate how to apply this testing framework locally, and then how to develop QIIME scripts so they can be used with this framework.

Running script usage tests
============================================

Interactive mode
----------------

After obtaining the ``qiime_test_data`` repository, you can ``cd`` to that directory. You'll see a script called script_usage_tests.py. You can run this from the ``qiime_test_data`` directory by calling::

	./script_usage_tests.py -i $PWD/ -l $HOME/all_qime_script_tests.log -v

This will run all of the tests which are currently defined in verbose mode. You can run specific tests by passing the names of those tests via the ``-t`` parameter. For example, to run only the tests for ``add_taxa.py`` and ``make_otu_table.py`` you can run the following::

	./script_usage_tests.py -i $PWD/ -l $HOME/qime_script_tests.log -t add_taxa,make_otu_table -v

These tests will print output to the screen because you've passed ``-v``. Regardless of whether you pass ``-v``, a summary of the test results will be written to the log path specified by ``-l``.

Automated mode
--------------

You can additionally set the script usage tests to run as part of ``all_tests.py`` by adding the full path to ``qiime_tests_data`` to your QIIME config file. For example::
	
	qiime_test_data_dir /home/ubuntu/qiime_test_data

``all_tests.py`` will detect this, and run all of the script usage tests that are currently defined.

How the script usage tests work
===============================
You'll see many sub-directories in ``qiime_test_data`` with names corresponding to the names of QIIME scripts. Each of these directories contains example input and output for the corresponding QIIME script. For example, the ``add_taxa`` directory contains the following test data for ``add_taxa.py``::

	ls add_taxa/
	otu_table_no_tax.biom	otu_table_w_score.biom	otu_table_w_tax_only.biom	tax.txt
	otu_table_w_alt_labeled_tax.biom	otu_table_w_tax.biom	score_only.txt	tax_only.txt

If you call ``add_taxa.py -h``, you'll see the following usage examples::

	Example: Given an input otu table with no metadata (otu_table_no_tax.biom) and a tab-separated 
	text file mapping OTU ids to taxonomic assignments and scores associated with those assignments
	(tax.txt), generate a new otu table that includes taxonomic assignments (otu_table_w_tax.biom).
	 add_taxa.py -i otu_table_no_tax.biom -o otu_table_w_tax.biom -t tax.txt

	Example: Given an input otu table with no metadata (otu_table_no_tax.biom) and a tab-separated
	text file mapping OTU ids to taxonomic assignments and scores associated with those assignments
	(tax.txt), generate a new otu table that includes taxonomic assignments (otu_table_w_tax.biom) 
	with alternate metadata identifiers.
	 add_taxa.py -i otu_table_no_tax.biom -o otu_table_w_alt_labeled_tax.biom -t tax.txt -l "Consensus Lineage,Score"

	Example: Given an input otu table with no metadata (otu_table_no_tax.biom) and a tab-separated 
	text file mapping OTU ids to some value, generate a new otu table that includes that metadata
	category labeled as "Score" (otu_table_w_score.biom).
	 add_taxa.py -i otu_table_no_tax.biom -o otu_table_w_score.biom -t score_only.txt -l "Score" --all_strings

What you'll notice is that the usage example input and output files correspond to the files in ``qiime_test_data/add_taxa``. The script interface testing works by copying all of the files in ``add_taxa`` to a temporary directory, changing into that directory, running each of the usage examples, and confirming that the script exited successfully (i.e., with an exit status of ``0``).

 .. warning:: Currently the script usage tests only test whether a script exits successfully: they do not check whether the results correspond to the example output. The reasoning is that that would duplicate the functionality of the unit tests (which isn't a bad thing, except that implementing this would be a lot of work). These are tests that the interfaces themselves are working.

If you don't see a directory corresponding to a script name in the ``qiime_test_data`` directory, that means that a script interface test has not been defined for the given script. We're currently working on extending this framework to cover all QIIME scripts.

Adding script interface testing for new scripts
===============================================

Adding new script interface tests is easy. All you do is create a new test directory in ``qiime_test_data``, where the name of the directory corresponds to the script's name. For example, if you're adding tests for ``my_script.py``, you'd add a directory called ``my_script``. In that directory you would create example input and output files for all of the script usage examples that are defined in your script. Make several usage examples that make use of different paths through your script. 

Full paths
----------
We recommend specifying full paths for many of QIIME scripts, and importantly for workflow and parallel scripts. To do this in your usage example, replace the full path with $PWD. For example (from ``pick_otus_through_otu_table.py``)::

	Simple example: The following command will start an analysis on seqs.fna (-i), which is a
	post-split_libraries fasta file. The sequence identifiers in this file should be of the form
	<sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data 
	preparation, are applied: Pick de novo OTUs at 97%; pick a representative sequence for each 
	OTU (the OTU centroid sequence); align the representative set with PyNAST; assign taxonomy 
	with RDP classifier; filter the alignment prior to tree building - remove positions which 
	are all gaps, and specified as 0 in the lanemask; build a phylogenetic tree with FastTree; 
	build an OTU table. All output files will be written to the directory specified by -o, and 
	subdirectories as appropriate. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented 
	here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).
	 pick_otus_through_otu_table.py -i $PWD/seqs.fna -o $PWD/otus/

Cleaning up output files
------------------------
Some scripts require that the user-specified output directory does not exist when the script runs, but we provide example output in the test directory. To automatically remove output directories prior to running the tests, add the ``script_usage_output_to_remove`` entry to your script info. For example, from ``pick_otus_through_otu_table.py``::

	script_info['script_usage_output_to_remove'] = ['$PWD/otus/']




