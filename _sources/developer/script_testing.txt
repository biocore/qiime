
Automated script testing
^^^^^^^^^^^^^^^^^^^^^^^^

QIIME is well-tested: we make automated software testing a high priority, as you can see by reviewing the ``Qiime/tests/`` directory `here <https://github.com/biocore/qiime/tree/master/tests>`_. After installing QIIME, users can run the automated test suite by changing to the ``Qiime/tests`` directory and running::

	./all_tests.py

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

To detect these types of errors, we've developed a script interface testing framework that makes use of the ``usage_examples`` in the ``script_info`` object that is implemented by all scripts. This framework relies on a set of example files that can be used with QIIME, as well as a script for running the script usage tests. This is all contained in the ``qiime_test_data`` directory in the QIIME repository. 

The following sections first illustrate how to apply this testing framework locally, and then how to develop QIIME scripts so they can be used with this framework.

Running script usage tests
============================================

After obtaining the ``qiime`` repository, you can ``cd`` to ``qiime/tests`` directory. You'll see a script called ``all_tests.py``. You can run this from the ``qiime_test_data`` directory by calling::

	./all_tests.py --suppress_unit_tests

This will run all of the tests which are currently defined in verbose mode. You can run specific tests by passing the names of those tests via the ``--script_usage_tests`` parameter. For example, to run only the tests for ``add_qiime_labels.py`` and ``make_otu_table.py`` you can run the following::

	./all_tests.py --suppress_unit_tests --script_usage_tests add_qiime_labels,make_otu_table

These tests will print output to the screen.

The recommended way of testing QIIME is to run QIIME's unit tests and script usage tests. You can do this simply by running::

	./all_tests.py

How the script usage tests work
===============================
You'll see many sub-directories in ``qiime_test_data`` with names corresponding to the names of QIIME scripts. Each of these directories contains example input files for the corresponding QIIME script. For example, the ``make_otu_table`` directory contains the following test data for ``make_otu_table.py``::

	ls make_otu_table/
	chimeric_seqs.txt			otu_map.txt				pynast_failures.fna		tax_assignments.txt

If you call ``make_otu_table.py -h``, you'll see the following usage examples::

	Example usage:
	Print help message and exit
	 make_otu_table.py -h

	Make OTU table: Make an OTU table from an OTU map (i.e., result from pick_otus.py) and a taxonomy assignment file (i.e., result from assign_taxonomy.py). Write the output file to otu_table.biom.
	 make_otu_table.py -i otu_map.txt -t tax_assignments.txt -o otu_table.biom

	Make OTU table, excluding OTU ids listed in a fasta file: Make an OTU table, excluding the sequences listed in pynast_failures.fna. Note that the file pass as -e must end with either '.fasta' or '.fna'.
	 make_otu_table.py -i otu_map.txt -o otu_table_no_pynast_failures.biom -e pynast_failures.fna

	Make OTU table, excluding a list of OTU ids: Make an OTU table, excluding the sequences listed in chimeric_seqs.txt
	 make_otu_table.py -i otu_map.txt -o otu_table_non_chimeric.biom -e chimeric_seqs.txt

What you'll notice is that the usage example input files correspond to the files in ``qiime/qiime_test_data/make_otu_table``. The script interface testing works by copying all of the files in ``make_otu_table`` to a temporary directory, changing into that directory, running each of the usage examples, and confirming that the script exited successfully (i.e., with an exit status of ``0``).

 .. warning:: Currently the script usage tests only test whether a script exits successfully: they do not check whether the results correspond to expected output. The reasoning is that that would duplicate the functionality of the unit tests (which isn't a bad thing, except that implementing this would be a lot of work). These are tests that the interfaces themselves are working.

If you don't see a directory corresponding to a script name in the ``qiime_test_data`` directory, that means that a script interface test has not been defined for the given script. We're currently working on extending this framework to cover all QIIME scripts.

Adding script interface testing for new scripts
===============================================

Adding new script interface tests is easy. All you need to do is create a new test directory in ``qiime_test_data``, where the name of the directory corresponds to the script's name. For example, if you're adding tests for ``my_script.py``, you'd add a directory called ``my_script``. In that directory you would create example input files for all of the script usage examples that are defined in your script. Make several usage examples that make use of different paths through your script. 

Full paths
----------
We recommend specifying full paths for many of QIIME scripts, and importantly for workflow and parallel scripts. To do this in your usage example, replace the full path with $PWD. For example (from ``pick_de_novo_otus.py``)::

	Simple example: The following command will start an analysis on seqs.fna (-i), which is a
	post-split_libraries fasta file. The sequence identifiers in this file should be of the form
	<sample_id>_<unique_seq_id>. The following steps, corresponding to the preliminary data 
	preparation, are applied: Pick de novo OTUs at 97%; pick a representative sequence for each 
	OTU (the OTU centroid sequence); align the representative set with PyNAST; assign taxonomy 
	with the uclust consensus taxonomy assigner; filter the alignment prior to tree building - remove positions which 
	are all gaps, and specified as 0 in the lanemask; build a phylogenetic tree with FastTree; 
	build an OTU table. All output files will be written to the directory specified by -o, and 
	subdirectories as appropriate. ALWAYS SPECIFY ABSOLUTE FILE PATHS (absolute path represented 
	here as $PWD, but will generally look something like /home/ubuntu/my_analysis/).
	 pick_de_novo_otus.py -i $PWD/seqs.fna -o $PWD/otus/
