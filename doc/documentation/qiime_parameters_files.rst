.. _qiime_parameter_files:

======================
QIIME parameters files
======================

QIIME parameters files are used with the QIIME workflow scripts, which plug together multiple QIIME commands. 

Motivation
==========

The QIIME workflow scripts plug together two or more QIIME commands to facilitate running common processes. For example, `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_ combines OTU picking, representative sequence picking, taxonomy assignment, sequence alignment, alignment filtering, and tree building in a single QIIME command. If you were to run each of these steps yourself, you would run about seven different commands (depending a little on what options you use). Because each of these individual commands may have many options, it's not feasible to make of these options accessible through the `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_  interface (to see why, run ``pick_otus.py -h`` to get list of options for just that single script). However, we still want users to have control over the different options available in the individual scripts. To achieve that, you can pass a QIIME parameters file to the workflow scripts to modify options to the component scripts of the workflow.

Format
======

The parameters file is a text file with one parameter setting per line. Blank lines or lines beginning with a ``#`` are ignored. A parameter setting is defined as ``script_name:parameter_name``, followed by a tab, and then the value. For example::
	
	align_seqs:alignment_method	pynast

This indicates that the ``--otu_picking_method`` will be set to ``uclust`` when calling ``pick_otus.py``. To get information on what a parameter in the ``qiime_parameters.txt`` file is, you should call the script name followed by ``-h`` to access the usage information for that script. In the above example, you could call::
	
	pick_otus.py -h

Flag options (i.e., those that don't take a value, like the ``--enable_rev_strand_match`` option to ``pick_otus.py``) are specified by passing ``True`` or ``False`` after the tab. For example::
	
	pick_otus:enable_rev_strand_match	True
	
If a parameter is not followed by an option, it will not be passed to the script resulting in the default value being used instead. For example::
	
	pick_otus:similarity

This results in no ``--similarity`` parameter being passed to ``pick_otus.py``. Alternatively, you can delete this line from your parameters file.

Note that when specifying parameters, you must specify their `long form` name, not their single-letter abbreviation. For example, if a parameter is defined as follows::

	-s SIMILARITY, --similarity=SIMILARITY
	                      Sequence similarity threshold (for blast, cdhit,
	                      uclust, uclust_ref, or usearch) [default: 0.97]

You would pass this as ``pick_otus:similarity``, NOT as ``pick_otus:s``. The latter will not work!


You can find information on the QIIME workflow scripts at:

	* `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_
	* `alpha_rarefaction.py <../scripts/alpha_rarefaction.html>`_
	* `beta_diversity_through_plots.py <../scripts/beta_diversity_through_plots.html>`_
	* `jackknifed_beta_diversity.py <../scripts/jackknifed_beta_diversity.html>`_




