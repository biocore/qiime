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
	
	pick_otus:otu_picking_method	uclust

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

Examples
========

Below are some commonly used parameter files.

Run OTU picking with uclust's default parameters (rather than QIIME's default parameters, which are more conservative but *much* slower). Also, enable reverse strand matching so reads in opposite orientations match one another (i.e., uclust's ``--rev`` parameter). This information can be saved to a text file and passed to `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_, `pick_reference_otus_through_otu_table.py <../scripts/pick_reference_otus_through_otu_table.html>`_, `pick_subsampled_reference_otus_through_otu_table.py <../scripts/pick_subsampled_reference_otus_through_otu_table.html>`_, or `core_qiime_analyses.py <../scripts/core_qiime_analyses.html>`_.
::
	
	pick_otus:enable_rev_strand_match True
	pick_otus:max_accepts 1
	pick_otus:max_rejects 8
	pick_otus:stepwords 8
	pick_otus:word_length 8

Run beta diversity with Unweighted UniFrac, Weighted UniFrac, Bray-Curtis, and Euclidean distance metrics. This can be used with `beta_diversity_through_plots.py <../scripts/beta_diversity_through_plots.html>`_ and `core_qiime_analyses.py <../scripts/core_qiime_analyses.html>`_.
::
	
	beta_diversity:metrics	bray_curtis,euclidean,unweighted_unifrac,weighted_unifrac

Workflow scripts
=================

You can find information on the QIIME workflow scripts at:

	* `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_
	* `pick_reference_otus_through_otu_table.py <../scripts/pick_reference_otus_through_otu_table.html>`_
	* `pick_subsampled_reference_otus_through_otu_table.py <../scripts/pick_subsampled_reference_otus_through_otu_table.html>`_
	* `alpha_rarefaction.py <../scripts/alpha_rarefaction.html>`_
	* `beta_diversity_through_plots.py <../scripts/beta_diversity_through_plots.html>`_
	* `summarize_taxa_through_plots.py <../scripts/summarize_taxa_through_plots.html>`_
	* `jackknifed_beta_diversity.py <../scripts/jackknifed_beta_diversity.html>`_
	* `core_qiime_analyses.py <../scripts/core_qiime_analyses.html>`_

Why is there no longer an example parameters file included with QIIME?
======================================================================

In the past we provided an example QIIME parameters file, which contained all of the parameters to all of the subcommands for all of the workflow scripts. This was really messy and hard to maintain: if we updated a default parameter in a QIIME script, we had to also update it in the parameters file, or using the example parameters file would default to an old parameter setting. Initially, the parameters files were also required options to the workflow scripts, but they are now optional, so we now recommend using the parameters files in a different way that we have in the past.

Previously we recommended starting with the example QIIME parameters file and modifying it change the parameters used during your workflow. We now recommend starting with no parameters file (the default) and creating one with only the parameters you want to modify, if you want to change the default behavior. We therefore no longer provide an example parameters file. To find the parameters that are available to modify for a given workflow, you can call the individual component scripts (i.e., the scripts that are called by the workflow scripts) with ``-h`` to see what parameters are available to be modified. We hope to make this even simpler in the future as we work on a general refactoring of the QIIME workflow functionality. 




