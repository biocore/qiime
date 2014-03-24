.. _doc_qiime_parameters:

Defining custom parameter files for interacting with the workflow scripts
==========================================================================

The QIIME 'workflow' scripts are designed to help automate some of the steps of QIIME. These rely on the user defined their analysis in a parameters file. An example is provided as ``Qiime/qiime_parameters.txt``. Working with this file is described in this document.

qiime_parameters.txt
--------------------
This file is used to give workflow script users control over the parameters to the individual scripts without having an extremely complex and hard to maintain interface to the workflow scripts. Users should copy the example ``qiime_parameters.txt`` script to the directory where they are performing their analysis, and edit the values in this file accordingly. This copy will be referred to as the user's working ``qiime_parameters.txt`` file. The parameters are defined as ``script_name:parameter_name``, followed by a tab, and then the value. For example::
	
	align_seqs:alignment_method	pynast

This indicates that the ``--alignment_method`` will be set to ``pynast`` when calling ``align_seqs.py``. To get information on what a parameter in the ``qiime_parameters.txt`` file is, you should call the script name followed by ``-h`` to access the usage information for that script. In the above example, you could call::
	
	python Qiime/qiime/align_seqs.py -h

Boolean options are specified by passing ``True`` or ``False`` after the tab. For example::
	
	parallel:retain_temp_files	False
	
When a parameter is not followed by an option, that indicates that it will not be passed to the script resulting in the default value being used instead. For example::
	
	align_seqs:blast_db

This results in no ``--blast_db`` parameter being passed to ``align_seqs.py``. Alternatively, you can delete this line from your working ``qiime_parameters.txt`` file.

You can find information on the QIIME workflow scripts at:

	* `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_
	* `alpha_rarefaction.py <../scripts/alpha_rarefaction.html>`_
	* `beta_diversity_through_3d_plots.py <../scripts/beta_diversity_through_3d_plots.html>`_
	* `jackknifed_upgma.py <../scripts/jackknifed_upgma.html>`_
