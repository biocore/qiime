.. _filtering_contamination_otus:

===============================================================
Filtering contaminant or category specific OTUs from OTU tables
===============================================================

This tutorial explains how to use several QIIME scripts to filter all OTUs that belong to a particular category of samples. This would be used, for example, if you ran blank control samples and want to remove any OTUs observed in these samples as likely contamination.

To accomplish this task several scripts are used to first generate an OTU table of target OTUs to remove, then filter those OTUs from the original OTU table, and finally to perform a clean-up step to remove the control samples that would now have zero sequences associated with them.

The OTU table and mapping file (generated from the QIIME tutorial data set) are available `here <https://s3.amazonaws.com/s3-qiime_tutorial_files/filtering_otus_tutorial_data.zip>`_.

Once these files are downloaded and extracted, open a terminal and change to the directory of the extracted files to begin processing.

Filtering out samples according to run
======================================

In this case, we are going to assume that multiple runs are present in an OTU table, and these are indicated in the ``Run_Number`` column in our mapping file.  As our example is removal of all OTUs from samples that should be blank control samples, we can assume that contamination will be limited to a single run, so we therefore want to begin by splitting the OTU table by the ``Run_Number`` field.  If there are not multiple runs to separate, this step can be skipped.

Splitting the OTU table by run
------------------------------

The runs are identified in the Run_Number column of the example mapping file.

.. note::

	* #SampleID	BarcodeSequence	LinkerPrimerSequence	Run_Number	Sample_Type	Description
	* #Modified tutorial mapping file to show procedure to filter OTUs from particular sample categories (ie. Control samples whose sequences are likely contamination across all samples)					
	* PC.354	AGCACGAGCCTA	YATGCTGCCTCCCGTAGGAGT	1	Control_Blank	Control_mouse_I.D._354
	* PC.355	AACTCGTCGATG	YATGCTGCCTCCCGTAGGAGT	1	Control_Blank	Control_mouse_I.D._355
	* PC.356	ACAGACCACTCA	YATGCTGCCTCCCGTAGGAGT	2	Control_Blank	Control_mouse_I.D._356
	* PC.481	ACCAGCGACTAG	YATGCTGCCTCCCGTAGGAGT	2	Control_Blank	Control_mouse_I.D._481
	* PC.593	AGCAGCACTTGT	YATGCTGCCTCCCGTAGGAGT	2	Control_Blank	Control_mouse_I.D._593
	* PC.607	AACTGTGCGTAC	YATGCTGCCTCCCGTAGGAGT	1	Test_Sample	Fasting_mouse_I.D._607
	* PC.634	ACAGAGTCGGCT	YATGCTGCCTCCCGTAGGAGT	1	Test_Sample	Fasting_mouse_I.D._634
	* PC.635	ACCGCAGAGTCA	YATGCTGCCTCCCGTAGGAGT	2	Test_Sample	Fasting_mouse_I.D._635
	* PC.636	ACGGTGAGTGTC	YATGCTGCCTCCCGTAGGAGT	2	Test_Sample	Fasting_mouse_I.D._636
	
To create per-run OTU tables containing, use the following command: ::

	split_otu_table.py -i otu_table.biom -m map.txt -f Run_Number -o split_otu_tables/

One can observe the initial sequences/sample in the run 1 OTU table: ::

	print_biom_table_summary.py -i split_otu_tables/otu_table_1.biom 

.. note::

	* Num samples: 4
	* Num otus: 419
	* Num observations (sequences): 595.0
	* 
	* Seqs/sample summary:
	*  Min: 147.0
	*  Max: 150.0
	*  Median: 149.0
	*  Mean: 148.75
	*  Std. dev.: 1.08972473589
	*  Median Absolute Deviation: 0.5
	*  Default even sampling depth in
 	*  core_qiime_analyses.py (just a suggestion): 149.0
	* 
	* Seqs/sample detail:
	*  PC.355: 147.0
	*  PC.354: 149.0
	*  PC.607: 149.0
	*  PC.634: 150.0

Filtering OTUs observed in the control blanks from the experimental samples
---------------------------------------------------------------------------------------

We will only examine run 1 in this example, but if you're working with multiple runs of data you would apply this step for each run. Note that if you don't have multiple runs, you would continue with ``otu_table.biom`` at this stage, rather than ``split_otu_tables/otu_table_1.biom``.

Create an OTU table with the just the blank control samples: ::

	filter_samples_from_otu_table.py -i split_otu_tables/otu_table_1.biom -o otu_table_run1_blank_samples.biom -m map.txt -s "Sample_Type:Control_Blank"
	
Filter out OTU ids that have zero counts, as we only want the OTUs with positive counts from the Control_Blank samples::

	filter_otus_from_otu_table.py -i otu_table_run1_blank_samples.biom -o filtered_otu_table_blank_samples.biom -n 1
	
Then create a tab separated version of this OTU table::

	biom convert -b -i filtered_otu_table_blank_samples.biom -o otus_to_remove.txt
	
Filter out OTU ids from the run 1 OTU table that were determined to be present in the Control_Blank samples::

	filter_otus_from_otu_table.py -i split_otu_tables/otu_table_1.biom -o otu_table_1_minus_contaminants.biom -e otus_to_remove.txt
	
The otu_table_1_minus_contaminants.biom file now has two samples with zero sequences associated with it.  These can be removed to get a final OTU table: ::

	filter_samples_from_otu_table.py -i otu_table_1_minus_contaminants.biom -o final_otu_table_1_minus_contaminants.biom -n 1
	
The final OTU table sequences/sample summary can be displayed now, sans OTUs from the Control_Blank samples: ::

	print_biom_table_summary.py -i final_otu_table_1_minus_contaminants.biom

.. note::

	* Num samples: 2
	* Num otus: 313
	* Num observations (sequences): 209.0
	* 
	* Seqs/sample summary:
	*  Min: 99.0
	*  Max: 110.0
	*  Median: 104.5
	*  Mean: 104.5
	*  Std. dev.: 5.5
	*  Median Absolute Deviation: 5.5
	*  Default even sampling depth in
	*   core_qiime_analyses.py (just a suggestion): 99.0
	* 
	* Seqs/sample detail:
	*  PC.607: 99.0
	*  PC.634: 110.0


If you apply this process to multiple runs, and then want to reassemble the final OTU tables into a single OTU table, you can use the ``merge_otu_tables.py`` command.
