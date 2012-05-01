.. _filtering_contamination_otus:


Filtering category specific OTUs from OTU tables
------------------------------------------------

This tutorial explains how to use several scripts of QIIME to filter all OTUs that belong to a particular category of samples (example:  blank control samples, where any sequences detected might be contamination affecting all samples).

To accomplish this task, several scripts are used, to first generate an OTU table of target OTUs to remove, followed by filtering the original OTU table of these OTUs and a clean-up step to remove the samples with zero sequences associated with them.

The sample OTU table and mapping file (generated from the QIIME tutorial data set) are available here:  <<<<< NEED URL >>>>> 

Once these files are downloaded and extracted, open a terminal and change to the directory of the extracted files to begin processing.

Filtering out samples according to run
======================================

In this case, we are going to assume that multiple runs are present in an OTU table.  As our example is removal of all OTUs from samples that should be blank control samples, we can assume that contamination will be limited to a single run.  If there are not multiple runs to separate, this step can be skipped.

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
	
To create an OTU table containing only samples from run 1, use the following command: ::

	filter_samples_from_otu_table.py -i otu_table.biom -o otu_table_run1.biom -m Sample_Mapping.txt -s 'Run_Number:1'

One can observe the initial sequences/sample in the run 1 OTU table: ::

	per_library_stats.py -i otu_table_run1.biom 

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



We will only examine run 1 in this example.  Next create an OTU table with the just the blank control samples: ::

	filter_samples_from_otu_table.py -i otu_table_run1.biom -o otu_table_run1_blank_samples.biom -m Sample_Mapping.txt -s "Sample_Type:Control_Blank"
	
Filter out OTU ids that have zero sequence counts, as we only want the OTUs with positive counts from the Control_Blank samples: ::

	filter_otus_from_otu_table.py -i otu_table_run1_blank_samples.biom -o filtered_otu_table_blank_samples.biom -n 1
	
Then create a tab separated version of this OTU table: ::

	convert_biom.py -b -i filtered_otu_table_blank_samples.biom -o otus_to_remove.txt
	
filter out OTU ids from the run 1 OTU table that are present in the Control_Blank samples: ::

	filter_otus_from_otu_table.py -i otu_table_run1.biom -o otus_removed_run1.biom -e otus_to_remove.txt
	
The otus_removed_run1.biom file now has two samples with zero sequences associated with it.  These can be removed to get a final OTU table: ::

	filter_samples_from_otu_table.py -i otus_removed_run1.biom -o final_filtered_otu_table.biom -n 1
	
The final OTU table sequences/sample summary can be displayed now, sans OTUs from the Control_Blank samples: ::

	per_library_stats.py -i final_filtered_otu_table.biom

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
