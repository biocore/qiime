.. _chaining_otu_pickers:

=======================
Multi-step OTU picking
=======================

This document describes how to perform chained or multi-step OTU picking, using the results from the QIIME Overview Tutorial. This is relevant, for example, when you have a very large collection of sequences, and you first want to apply a fast but rough OTU picker (e.g., PrefixSuffix) and then want to apply a slow but better OTU picker (e.g., cdhit). 

This example illustrates how to chain two OTU pickers, which is probably the most common usage. It is possible however to chain an arbitrary number of OTU pickers.

Step 1. Pick OTUs with the fast method using your full input sequence collection.
------------------------------------------------------------------------------------------
::
	
	pick_otus.py -m prefix_suffix -u 0 -i split_library_output/seqs.fna -o prefix_picked_otus
	
The resulting OTU map will look something like:
::
	
	0	PC.634_143	PC.634_196	PC.634_211
	1	PC.481_611


where 0 and 1 are OTU IDs, and PC.* are sequence IDs.
	
Step 2. Pick a representative set of sequences for the resulting OTU map.
------------------------------------------------------------------------------------------
::
	
	pick_rep_set.py -i prefix_picked_otus/seqs_otus.txt -f split_library_output/seqs.fna -o prefix_picked_otus/rep_set.fasta
	
The resulting fasta file will look something like:

::
	
	>0 PC.634_143
	TTGGGCCGTGTCTCAGTCCCAATGTGGCCGTTTACCCTCTCAGGCCGGCTACGCATCATCGCCTTGGTGGGCCGTTACCTCACCAACTAGCTAATGCGCCGCAGGTCCATCCATGTTCACGCCTTGATGGGCGCTTTAATATACTGAGCATGCGCTCTGTATACCTATCCGGTTTTAGCTACCGTTTCCAGCAGTTATCCCGGACACATGGGCAGGTT
	>1 PC.481_611
	TTGGTCCGTGTCTCAGTACCAATGTGGGGGTTAACCTCTCAGTCCCCTATGTATCGTCGCCTTGGTGAGCCGTTACCTCACCAACCAGCTAATACAACGCATGCCCATCCATAACCACCGGAGTTTTCAATCAAAAGGGATGCCCCTCTTGATGTTATGGGATATTAGTACCGATTTCTCAGTGTTATCCCCCTGTTATGGGTAGTTGCATACGCGTTACGCACCCGTGCGCCGGTCG


Step 3. Pick OTUs with the slow method using the representative set.
------------------------------------------------------------------------------------------
::
	
	pick_otus.py -m cdhit -i prefix_picked_otus/rep_set.fasta -o prefix_picked_otus/cdhit_picked_otus/
	
The resulting OTU map will look something like:
::
	
	0	39
	1	65	103 163

where the first column (containing 0 and 1) are the OTU IDs and the remaining values on each line are the OTU IDs from the first round of OTU picking.

Step 4. Next, you must merge your OTU maps, so your final OTU IDs map to input sequence identifiers. 
--------------------------------------------------------------------------------------------------------
You'll need to provide the OTU map filepaths generated above *in the order that they were generated!* The files paths (passed via -i) should be comma-separated, with no spaces.

::
	
	merge_otu_maps.py -i prefix_picked_otus/seqs_otus.txt,prefix_picked_otus/cdhit_picked_otus/rep_set_otus.txt -o otus.txt 

The resulting OTU map will look something like:
::
	
	133 PC.355_971
	132 PC.607_1118
	131 PC.354_823	PC.593_1312

where 133, 132 and 131 are the OTU IDs from step 3 (i.e., the second round of OTU picking), and PC.* are the original sequence identifiers.
	
Step 5. Finally, you can pick your representative set. 
------------------------------------------------------------------------------------------
At this stage, you'll pass the final OTU map (otus.txt) and the *original* sequence collection (split_library_output/seqs.fna).
::
	
	pick_rep_set.py -i otus.txt -f split_library_output/seqs.fna -o repr_set.fasta