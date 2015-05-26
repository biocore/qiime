.. _chaining_otu_pickers:

=======================
Multi-step OTU picking
=======================

This document describes how to perform chained or multi-step OTU picking. This is relevant, for example, when you have a very large collection of sequences, and you first want to apply a fast but rough OTU picker (e.g., PrefixSuffix) and then want to apply a slow but better OTU picker (e.g., cdhit). 

This example illustrates how to chain two OTU pickers, which is probably the most common usage. It is possible however to chain an arbitrary number of OTU pickers.

Step 1. Pick OTUs with the fast method using your full input sequence collection.
------------------------------------------------------------------------------------------
::
	
	pick_otus.py -m prefix_suffix -u 0 -i seqs.fasta -o prefix_picked_otus
	
The resulting OTU map will look something like:
::
	
	0   seq1    seq2    seq5
	1   seq3    seq4    
	2   seq6    seq7    seq8

where 0, 1, and 2 are OTU ids, and seq* are sequence ids.
	
Step 2. Pick a representative set of sequences for the resulting OTU map.
------------------------------------------------------------------------------------------
::
	
	pick_rep_set.py -i prefix_picked_otus/seqs_otus.txt -f seqs.fasta -o prefix_picked_otus/repr_set.fasta
	
The resulting fasta file will look something like:

::
	
	>0
	ACCGGTAGAGAGATTAG
	>1
	ACGGTGCAGGGA
	>2
	CCTGAGGGGTTTGGAAAAAAG
	
Step 3. Pick OTUs with the slow method using the representative set.
------------------------------------------------------------------------------------------
::
	
	pick_otus.py -m cdhit -i prefix_picked_otus/repr_set.fasta -o prefix_picked_otus/cdhit_picked_otus/
	
The resulting OTU map will look something like:
::
	
	0 0   2
	1 1

where the first column (containing 0 and 1) are the OTU ids and the remaining values on each line are the OTU ids from the first round of OTU picking.

Step 4. Next, you must merge your OTU maps, so your final OTU ids map to input sequence identifiers. 
--------------------------------------------------------------------------------------------------------
You'll need to provide the OTU map filepaths generated above *in the order that they were generated!* The files paths (passed via -i) should be comma-separated, with no spaces.

::
	
	merge_otu_maps.py -i prefix_picked_otus/seqs_otus.txt,prefix_picked_otus/cdhit_picked_otus/rep_set_otus.txt -o otus.txt 

The resulting OTU map will look something like:
::
	
	0 seq1    seq2    seq5    seq6    seq7    seq8
	1 seq3    seq4

where 0 and 1 are the OTU ids from step 3 (i.e., the second round of OTU picking), and seq* are the original sequence identifiers.
	
Step 5. Finally, you can pick your representative set. 
------------------------------------------------------------------------------------------
At this stage, you'll pass the final OTU map (otus.txt) and the *original* sequence collection (seqs.fasta).
::
	
	pick_rep_set.py -i otus.txt -f seqs.fasta -o repr_set.fasta