.. _pick_otus:

.. index:: pick_otus.py

*pick_otus.py* -- OTU picking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The OTU picking step assigns similar sequences to operational taxonomic units, or OTUs, by clustering sequences based on a user-defined similarity threshold. Sequences which are similar at or above the threshold level are taken to represent the presence of a taxonomic unit (e.g., a genus, when the similarity threshold is set at 0.94) in the sequence collection.

Currently, the following clustering methods have been implemented in QIIME:

1. cd-hit (Li & Godzik, 2006; Li, Jaroszewski, & Godzik, 2001), which applies a "longest-sequence-first list removal algorithm" to cluster sequences.  

2. blast (Altschul, Gish, Miller, Myers, & Lipman, 1990), which compares and clusters each sequence against a reference database of sequences.

3. Mothur (Schloss et al., 2009), which requires an input file of aligned sequences.  The input file of aligned sequences may be generated from an input file like the one described below by running `align_seqs.py <./align_seqs.html>`_.  For the Mothur method, the clustering algorithm may be specified as nearest-neighbor, furthest-neighbor, or average-neighbor.  The default algorithm is furthest-neighbor.

4. prefix/suffix [Qiime team, unpublished], which will collapse sequences which are identical in their first and/or last bases (i.e., their prefix and/or suffix). The prefix and suffix lengths are provided by the user and default to 50 each.

5. Trie [Qiime team, unpublished], which collapsing identical sequences and sequences which are subsequences of other sequences.

6. uclust (Edgar, RC 2010), creates "seeds" of sequences which generate clusters based on percent identity.

7. usearch (Edgar, RC 2010), creates "seeds" of sequences which generate clusters based on percent identity, filters low abundance clusters, performs de novo and reference based chimera detection.

The primary inputs for `pick_otus.py <./pick_otus.html>`_ are:

1. A FASTA file containing sequences to be clustered

2. An OTU threshold (default is 0.97, roughly corresponding to species-level OTUs);

3. The method to be applied for clustering sequences into OTUs.

`pick_otus.py <./pick_otus.html>`_ takes a standard fasta file as input.




**Usage:** :file:`pick_otus.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_seqs_filepath
		Path to input sequences file
	
	**[OPTIONAL]**
		
	-m, `-`-otu_picking_method
		Method for picking OTUs.  Valid choices are: usearch, usearch_ref, prefix_suffix, mothur, trie, blast, uclust_ref, cdhit, uclust. The mothur method requires an input file of aligned sequences.  usearch will enable the usearch quality filtering pipeline. [default: uclust]
	-c, `-`-clustering_algorithm
		Clustering algorithm for mothur otu picking method.  Valid choices are: furthest, nearest, average. [default: furthest]
	-M, `-`-max_cdhit_memory
		Maximum available memory to cd-hit-est (via the program's -M option) for cdhit OTU picking method (units of Mbyte) [default: 400]
	-o, `-`-output_dir
		Path to store result file [default: ./<OTU_METHOD>_picked_otus/]
	-r, `-`-refseqs_fp
		Path to reference sequences to search against when using -m blast, -m uclust_ref, or -m usearch_ref [default: None]
	-b, `-`-blast_db
		Pre-existing database to blast against when using -m blast [default: None]
	`-`-min_aligned_percent
		Minimum percent of query sequence that can be aligned to consider a hit (BLAST OTU picker only) [default: 0.5]
	-s, `-`-similarity
		Sequence similarity threshold (for blast, cdhit, uclust, uclust_ref, or usearch) [default: 0.97]
	-e, `-`-max_e_value
		Max E-value when clustering with BLAST [default: 1e-10]
	-q, `-`-trie_reverse_seqs
		Reverse seqs before picking OTUs with the Trie OTU picker for suffix (rather than prefix) collapsing [default: False]
	-n, `-`-prefix_prefilter_length
		Prefilter data so seqs with identical first prefix_prefilter_length are automatically grouped into a single OTU.  This is useful for large sequence collections where OTU picking doesn't scale well [default: None; 100 is a good value]
	-t, `-`-trie_prefilter
		Prefilter data so seqs which are identical prefixes of a longer seq are automatically grouped into a single OTU; useful for large sequence collections where OTU picking doesn't scale well [default: False]
	-p, `-`-prefix_length
		Prefix length when using the prefix_suffix otu picker; WARNING: CURRENTLY DIFFERENT FROM prefix_prefilter_length (-n)! [default: 50]
	-u, `-`-suffix_length
		Suffix length when using the prefix_suffix otu picker [default: 50]
	-z, `-`-enable_rev_strand_match
		Enable reverse strand matching for uclust otu picking, will double the amount of memory used. [default: False]
	-D, `-`-suppress_presort_by_abundance_uclust
		Suppress presorting of sequences by abundance when picking OTUs with uclust or uclust_ref [default: False]
	-A, `-`-optimal_uclust
		Pass the --optimal flag to uclust for uclust otu picking. [default: False]
	-E, `-`-exact_uclust
		Pass the --exact flag to uclust for uclust otu picking. [default: False]
	-B, `-`-user_sort
		Pass the --user_sort flag to uclust for uclust otu picking. [default: False]
	-C, `-`-suppress_new_clusters
		Suppress creation of new clusters using seqs that don't match reference when using -m uclust_ref or -m usearch_ref [default: False]
	`-`-max_accepts
		Max_accepts value to uclust and uclust_ref [default: 20]
	`-`-max_rejects
		Max_rejects value to uclust and uclust_ref [default: 500]
	`-`-stepwords
		Stepwords value to uclust and uclust_ref [default: 20]
	`-`-word_length
		W value to usearch, uclust, and uclust_ref.  Set to 64 for usearch. [default: 12]
	`-`-uclust_otu_id_prefix
		OTU identifier prefix (string) for the de novo uclust OTU picker and for new clusters when uclust_ref is used without -C [default: None, OTU ids are ascending integers]
	`-`-suppress_uclust_stable_sort
		Don't pass --stable-sort to uclust [default: False]
	`-`-suppress_uclust_prefilter_exact_match
		Don't collapse exact matches before calling uclust [default: False]
	-d, `-`-save_uc_files
		Enable preservation of intermediate uclust (.uc) files that are used to generate clusters via uclust.  Also enables preservation of all intermediate files created by usearch (usearch_qf). [default: True]
	-j, `-`-percent_id_err
		Percent identity threshold for cluster error detection with usearch_qf. [default: 0.97]
	-g, `-`-minsize
		Minimum cluster size for size filtering with usearch_qf. [default: 4]
	-a, `-`-abundance_skew
		Abundance skew setting for de novo chimera detection with usearch_qf. [default: 2.0]
	-f, `-`-db_filepath
		Reference database of fasta sequences for reference based chimera detection with usearch_qf. [default: None]
	`-`-perc_id_blast
		Percent ID for mapping OTUs created by usearch_qf back to original sequence IDs [default: 0.97]
	`-`-de_novo_chimera_detection
		Deprecated:  de novo chimera detection performed by default, pass --suppress_de_novo_chimera_detection to disable. [default: None]
	-k, `-`-suppress_de_novo_chimera_detection
		Suppress de novo chimera detection in usearch_qf. [default: False]
	`-`-reference_chimera_detection
		Deprecated:  Reference based chimera detection performed by default, pass --supress_reference_chimera_detection to disable [default: None]
	-x, `-`-suppress_reference_chimera_detection
		Suppress reference based chimera detection in usearch_qf. [default: False]
	`-`-cluster_size_filtering
		Deprecated, cluster size filtering enabled by default, pass --disable_cluster_size_filtering to disable.  [default: None]
	-l, `-`-suppress_cluster_size_filtering
		Suppress cluster size filtering in usearch_qf.  [default: False]
	`-`-remove_usearch_logs
		Disable creation of logs when usearch is called.  Up to nine logs are created, depending on filtering steps enabled.  [default: False]
	`-`-derep_fullseq
		Dereplication of full sequences, instead of subsequences. Faster than the default --derep_subseqs in usearch. [default: False]
	-F, `-`-non_chimeras_retention
		Selects subsets of sequences detected as non-chimeras to retain after de novo and refernece based chimera detection.  Options are intersection or union.  union will retain sequences that are flagged as non-chimeric from either filter, while intersection will retain only those sequences that are flagged as non-chimeras from both detection methods. [default: union]
	`-`-minlen
		Minimum length of sequence allowed for usearch. [default: 64]


**Output:**

The output consists of two files (i.e. seqs_otus.txt and seqs_otus.log). The .txt file is composed of tab-delimited lines, where the first field on each line corresponds to an (arbitrary) cluster identifier, and the remaining fields correspond to sequence identifiers assigned to that cluster. Sequence identifiers correspond to those provided in the input FASTA file.  Usearch (i.e. usearch quality filter) can additionally have log files for each intermediate call to usearch.

Example lines from the resulting .txt file:

=   ====    ====    ====
0   seq1    seq5        
1   seq2                
2   seq3                
3   seq4    seq6    seq7
=   ====    ====    ====

This result implies that four clusters were created based on 7 input sequences. The first cluster (cluster id 0) contains two sequences, sequence ids seq1 and seq5; the second cluster (cluster id 1) contains one sequence, sequence id seq2; the third cluster (cluster id 2) contains one sequence, sequence id seq3, and the final cluster (cluster id 3) contains three sequences, sequence ids seq4, seq6, and seq7.

The resulting .log file contains a list of parameters passed to the `pick_otus.py <./pick_otus.html>`_ script along with the output location of the resulting .txt file.


**Example (uclust method, default):**

Using the seqs.fna file generated from `split_libraries.py <./split_libraries.html>`_ and outputting the results to the directory "picked_otus_default/", while using default parameters (0.97 sequence similarity, no reverse strand matching):

::

	pick_otus.py -i seqs.fna -o picked_otus_default

To change the percent identity to a lower value, such as 90%, and also enable reverse strand matching, the command would be the following:

::

	pick_otus.py -i seqs.fna -o picked_otus_90_percent_rev/ -s 0.90 -z

**Uclust Reference-based OTU picking example:**

uclust_ref can be passed via -m to pick OTUs against a reference set where sequences within the similarity threshold to a reference sequence will cluster to an OTU defined by that reference sequence, and sequences outside of the similarity threshold to a reference sequence will form new clusters. OTU identifiers will be set to reference sequence identifiers when sequences cluster to reference sequences, and 'qiime_otu_<integer>' for new OTUs. Creation of new clusters can be suppressed by passing -C, in which case sequences outside of the similarity threshold to any reference sequence will be listed as failures in the log file, and not included in any OTU.

::

	pick_otus.py -i seqs.fna -r refseqs.fasta -m uclust_ref --uclust_otu_id_prefix qiime_otu_

**Example (cdhit method):**

Using the seqs.fna file generated from `split_libraries.py <./split_libraries.html>`_ and outputting the results to the directory "cdhit_picked_otus/", while using default parameters (0.97 sequence similarity, no prefix filtering):

::

	pick_otus.py -i seqs.fna -m cdhit -o cdhit_picked_otus/

Currently the cd-hit OTU picker allows for users to perform a pre-filtering step, so that highly similar sequences are clustered prior to OTU picking. This works by collapsing sequences which begin with an identical n-base prefix, where n is specified by the -n parameter. A commonly used value here is 100 (e.g., -n 100). So, if using this filter with -n 100, all sequences which are identical in their first 100 bases will be clustered together, and only one representative sequence from each cluster will be passed to cd-hit. This is used to greatly decrease the run-time of cd-hit-based OTU picking when working with very large sequence collections, as shown by the following command:

::

	pick_otus.py -i seqs.fna -m cdhit -o cdhit_picked_otus_filter/ -n 100

Alternatively, if the user would like to collapse identical sequences, or those which are subsequences of other sequences prior to OTU picking, they can use the trie prefiltering ("-t") option as shown by the following command.

Note: It is highly recommended to use one of the prefiltering methods when analyzing large datasets (>100,000 seqs) to reduce run-time.

::

	pick_otus.py -i seqs.fna -m cdhit -o cdhit_picked_otus_trie_prefilter/ -t

**BLAST OTU-Picking Example:**

OTUs can be picked against a reference database using the BLAST OTU picker. This is useful, for example, when different regions of the SSU RNA have sequenced and a sequence similarity based approach like cd-hit therefore wouldn't work. When using the BLAST OTU picking method, the user must supply either a reference set of sequences or a reference database to compare against. The OTU identifiers resulting from this step will be the sequence identifiers in the reference database. This allows for use of a pre-existing tree in downstream analyses, which again is useful in cases where different regions of the 16s gene have been sequenced.

The following command can be used to blast against a reference sequence set, using the default E-value and sequence similarity (0.97) parameters:

::

	pick_otus.py -i seqs.fna -o blast_picked_otus/ -m blast -r refseqs.fasta

If you already have a pre-built BLAST database, you can pass the database prefix as shown by the following command:

::

	pick_otus.py -i seqs.fna -o blast_picked_otus_prebuilt_db/ -m blast -b refseqs.fasta

If the user would like to change the sequence similarity ("-s") and/or the E-value ("-e") for the blast method, they can use the following command:

::

	pick_otus.py -i seqs.fna -o blast_picked_otus_90_percent/ -m blast -r refseqs.fasta -s 0.90 -e 1e-30

**Prefix-suffix OTU Picking Example:**

OTUs can be picked by collapsing sequences which begin and/or end with identical bases (i.e., identical prefixes or suffixes).  This OTU picker is currently likely to be of limited use on its own, but will be very useful in collapsing very similar sequences in a chained OTU picking strategy that is currently in development. For example, the user will be able to pick OTUs with this method, followed by representative set picking, and then re-pick OTUs on their representative set. This will allow for highly similar sequences to be collapsed, followed by running a slower OTU picker. This ability to chain OTU pickers is not yet supported in QIIME. The following command illustrates how to pick OTUs by collapsing sequences which are identical in their first 50 and last 25 bases:

::

	pick_otus.py -i seqs.fna -o prefix_suffix_picked_otus/ -m prefix_suffix -p 50 -u 25

**Mothur OTU Picking Example:**

The Mothur program (http://www.mothur.org/) provides three clustering algorithms for OTU formation: furthest-neighbor (complete linkage), average-neighbor (group average), and nearest-neighbor (single linkage). Details on the algorithms may be found on the Mothur website and publications (Schloss et al., 2009). However, the running times of Mothur's clustering algorithms scale with the number of sequences squared, so the program may not be feasible for large data sets.

The following command may be used to create OTUs based on a furthest-neighbor algorithm (the default setting) using aligned sequences as input:

::

	pick_otus.py -i seqs.aligned.fna -o mothur_picked_otus/ -m mothur

If you prefer to use a nearest-neighbor algorithm instead, you may specify this with the '-c' flag:

::

	pick_otus.py -i seqs.aligned.fna -o mothur_picked_otus_nn/ -m mothur -c nearest

The sequence similarity parameter may also be specified. For example, the following command may be used to create OTUs at the level of 90% similarity:

::

	pick_otus.py -i seqs.aligned.fna -o mothur_picked_otus_90_percent/ -m mothur -s 0.90

**Usearch_qf ('usearch quality filter'):**

Usearch (http://www.drive5.com/usearch/) provides clustering, chimera checking, and quality filtering. The following command specifies a minimum cluster size of 2 to be used during cluster size filtering:

::

	pick_otus.py -i seqs.fna -m usearch --word_length 64 --db_filepath refseqs.fasta -o usearch_qf_results/ --minsize 2

**Usearch (usearch_qf) example where reference-based chimera detection is disabled, and minimum cluster size filter is reduced from default (4) to 2:**

::

	pick_otus.py -i seqs.fna -m usearch --word_length 64 --suppress_reference_chimera_detection --minsize 2 -o usearch_qf_results_no_ref_chim_detection/


