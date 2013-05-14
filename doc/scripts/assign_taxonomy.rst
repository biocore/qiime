.. _assign_taxonomy:

.. index:: assign_taxonomy.py

*assign_taxonomy.py* -- Assign taxonomy to each sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Contains code for assigning taxonomy, using several techniques.

Given a set of sequences, assign_taxonomy.py attempts to assign the taxonomy of each sequence. Currently there are three methods implemented: assignment with BLAST, assignment with the RDP classifier, and assignment with the RTAX classifier. The output of this step is a mapping of input sequence identifiers (1st column of output file) to taxonomy (2nd column) and quality score (3rd column). The sequence identifier of the best BLAST hit is also included if the blast method is used (4th column).

Example reference data sets and id_to_taxonomy maps can be found in the Greengenes OTUs. To get the latest build of those click the "Most recent Greengenes OTUs" link on the top right of http://blog.qiime.org. After downloading and unzipping you can use the following following files as -r and -t. As of this writing the latest build was gg_otus_4feb2011, but that portion of path to these files will change with future builds. Modify these paths accordining when calling assign_taxonomy.py.

-r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta
-t gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train.txt (best for retraining the RDP classifier)
-t gg_otus_4feb2011/taxonomies/greengenes_tax.txt (best for BLAST taxonomy assignment)



**Usage:** :file:`assign_taxonomy.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-t, `-`-id_to_taxonomy_fp
		Path to tab-delimited file mapping sequences to assigned taxonomy. Each assigned taxonomy is provided as a semicolon-separated list. For assignment with rdp, each assigned taxonomy must be exactly 6 levels deep. [default: /Users/caporaso/data/gg_12_10_otus/taxonomy/97_otu_taxonomy.txt; REQUIRED when method is blast]
	-r, `-`-reference_seqs_fp
		Path to reference sequences.  For assignment with blast, these are used to generate a blast database. For assignment with rdp, they are used as training sequences for the classifier. [default: /Users/caporaso/data/gg_12_10_otus/rep_set/97_otus.fasta; REQUIRED if -b is not provided when method is blast]
	-p, `-`-training_data_properties_fp
		Path to ".properties" file in pre-compiled training data for the RDP Classifier.  This option is overridden by the -t and -r options. [default: None]
	`-`-read_1_seqs_fp
		Path to fasta file containing the first read from paired-end sequencing, prior to OTU clustering (used for RTAX only). [default: None]
	`-`-read_2_seqs_fp
		Path to fasta file containing a second read from paired-end sequencing, prior to OTU clustering (used for RTAX only). [default: None]
	`-`-single_ok
		When classifying paired ends, allow fallback to single-ended classification when the mate pair is lacking (used for RTAX only). [default: False]
	`-`-no_single_ok_generic
		When classifying paired ends, do not allow fallback to single-ended classification when the mate pair is overly generic (used for RTAX only). [default: False]
	`-`-read_id_regex
		Used to parse the result of OTU clustering, to get the read_1_id for each clusterID.  (used for RTAX only). [default: \S+\s+(\S+)]
	`-`-amplicon_id_regex
		Used to parse the result of split_libraries, to get the ampliconID for each read_1_id.  Two groups capture read_1_id and ampliconID, respectively.  (used for RTAX only). [default: (\S+)\s+(\S+?)\/]
	`-`-header_id_regex
		Used to choose the part of the header in the OTU clustering file that Rtax reports back as the ID.  The default uses the amplicon ID, not including /1 or /3, as the primary key for the query sequences. (used for RTAX only). [default: \S+\s+(\S+?)\/]
	-m, `-`-assignment_method
		Taxon assignment method, must be one of rdp, blast, rtax, mothur, tax2tree [default: rdp]
	-b, `-`-blast_db
		Database to blast against.  Must provide either --blast_db or --reference_seqs_db for assignment with blast [default: None]
	-c, `-`-confidence
		Minimum confidence to record an assignment, only used for rdp and mothur methods [default: 0.8]
	`-`-rdp_max_memory
		Maximum memory allocation, in MB, for Java virtual machine when using the rdp method.  Increase for large training sets [default: 1500]
	-e, `-`-e_value
		Maximum e-value to record an assignment, only used for blast method [default: 0.001]
	`-`-tree_fp
		The filepath to a prebuilt tree containing both the representative and reference sequences. Required for Tax2Tree assignment.
	-o, `-`-output_dir
		Path to store result file [default: <ASSIGNMENT_METHOD>_assigned_taxonomy]


**Output:**

The consensus taxonomy assignment implemented here is the most detailed lineage description shared by 90% or more of the sequences within the OTU (this level of agreement can be adjusted by the user). The full lineage information for each sequence is one of the output files of the analysis. In addition, a conflict file records cases in which a phylum-level taxonomy assignment disagreement exists within an OTU (such instances are rare and can reflect sequence misclassification within the greengenes database).


**Sample Assignment with BLAST:**


Taxonomy assignments are made by searching input sequences against a blast database of pre-assigned reference sequences. If a satisfactory match is found, the reference assignment is given to the input sequence. This method does not take the hierarchical structure of the taxonomy into account, but it is very fast and flexible. If a file of reference sequences is provided, a temporary blast database is built on-the-fly. The quality scores assigned by the BLAST taxonomy assigner are e-values.

To assign the sequences to the representative sequence set, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "blast_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt

Optionally, the user could changed the E-value ("-e"), using the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -e 0.01

**Assignment with the RDP Classifier:**

The RDP Classifier program (Wang, Garrity, Tiedje, & Cole, 2007) assigns taxonomies by matching sequence segments of length 8 to a database of previously assigned sequences. It uses a naive bayesian algorithm, which means that for each potential assignment, it attempts to calculate the probability of the observed matches, assuming that the assignment is correct and that the sequence segments are completely independent. The RDP Classifier is distributed with a pre-built database of assigned sequence, which is used by default. The quality scores provided by the RDP classifier are confidence values.

Note: If a reference set of sequences and taxonomy to id assignment file are provided, the script will use them to generate a new training dataset for the RDP Classifier on-the-fly.  Because of the RDP Classifier's implementation, all lineages in the training dataset must contain the same number of ranks.

To assign the representative sequence set, where the output directory is "rdp_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -m rdp

Alternatively, the user could change the minimum confidence score ("-c"), using the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -m rdp -c 0.85

**Sample Assignment with RTAX:**


Taxonomy assignments are made by searching input sequences against a fasta database of pre-assigned reference sequences. All matches are collected which match the query within 0.5% identity of the best match.  A taxonomy assignment is made to the lowest rank at which more than half of these hits agree.  Note that both unclustered read fasta files are required as inputs in addition to the representative sequence file.

To make taxonomic classifications of the representative sequences, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "rtax_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i rtax_repr_set_seqs.fasta -m rtax --read_1_seqs_fp read_1.seqs.fna --read_2_seqs_fp read_2.seqs.fna -r rtax_ref_seq_set.fna -t rtax_id_to_taxonomy.txt

**Sample Assignment with Mothur:**

The Mothur software provides a naive bayes classifier similar to the RDP Classifier.  A set of training sequences and id-to-taxonomy assignments must be provided.  Unlike the RDP Classifier, sequences in the training set may be assigned at any level of the taxonomy.

To make taxonomic classifications of the representative sequences, where the results are output to default directory "mothur_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i mothur_repr_set_seqs.fasta -m mothur -r mothur_ref_seq_set.fna -t mothur_id_to_taxonomy.txt


