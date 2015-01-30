.. _assign_taxonomy:

.. index:: assign_taxonomy.py

*assign_taxonomy.py* -- Assign taxonomy to each sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Contains code for assigning taxonomy, using several techniques.

Given a set of sequences, assign_taxonomy.py attempts to assign the taxonomy of each sequence. Currently the methods implemented are assignment with BLAST, the RDP classifier, RTAX, mothur, and uclust. The output of this step is an observation metadata mapping file of input sequence identifiers (1st column of output file) to taxonomy (2nd column) and quality score (3rd column). There may be method-specific information in subsequent columns.

Reference data sets and id-to-taxonomy maps for 16S rRNA sequences can be found in the Greengenes reference OTU builds. To get the latest build of the Greengenes OTUs (and other marker gene OTU collections), follow the "Resources" link from http://qiime.org. After downloading and unzipping you can use the following files as -r and -t, where <otus_dir> is the name of the new directory after unzipping the reference OTUs tgz file.

-r <otus_dir>/rep_set/97_otus.fasta
-t <otus_dir>/taxonomy/97_otu_taxonomy.txt



**Usage:** :file:`assign_taxonomy.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-t, `-`-id_to_taxonomy_fp
		Path to tab-delimited file mapping sequences to assigned taxonomy. Each assigned taxonomy is provided as a semicolon-separated list. For assignment with rdp, each assigned taxonomy must be exactly 6 levels deep. [default: /Users/jairideout/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt]
	-r, `-`-reference_seqs_fp
		Path to reference sequences.  For assignment with blast, these are used to generate a blast database. For assignment with rdp, they are used as training sequences for the classifier. [default: /Users/jairideout/.virtualenvs/qiime/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta]
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
		Used to parse the result of OTU clustering, to get the read_1_id for each clusterID.  The clusterID itself is assumed to be the first field, and is not captured by the regex.  (used for RTAX only). [default: \S+\s+(\S+)]
	`-`-amplicon_id_regex
		Used to parse the result of split_libraries, to get the ampliconID for each read_1_id.  Two groups capture read_1_id and ampliconID, respectively.  (used for RTAX only). [default: (\S+)\s+(\S+?)\/]
	`-`-header_id_regex
		Used to parse the result of split_libraries, to get the portion of the header that RTAX uses to match mate pairs.  The default uses the amplicon ID, not including /1 or /3, as the primary key for the query sequences.  Typically this regex will be the same as amplicon_id_regex, except that only the second group is captured.  (used for RTAX only). [default: \S+\s+(\S+?)\/]
	-m, `-`-assignment_method
		Taxon assignment method, must be one of rdp, blast, rtax, mothur, uclust, sortmerna [default: uclust]
	`-`-sortmerna_db
		Pre-existing database to search against when using sortmerna [default: None]
	`-`-sortmerna_e_value
		Maximum E-value when clustering [default = 1.0]
	`-`-sortmerna_coverage
		Mininum percent query coverage (of an alignment) to consider a hit, expressed as a fraction between 0 and 1 [default: 0.9]
	`-`-sortmerna_best_N_alignments
		This option specifies how many best alignments per read will be written [default: 5]
	`-`-sortmerna_threads
		Specify number of threads to be used for sortmerna mapper which utilizes multithreading. [default: 1]
	-b, `-`-blast_db
		Database to blast against.  Must provide either --blast_db or --reference_seqs_db for assignment with blast [default: None]
	-c, `-`-confidence
		Minimum confidence to record an assignment, only used for rdp and mothur methods [default: 0.5]
	`-`-min_consensus_fraction
		Minimum fraction of database hits that must have a specific taxonomic assignment to assign that taxonomy to a query, only used for sortmerna and uclust methods [default: 0.51]
	`-`-similarity
		Minimum percent similarity (expressed as a fraction between 0 and 1) to consider a database match a hit, only used for sortmerna and uclust methods [default: 0.9]
	`-`-uclust_max_accepts
		Number of database hits to consider when making an assignment, only used for uclust method [default: 3]
	`-`-rdp_max_memory
		Maximum memory allocation, in MB, for Java virtual machine when using the rdp method.  Increase for large training sets [default: 4000]
	-e, `-`-blast_e_value
		Maximum e-value to record an assignment, only used for blast method [default: 0.001]
	-o, `-`-output_dir
		Path to store result file [default: <ASSIGNMENT_METHOD>_assigned_taxonomy]


**Output:**

The consensus taxonomy assignment implemented here is the most detailed lineage description shared by 90% or more of the sequences within the OTU (this level of agreement can be adjusted by the user). The full lineage information for each sequence is one of the output files of the analysis. In addition, a conflict file records cases in which a phylum-level taxonomy assignment disagreement exists within an OTU (such instances are rare and can reflect sequence misclassification within the greengenes database).


**Assign taxonomy with the uclust consensus taxonomy assigner (default):**

Perform database search with uclust to retrive up to uclust_max_accepts hits for each query sequence. Then assign the most specific taxonomic label that is associated with at least min_consensus_fraction of the hits.

::

	assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt

**Assignment with SortMeRNA:**

Perform database search with sortmerna to retrieve up to sortmerna_best_N_alignments hits for each query sequence. Then assign the most specific taxonomic label that is associated with at least min_consensus_fraction of the hits. 

::

	assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -m sortmerna

**Assignment with BLAST:**


Taxonomy assignments are made by searching input sequences against a blast database of pre-assigned reference sequences. If a satisfactory match is found, the reference assignment is given to the input sequence. This method does not take the hierarchical structure of the taxonomy into account, but it is very fast and flexible. If a file of reference sequences is provided, a temporary blast database is built on-the-fly. The quality scores assigned by the BLAST taxonomy assigner are e-values.

To assign the sequences to the representative sequence set, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "blast_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -m blast

Optionally, the user could changed the E-value ("-e"), using the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -r ref_seq_set.fna -t id_to_taxonomy.txt -e 0.01 -m blast

**Assignment with the RDP Classifier:**

The RDP Classifier (Wang, Garrity, Tiedje, & Cole, 2007) assigns taxonomies using Naive Bayes classification. By default, the classifier is retrained using the values provided for --id_to_taxonomy_fp and --reference_seqs_fp.

::

	assign_taxonomy.py -i repr_set_seqs.fasta -m rdp

Assignment with the RDP Classifier using an alternative minimum confidence score by passing -c:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -m rdp -c 0.80

**Assignment with RTAX:**


Taxonomy assignments are made by searching input sequences against a fasta database of pre-assigned reference sequences. All matches are collected which match the query within 0.5% identity of the best match.  A taxonomy assignment is made to the lowest rank at which more than half of these hits agree.  Note that both unclustered read fasta files are required as inputs in addition to the representative sequence file.

To make taxonomic classifications of the representative sequences, using a reference set of sequences and a taxonomy to id assignment text file, where the results are output to default directory "rtax_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i rtax_repr_set_seqs.fasta -m rtax --read_1_seqs_fp read_1.seqs.fna --read_2_seqs_fp read_2.seqs.fna -r rtax_ref_seq_set.fna -t rtax_id_to_taxonomy.txt

**Assignment with Mothur:**

The Mothur software provides a naive bayes classifier similar to the RDP Classifier.A set of training sequences and id-to-taxonomy assignments must be provided.  Unlike the RDP Classifier, sequences in the training set may be assigned at any level of the taxonomy.

To make taxonomic classifications of the representative sequences, where the results are output to default directory "mothur_assigned_taxonomy", you can run the following command:

::

	assign_taxonomy.py -i mothur_repr_set_seqs.fasta -m mothur -r mothur_ref_seq_set.fna -t mothur_id_to_taxonomy.txt


