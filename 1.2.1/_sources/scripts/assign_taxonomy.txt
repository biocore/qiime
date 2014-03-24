.. _assign_taxonomy:

.. index:: assign_taxonomy.py

*assign_taxonomy.py* -- Assign taxonomy to each sequence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

Contains code for assigning taxonomy, using several techniques.

Given a set of sequences, assign_taxonomy attempts to assign the taxonomy of each sequence. Currently there are two methods implemented: assignment with BLAST and assignment with the RDP classifier. The output of this step is a mapping of input sequence identifiers (1st column of output file) to taxonomy (2nd column) and quality score (3rd column). The sequence identifier of the best BLAST hit is also included if the blast method is used (4th column). 


**Usage:** :file:`assign_taxonomy.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_fasta_fp
		Path to the input fasta file
	
	**[OPTIONAL]**
		
	-t, `-`-id_to_taxonomy_fp
		Path to tab-delimited file mapping sequences to assigned taxonomy. Each assigned taxonomy is provided as a semicolon-separated list. For assignment with rdp, each assigned taxonomy must be exactly 6 levels deep. [default: None; REQUIRED when method is blast]
	-r, `-`-reference_seqs_fp
		Path to reference sequences.  For assignment with blast, these are used to generate a blast database. For assignment with rdp, they are used as training sequences for the classifier.[default: None; REQUIRED if -b is not provided when method is blast]
	-p, `-`-training_data_properties_fp
		Path to ".properties" file in pre-compiled training data for the RDP Classifier.  This option is overridden by the -t and -r options. [default: None]
	-m, `-`-assignment_method
		Taxon assignment method [default:rdp]
	-b, `-`-blast_db
		Database to blast against.  Must provide either --blast_db or --reference_seqs_db for assignment with blast [default: None]
	-c, `-`-confidence
		Minimum confidence to record an assignment, only used for rdp method [default: 0.8]
	-e, `-`-e_value
		Maximum e-value to record an assignment, only used for blast method [default: 0.001]
	-o, `-`-output_dir
		Path to store result file [default: <ASSIGNMENT_METHOD>_assigned_taxonomy]


**Output:**

The consensus taxonomy assignment implemented here is the most detailed lineage description shared by 90% or more of the sequences within the OTU (this level of agreement can be adjusted by the user). The full lineage information for each sequence is one of the output files of the analysis. In addition, a conflict file records cases in which a phylum-level taxonomy assignment disagreement exists within an OTU (such instances are rare and can reflect sequence misclassification within the greengenes database).



Example of consensus lineage: 

The OTU containing 5 sequences annotated as shown below would be assigned to the "Desulfovibrionaceae" level because only 80% of sequences agree with the "LE30" annotation.

* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; 
* Bacteria; Proteobacteria; Desulfovibrionales; Desulfovibrionaceae; LE30

Assignments are provided in a two column tab-delimited format, which maps input sequence identifiers to assignments. Each assignment is specified as a list of taxa separated by a ';' character.

Example of an assignment output file:

======== =================================================================
AY800210 Archaea;Euryarchaeota;Halobacteriales;uncultured 
EU883771 Archaea;Euryarchaeota;Methanomicrobiales;Methanomicrobium et rel.
EF503699 Archaea;Crenarchaeota;uncultured;uncultured 
DQ260310 Archaea;Euryarchaeota;Methanobacteriales;Methanobacterium 
EF503697 Archaea;Crenarchaeota;uncultured;uncultured
======== =================================================================


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

To assign the representative sequence set, where the output directory is "rdp_assigned_taxonomy", the you can run the following command:


::

	assign_taxonomy.py -i repr_set_seqs.fasta -m rdp

Alternatively, the user could change the minimum confidence score ("-c"), using the following command:

::

	assign_taxonomy.py -i repr_set_seqs.fasta -m rdp -c 0.85

Note: If a reference set of sequences and taxonomy to id assignment file are provided, the script will use them to generate a new training dataset for the RDP Classifier on-the-fly. Due to limitations in the generation of a training set, each provided assignment must contain exactly 6 taxa in the following order: domain (level=2), phylum (level=3), class (level=4), order (5), family (level=6), and genus (level=7). Additionally, each genus name must be unique, due to the internal algorithm used by the RDP Classifier.



