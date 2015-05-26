.. _pick_rep_set:

.. index:: pick_rep_set.py

*pick_rep_set.py* -- Pick representative set of sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

After picking OTUs, you can then pick a representative set of sequences. For each OTU, you will end up with one sequence that can be used in subsequent analyses. By default, the representative sequence for an OTU is chosen as the most abundant sequence showing up in that OTU. This is computed by collapsing identical sequences, and choosing the one that was read the most times as the representative sequence (note that each of these would have a different sequence identifier in the FASTA provided as input).


**Usage:** :file:`pick_rep_set.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_file
		Path to input otu mapping file [REQUIRED]
	
	**[OPTIONAL]**
		
	-f, `-`-fasta_file
		Path to input fasta file [REQUIRED if not picking against a reference set; default: None]
	-m, `-`-rep_set_picking_method
		Method for picking representative sets.  Valid choices are random, longest, most_abundant, first [default: first (first chooses cluster seed when picking otus with uclust)]
	-o, `-`-result_fp
		Path to store result file [default: <input_sequences_filepath>_rep_set.fasta]
	-l, `-`-log_fp
		Path to store log file [default: No log file created.]
	-s, `-`-sort_by
		Sort by otu or seq_id [default: otu]
	-r, `-`-reference_seqs_fp
		Collection of preferred representative sequences [default: None]


**Output:**

The output from `pick_rep_set.py <./pick_rep_set.html>`_ is a single FASTA file containing one sequence per OTU. The FASTA header lines will be the OTU identifier (from here on used as the unique sequence identifier) followed by a space, followed by the sequence identifier originally associated with the representative sequence. The name of the output FASTA file will be <input_sequences_filepath>_rep_set.fasta by default, or can be specified via the "-o" parameter.



**Simple example: picking a representative set for de novo-picked OTUs:**

The script `pick_rep_set.py <./pick_rep_set.html>`_ takes as input an 'OTU map' (via the "-i" parameter) which maps OTU identifiers to sequence identifiers. Typically, this will be the output file provided by `pick_otus.py <./pick_otus.html>`_. Additionally, a FASTA file is required, via "-f", which contains all of the sequences whose identifiers are listed in the OTU map.  By default, a representative sequence will be chosen as the most abundant sequence in the OTU. This can be changed to, for example, choose the first sequence listed in each OTU by passing -m first.

::

	pick_rep_set.py -i seqs_otus.txt -f seqs.fna -o rep_set1.fna

**Picking OTUs with "preferred representative" sequences:**

Under some circumstances you may have a fasta file of "preferred representative" sequences. An example of this is if you were to pick OTUs against a reference collection with uclust_ref. In this case you may want your representative sequences to be the sequences from the reference collection, rather than the sequences from your sequencing run. To achieve this, you can pass the original reference collection via -r. If you additionally allowed for new clusters (i.e., sequences which don't match a reference sequence are used as seeds for new OTUs) you'll also need to pass the original sequence collection to pick a representative sequence from the sequencing run in that case.

::

	pick_rep_set.py -i seqs_otus.txt -f seqs.fna -r refseqs.fasta -o rep_set2.fna


