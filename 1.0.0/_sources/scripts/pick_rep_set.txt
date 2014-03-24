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
	-f, `-`-fasta_file
		Path to input fasta file [REQUIRED]
	
	**[OPTIONAL]**
		
	-m, `-`-rep_set_picking_method
		Method for picking representative sets.  Valid choices are random, longest, most_abundant, first[default: most_abundant]
	-o, `-`-result_fp
		Path to store result file [default: <input_sequences_filepath>_rep_set.fasta]
	-l, `-`-log_fp
		Path to store log file [default: No log file created.]
	-s, `-`-sort_by
		Sort by otu or seq_id [default: otu]


**Output:**

The output from `pick_rep_set.py <./pick_rep_set.html>`_ is a single FASTA file containing one sequence per OTU. The FASTA header lines will be the OTU identifier (from here on used as the unique sequence identifier) followed by a space, followed by the sequence identifier originally associated with the representative sequence. The name of the output FASTA file will be <input_sequences_filepath>_rep_set.fasta by default, or can be specified via the "-o" parameter.



**Simple example:**

The script `pick_rep_set.py <./pick_rep_set.html>`_ takes as input an 'OTU file' (via the "-i" parameter) which maps OTU identifiers to sequence identifiers. Typically, this will be the output file provided by `pick_otus.py <./pick_otus.html>`_. Additionally, a FASTA file is required, via "-f", which contains all of the sequences whose identifiers are listed in the OTU file. The following command shows an example of this where the resulting file is output to the directory "repr_set/" and default parameters were used (choose most abundant, sort by OTU id and do not write a log file):

::

	pick_rep_set.py -i seqs_otus.txt -f seqs.fna -o repr_set/

**Random selection example:**

Alternatively, if the user would like to choose the sequence by random "-m random" and then sort by the sequence identifier ("-s seq_id"), they could use the following command:

::

	pick_rep_set.py -i seqs_otus.txt -f seqs.fna -o repr_set/ -m random -s seq_id


