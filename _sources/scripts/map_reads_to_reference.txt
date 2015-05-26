.. _map_reads_to_reference:

.. index:: map_reads_to_reference.py

*map_reads_to_reference.py* --  Script for performing assignment of reads against a reference database 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

 


**Usage:** :file:`map_reads_to_reference.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_seqs_filepath
		Path to input sequences file
	-r, `-`-refseqs_fp
		Path to reference sequences to search against [default: None]
	
	**[OPTIONAL]**
		
	-m, `-`-assignment_method
		Method for picking OTUs.  Valid choices are: bwa-short, usearch, bwa-sw, blat, blat-nt. [default: usearch]
	-t, `-`-observation_metadata_fp
		Path to observation metadata (e.g., taxonomy, EC, etc) [default: None]
	-o, `-`-output_dir
		Path to store result file [default: ./<METHOD>_mapped/]
	-e, `-`-evalue
		Max e-value to consider a match [default: 1e-10]
	-s, `-`-min_percent_id
		Min percent id to consider a match, expressed as a fraction between 0 and 1 [default: 0.75]
	`-`-genetic_code
		ID of genetic code to use for DNA translations (please see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) Only valid with -m blat. [default: 11]
	`-`-max_diff
		MaxDiff to consider a match (applicable for -m bwa-short) -- see the aln section of "man bwa" for details [default (defined by bwa): 0.04]
	`-`-queryalnfract
		Min percent of the query seq that must match to consider a match, expressed as a fraction between 0 and 1 (usearch only) [default: 0.35]
	`-`-targetalnfract
		Min percent of the target/reference seq that must match to consider a match, expressed as a fraction between 0 and 1 (usearch only) [default: 0.0]
	`-`-max_accepts
		Max_accepts value (usearch only) [default: 1]
	`-`-max_rejects
		Max_rejects value to (usearch only) [default: 32]


**Output:**

 


Run assignment with usearch using default parameters

::

	map_reads_to_reference.py -i query_nt.fasta -r refseqs_pr.fasta

Run nucleotide versus protein BLAT using default parameters

::

	map_reads_to_reference.py -i query_nt.fasta -r refseqs_pr.fasta -m blat

Run nucleotide versus protein BLAT using scricter e-value threshold

::

	map_reads_to_reference.py -i query_nt.fasta -r refseqs_pr.fasta -o blat_mapped_strict/ -e 1e-70  -m blat

Run nucleotide versus nucleotide BLAT with default parameters

::

	map_reads_to_reference.py -i query_nt.fasta -r refseqs_nt.fasta -m blat-nt

Run assignment with bwa-short using default parameters. bwa-short is intended to be used for reads up to 200bp. WARNING: reference sequences must be dereplicated! No matches will be found to reference sequences which show up multiple times (even if their sequence identifiers are different)!

::

	map_reads_to_reference.py -i query_nt.fasta -r refseqs_nt.fasta -m bwa-short

Run assignment with bwa-sw using default parameters.  WARNING: reference sequences must be dereplicated! No matches will be found to reference sequences which show up multiple times (even if their sequence identifiers are different)!

::

	map_reads_to_reference.py -i query_nt.fasta -r refseqs_nt.fasta -m bwa-sw


