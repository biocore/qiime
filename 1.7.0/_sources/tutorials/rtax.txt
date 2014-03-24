.. _rtax:

=======================================================================
Taxonomic classifications of single- and paired-end sequences with RTAX
=======================================================================

The RTAX procedure takes advantage of mate-pair information when performing taxonomic classification.  The additional information from a second read may allow a more precise taxonomy assignment to be made.  Conversely, the second read may reveal greater uncertainty in the assignment than would have been inferred from one read alone, leading to a less precise (but more accurate) assignment.

RTAX may also be used to classify single reads.  Here too it attempts to avoid overprediction by considering a range of best hits (something like a k-nearest-neighbor classifier), in contrast to the single best hit approach used by the Blast classifier.

Taxonomic classifications of paired-end sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To obtain taxonomic classifications of paired-end sequences with RTAX, follow these steps:

Obtain fasta files for each read
--------------------------------

Run `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ on each read individually (making sure to reverse-complement the reverse read)::

    split_libraries_fastq.py -i s_1_1.seqs.fastq -o sl_out.1/ -b s_1_1.barcodes.fastq -m s_1_map.txt

    split_libraries_fastq.py -i s_1_3.seqs.fastq -o sl_out.2/ -b s_1_3.barcodes.fastq -m s_1_map.txt --rev_comp

Note that while most amplicons should be represented in both the forward and the reverse read files (i.e., with mate pairs matchable by amplicon ID), it may happen that one read is eliminated by quality filters while the other is not.

Normally the fasta files resulting from `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ will contain millions of sequences; for tutorial purposes we provide tiny samples taken from (`Caporaso et al. 2010 PNAS <http://www.ncbi.nlm.nih.gov/pubmed/20534432>`_).  These are constructed to contain 10,000 mate-pairs, as well as 1000 sequences represented only by the forward read and 2000 represented only by the reverse read.  The sample forward and reverse read files thus contain 11,000 and 12,000 sequences, respectively.

Perform OTU picking
-------------------

Our procedure is to perform OTU picking on one read only, but then to obtain additional information from the second read at the taxonomic classification step.  So, start by running the standard one-read pipeline::

    pick_de_novo_otus.py -i rtax_sample_data/forward_read.11k.fna -o pick_otus.out -v

In addition to choosing representative sequences, that produces taxonomy information using just the given read and the default classifier; this may be interesting for comparison with the RTAX results. (The PyNast part may fail on short sequences, given the default options, but the OTU representatives will be written to pick_otus.out/rep_set, which is really all we need). You can choose OTUs based on the forward or the reverse read here, or do it both ways to compare; the downstream RTAX classification doesn't care which read is which.

The OTU-picking procedure is completely unaware of the reverse read, and so may choose representative sequences for which only the forward read was available.  In this example, there are 2301 cluster representatives, of which 2079 turn out to have a mate pair and 222 do not.

Classify cluster representatives
--------------------------------

To classify the cluster representatives, pass them as the main input to assign_taxonomy as usual, but also provide the complete unclustered fasta files for both reads::

    assign_taxonomy.py -i pick_otus.out/rep_set/forward_read.11k_rep_set.fasta -m rtax --read_1_seqs_fp rtax_sample_data/forward_read.11k.fna --read_2_seqs_fp rtax_sample_data/reverse_read.12k.fna -r  /software/gg_otus-4feb2011-release/rep_set/gg_97_otus_4feb2011.fasta -t  /software/gg_otus-4feb2011-release/taxonomies/greengenes_tax.txt -v

The resulting taxonomy file will contain the same FASTA headers that were present in the representative sequence file, typically a cluster ID and the read ID produced by the `split_libraries_fastq.py <../scripts/split_libraries_fastq.html>`_ step.

RTAX does not presently compute a confidence value for each assignment.  Because Qiime expects one, RTAX simply reports 1.0 for every sequence-- but this is not meaningful.

Impact of reference databases
-----------------------------

The choice of reference database can have a large impact on the results.  The 97%-clustered database provided with QIIME is used in this example.  For 99% OTUs, simply replace "97" with "99" in the -r option.

An alternative set of reference databases--those used in (`Soergel et al. 2012 <http://www.ncbi.nlm.nih.gov/pubmed/22237546>`_)--is available `here <http://dev.davidsoergel.com/rtax>`_.  To use these, the command line options would include ``-t rtax.greengenes.20110311/gg.nr.taxonomy -r rtax.greengenes.20110311/gg.nr.fasta``.

Fallback from paired-end to single-ended classification
-------------------------------------------------------

By default, RTAX classifies only those sequences for which mate pairs are available.  In this case, 2079 of the 2301 cluster representatives have an associated reverse read and thus are classified, while the remaining 222 are ignored.

If you wish to classify the remaining sequences using the single-ended classification procedure, pass the "--single_ok" option::

	assign_taxonomy.py -i pick_otus.out/rep_set/forward_read.11k_rep_set.fasta -m rtax --read_1_seqs_fp rtax_sample_data/forward_read.11k.fna --read_2_seqs_fp rtax_sample_data/reverse_read.12k.fna -r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -t gg_otus_4feb2011/taxonomies/greengenes_tax.txt -v --single_ok

This classifies all of the cluster representative sequences, with the caveat that the taxonomic predictions for the single-ended representatives may be less accurate than those for the paired-end cases.

Taxonomic classifications of single-ended sequences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RTAX can be used for single-ended classification as well; simply omit the "--read_2_seqs_fp" option.  The original read file for the single read must nonetheless be passed via "--read_1_seqs_fp"; it is needed to establish the mapping between cluster ID, read representative ID, and amplicon ID.::

	assign_taxonomy.py -i pick_otus.out/rep_set/forward_read.11k_rep_set.fasta -m rtax --read_1_seqs_fp rtax_sample_data/forward_read.11k.fna -r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -t gg_otus_4feb2011/taxonomies/greengenes_tax.txt -v
