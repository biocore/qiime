.. _rtax:

===========================================================
Taxonomic classifications of paired-end sequences with RTAX
===========================================================

To obtain taxonomic classifications of paired-end sequences with RTAX, follow these steps:


Run split_libraries_fastq on each read individually (making sure to reverse-complement the reverse read)::

    split_libraries_fastq.py -i s_1_1.seqs.fastq -o sl_out.1/ -b s_1_1.barcodes.fastq -m s_1_map.txt

    split_libraries_fastq.py -i s_1_3.seqs.fastq -o sl_out.2/ -b s_1_3.barcodes.fastq -m s_1_map.txt

Our procedure is to perform OTU picking on one read only, but then to obtain additional information from the second read at the taxonomic classification step.
So, start by running the standard one-read pipeline::

    pick_otus_through_otu_table.py -i sl_out.1/s_1_1.seqs.fna -o pick_otus.out -v

In addition to choosing representative sequences, that produces taxonomy information using just the given read and the default classifier; this may be interesting for comparison with the RTAX results.
(The PyNast part may fail on short sequences, given the default options, but the OTU representatives will be written to pick_otus.out/rep_set, which is really all we need).
You can choose OTUs based on the forward or the reverse read here, or do it both ways to compare; the downstream RTAX classification doesn't care which read is which.

Then pass the cluster representatives as the main input to assign_taxonomy, but also provide the complete unclustered fasta files for both reads::

    assign_taxonomy.py -i pick_otus.out/rep_set/s_1_1.seqs_rep_set.fasta -m rtax --read_1_seqs_fp sl_out.1/s_1_1.seqs.fna --read_2_seqs_fp sl_out.2/s_1_3.seqs.fna -r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -t gg_otus_4feb2011/taxonomies/greengenes_tax.txt -v

The choice of reference database can have a large impact on the results.  The 97%-clustered database provided with QIIME is used in this example.  For 99% OTUs, simply replace "97" with "99" in the -r option.

An alternative set of reference databases--those used in (Soergel et al. 2012)--is available `here <http://dev.davidsoergel.com/rtax`>.
To use these, the command line options would include ``-t rtax.greengenes.20110311/gg.nr.taxonomy -r rtax.greengenes.20110311/gg.nr.fasta``.

The resulting taxonomy file will contain the same FASTA headers that were present in the representative sequence file, typically a cluster ID and the read ID produced by the split_libraries_fastq step.





