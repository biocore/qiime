.. _retraining_rdp:

==============================
Re-training the RDP Classifier
==============================

Retraining the RDP classifier and assign taxonomy
=================================================

This tutorial covers how to retrain the RDP Classifier with an alternate taxonomy to use the RDP Classifier with arbitrary taxonomies. This is useful, for example, to assign greengenes taxonomy strings to your sequences, or to assign taxonomy to eukaryotic sequences using the Silva database.

This tutorial will illustrate an example where you've run the `QIIME Overview Tutorial <../tutorials/tutorial.html>`_, and then want to re-assign taxonomy using the greengenes taxonomy. To do this you'll need the greengenes reference OTUs. This is covered in the first step.

This tutorial assumes that you've already run the `QIIME Overview Tutorial <../tutorials/tutorial.html>`_, and that you're working in the directory where you ran the tutorial commands.

Getting the greengenes reference OTUs
-------------------------------------

The most recent version of the greengenes OTUs is always listed on the top right corner of the `QIIME Blog homepage <http://blog.qiime.org>`_ (click the "Most recent Greengenes OTUs" link on that page). As of this writing that is the ``4feb2011`` version, so we'll illustrate commands working with that. 

Download and unzip the greengenes reference OTUs::

	wget http://greengenes.lbl.gov/Download/Sequence_Data/Fasta_data_files/Caporaso_Reference_OTUs/gg_otus_4feb2011.tgz
	tar -xvzf gg_otus_4feb2011.tgz

Re-train the RDP Classifier and classify the representative sequences
---------------------------------------------------------------------

Next you'll retrain the RDP classifier and classify your sequences. You can use either `assign_taxonomy.py <../scripts/assign_taxonomy.html>`_ or `parallel_assign_taxonomy_rdp.py <../scripts/parallel_assign_taxonomy_rdp.html>`_ for this.

::

	assign_taxonomy.py -i otus/rep_set/seqs_rep_set.fasta -t gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train.txt -r gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta -o otus/rdp_assigned_taxonomy_gg/
	
Integrate the taxonomy assignments into a new OTU table
-------------------------------------------------------

Next, you'll rebuild the OTU table with the new taxonomic information.

::

	make_otu_table.py -i otus/uclust_picked_otus/seqs_otus.txt -t otus/rdp_assigned_taxonomy_gg/seqs_rep_set_tax_assignments.txt -o otus/otu_table_gg.biom

That's it. The resulting OTU table (``otu_table_gg.biom``) can now be used in downstream analyses, such as `summarize_taxa_through_plots.py <../scripts/summarize_taxa_through_plots.html>`_.


Retraining RDP using a custom parameters file
=============================================

If you want to integrate retraining of the RDP classifier into your QIIME workflows, you can create a custom parameters file that can be used with the `pick_de_novo_otus.py <../scripts/pick_de_novo_otus.html>`_ workflow script. If the ``gg_otus_4feb2011`` directory is in ``$HOME/``, the values in your custom parameters file would be::

	assign_taxonomy:reference_seqs_fp	$HOME/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta
	assign_taxonomy:id_to_taxonomy_fp	$HOME/gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train.txt


Defining alternate training files
=================================

Training files can be defined by users for other taxonomies. The format is the same as the ``id_to_taxonomy_map`` used by the BLAST taxonomy assigner, defined `here <../documentation/file_formats.html#sequence-id-to-taxonomy-mapping-files>`_. You must provide this file as well as a fasta file of reference sequences where the identifiers correspond to the ids in the ``id_to_taxonomy_map``.

The RDP Classifier has several requirements about its taxonomy strings for retraining.  The first column of this tab separated file is the sequence identifiers (see the reference sequence file below).  The second column is the taxonomy strings in descending order of taxonomic specification, separated by semicolons.  The number of taxonomic levels must be equal for every line.

An example set of four lines in the ``4feb2011`` greengenes OTUs that are valid for RDP retraining are::

	573145	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia;s__
	89440	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia;s__
	222043	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Raoultella;s__Raoultellaornithinolytica
	430240	k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Serratia;s__Serratiamarcescens
	
The reference sequence file should have fasta labels that match all of the labels as listed in the taxonomy mapping file.  Orientation of the sequences does not matter for RDP, but as this can impact other software such as uclust, so it is suggested that the sequences be in the same orientation to avoid complications with other QIIME scripts.

An example fasta file (with truncated nucleotide sequences) that matches the above taxonomy strings is::

	>573145
	AGAGTTTGATCATGGCTCAGATTGAACGCAGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCT
	>89440
	AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGC
	>222043
	AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGTAGCACAGAAAGCTTACTC
	>101567
	TGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGAGGAGGAAGGCATTAAGGTTAATAACCTTAGTGATTGA
	
Tips for building and troubleshooting custom RDP retraining files
=================================================================

1.  Always use a plain-text editor when modifying taxonomy mapping or reference sequence files to avoid the addition of unwanted characters.  If errors occur, it may be necessary to check a file using a command such as ``less`` in the terminal to check for hidden characters.
2.  Always have the same number of taxonomic levels (separated by semicolons) in the taxonomy mapping file.
3.  The memory requirements can change when alternative files are used.  Additionally memory can be allocated with the --rdp_max_memory parameter when calling assign_taxonomy.py.
4.  Avoid having white space in the taxonomy mapping data.  In the above example taxonomy mapping data, if the first line had an extra space before the class level, an error would be raised.
5.  Avoid empty levels in the taxonomy mapping with double semicolons.  At an unknown taxonomic level, use a consistent naming convention (such as the `s__` listed above for unknown species), rather than leaving a level empty.
6.  Used consistent capitalization.  For instance, do not have one class level named c__Gammaproteobacteria and another named c__gammaproteobacteria, as these would be considered distinct taxonomies.
7.  Reference sequences should be unaligned (no gap or leading/trailing characters such as . or -) with nucleotide characters only.
8.  As the exact line causing an error is sometimes difficult to detect, it may be advisable to break apart the taxonomy mapping file and reference sequences into smaller subsets to determine which part(s) are causing error(s).
