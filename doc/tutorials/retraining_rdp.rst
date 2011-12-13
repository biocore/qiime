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

	make_otu_table.py -i otus/uclust_picked_otus/seqs_otus.txt -t otus/rdp_assigned_taxonomy_gg/seqs_rep_set_tax_assignments.txt -o otus/otu_table_gg.txt

That's it. The resulting OTU table (``otu_table_gg.txt``) can now be used in downstream analyses, such as `summarize_taxa_through_plots.py <../scripts/summarize_taxa_through_plots.html>`_.


Retraining RDP using a custom parameters file
=============================================

If you want to integrate retraining of the RDP classifier into your QIIME workflows, you can create a custom parameters file that can be used with the `pick_otus_through_otu_table.py <../scripts/pick_otus_through_otu_table.html>`_ workflow script. If the ``gg_otus_4feb2011`` directory is in ``/home/qiime/``, the values in your custom parameters file would be::

	assign_taxonomy:reference_seqs_fp	/home/qiime/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta
	assign_taxonomy:id_to_taxonomy_fp	/home/qiime/gg_otus_4feb2011/taxonomies/greengenes_tax_rdp_train.txt


Defining alternate training files
=================================

Training files can be defined by users for other taxonomies. The format is the same as the ``id_to_taxonomy_map`` used by the BLAST taxonomy assigner, defined `here <../documentation/file_formats.html#sequence-id-to-taxonomy-mapping-files>`_. You must provide this file as well as a fasta file of reference sequences where the identifiers correspond to the ids in the ``id_to_taxonomy_map``.

The RDP Classifier has several requirements about its taxonomy strings for retraining. The first entry in each taxonomy string must be ``Root``, and there must be exactly six levels, including ``Root``. For example, the first four lines in the ``4feb2011`` greengenes OTUs are::

	573145	Root;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae
	89440	Root;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae
	452783	Root;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae
	430240	Root;k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae

