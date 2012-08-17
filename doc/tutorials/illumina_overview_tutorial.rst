.. _illumina_overview_tutorial:

==========================
Illumina Overview Tutorial
==========================

This tutorial covers a variety of QIIME analyses beginning with Illumina sequencing data. It expands on the `Processing Illumina Data <>`_ tutorial, which focused primarily on demultiplexing Illumina data. This tutorial is is intended to be quick to run, and as a result is based on a small subset of a full Illumina Genome Analyzer II (GAIIx) run. 

For the sake of runtime we additionally take one big shortcut. In this tutorial, we perform *closed-reference* OTU picking, where reads are searched against a reference collection (Greengenes in this case), and reads that do not hit the reference collection are discarded. This works well for well-characterized environments such as human-associated habitats, but not as well for environments that harbor a lot of novel diversity with respect to the reference collection (e.g., the `Guerrero Negro microbial mats <>`_ are an example of where this does not work well). In general, we *highly* recommend that you use our `subsampled open-reference OTU picking protocol <>`_ for processing Illumina data.

The data used in this analysis is derived from the `Moving Pictures of the Human Microbiome <>`_ study, where two human subjects collected daily samples from four body sites: the tongue, the palm of the left hand, the palm of the right hand, and the gut (via fecal samples obtained by swapping used toilet paper). This data was sequenced across six lanes of an Illumina GAIIx, using the barcoding amplicon sequencing protocol described in `Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample <>`_. A more recent version of this protocol that can be used with the Illumina HiSeq 2000 and MiSeq can be found `here <>`_. 

Obtaining the data
------------------

The data used in this tutorial is available for download `here <>`_. You can pull this data to your system by running the following command::

	wget XXX

.. note:: MacOS does not come with ``wget`` pre-installed. You can either install this using your preferred method, or just download the dataset directly by clicking the above link. 

Once you have this data on your system you can unzip it as follows::

	tar -xzf moving_pictures_tutorial.tgz

Next, change to the directory that was just unzipped. We'll work there for this tutorial.
::
	
	cd moving_pictures_tutorial

Demultiplexing and quality filtering
------------------------------------

We start by demultiplexing our sequences (or assigning barcoded reads to the samples they are derived from). In general, you should get fastq files for your sequence and barcode reads as separate files. Additionally, on the multiple-lane Illumina platforms, we typically reuse barcodes across lanes, so we must demultiplex each lane independently. To do that, run the following commands::

	split_libraries_fastq.py -o slout/ -i subsampled_fastq/subsampled_s_1_sequence.fastq,subsampled_fastq/subsampled_s_2_sequence.fastq,subsampled_fastq/subsampled_s_3_sequence.fastq,subsampled_fastq/subsampled_s_4_sequence.fastq,subsampled_fastq/subsampled_s_5_sequence.fastq,subsampled_fastq/subsampled_s_6_sequence.fastq -b subsampled_fastq/subsampled_s_1_sequence_barcodes.fastq,subsampled_fastq/subsampled_s_2_sequence_barcodes.fastq,subsampled_fastq/subsampled_s_3_sequence_barcodes.fastq,subsampled_fastq/subsampled_s_4_sequence_barcodes.fastq,subsampled_fastq/subsampled_s_5_sequence_barcodes.fastq,subsampled_fastq/subsampled_s_6_sequence_barcodes.fastq -m filtered_mapping_l1.txt,filtered_mapping_l2.txt,filtered_mapping_l3.txt,filtered_mapping_l4.txt,filtered_mapping_l5.txt,filtered_mapping_l6.txt

This is a big command, but it's relatively straight-forward. We're telling QIIME that we have six lanes of sequence data (specified as a comma-separated list of files passed as ``-i``), six lanes of barcode data (specified as a comma-separated list of files passed as ``-b``), and a metadata mapping file corresponding to each lane (specified as a comma-separated list of files passed as ``-m``). The metadata mapping file contains the sample-to-barcode mapping that we need for demultiplexing. **Important**: The order of files passed for ``-m``, ``-b``, and ``-i`` must be consistent, so if you pass the lane 1 sequence data first for ``-i``, you must pass the lane 1 barcode data first for ``-b``, and the lane 1 metadata mapping file first as ``-m``. The only other parameter here is the output directory, which we'll call ``slout``, for *split libraries output*.

When this completes you can view the files that were created by calling ``ls`` on the output directory::
	
	ls slout/

You'll see that we now have three files. The ones we care the most about are ``split_library_log.txt``, and ``seqs.fna``. The former provides a summary of what was filtered during quality filtering, and the latter is the demultiplexed sequence data, combined across all lanes.

Review the ``split_library_log.txt`` file as follows::

	less split_library_log.txt

The command that we're using here, ``less``, allows us to view files in read-only mode. Use the up and down arrows to scroll through the file. When you're done, hit the ``q`` key to exit from ``less``. 

If you want to see how many sequences remain after demultiplexing and quality filtering you can call the following QIIME command::
	
	count_seqs.py -i slout/seqs.fna

You can use ``count_seqs.py`` to count the sequences and summarize the sequence lengths in any ``fasta`` or ``fastq`` file.

OTU picking
-----------

Now that we have demultiplexed sequences, we're ready to cluster these sequences into OTUs. As described above, in the interest of providing a tutorial that can be run quickly for educational purposes, we're using a closed-reference OTU picking protocol here (in general, you'll want to use open-reference OTU picking, as discussed `here <>`_). For closed-reference OTU picking we use `pick_reference_otus_through_otu_table.py`::

	pick_reference_otus_through_otu_table.py -o ucrC_fast/ -i slout/seqs.fna -r $reference_seqs -t $reference_tax -p ucrC_fast_params.txt

Note that this command takes the ``seqs.fna`` file that was generated in the previous step, as well as the reference fasta file (``$reference_seqs`` here) and the taxonomies associated with the reference sequences (``$reference_tax`` here). We're also taking on additional shortcut here for the sake of reduced run time. We're using the *fast uclust* parameters. As with most tools for sequence clustering, it is possible to trade performance (in terns of runtime) for accuracy. To allow this to run in a just a couple of minutes, we're using parameters that are optimized for reduced runtime. These correspond to ``uclust``'s default parameters. QIIME uses slightly more stringent parameter settings by default. These parameters are specified the the *parameters file* which is passes as ``-p``. You can find information on defining parameters files `here <>`_.

The primary output that we can about from this command is the *OTU table*, or the counts of the number of times each OTU is observed in each sample. QIIME uses the Genomics Standards Consortium *candidate standard* Biological Observation Matrix (BIOM) format for representing these files. You can find additional information on the `BIOM format here <http://www.biom-format.org>`_, and information on converting this files to tab-separated text that can be view in spreadsheet programs `here <>`_. 

To see some summary statistics of the OTU table we can run the following command::

	per_library_stats.py -i ucrC_fast/otu_table.biom

Because we started with six lanes of data but have now summarized these in a single OTU table, we need to merge the per-lane mapping files into a single *combined* mapping file that represents all six lanes, and therefore all of our data. Note that we will have duplicated barcodes in our mapping file, but that's OK as we've already demultiplexed our reads. We don't use the barcodes again. We can merge the six mapping files as follows::

	merge_mapping_files.py -o combined_mapping_file.txt -m filtered_mapping_l1.txt,filtered_mapping_l2.txt,filtered_mapping_l3.txt,filtered_mapping_l4.txt,filtered_mapping_l5.txt,filtered_mapping_l6.txt

From this point on, we'll work with ``combined_mapping_file.txt``.

The OTU table is a key piece of data, and essentially all of the additional analyses that you'll want to do with QIIME use that as input. We'll now explore some of the additional analyses. As these all branch from the OTU table, it's not necessary to run this in order. 

Comparing microbial communities: Beta diversity
-----------------------------------------------

Now that we have an OTU table, and we're working with a reference phylogenetic tree for our analysis, we can compute UniFrac distances between our samples. To do this, we use one of QIIME's *workflow* scripts, which computes a beta diversity distance matrix containing distances between all samples, summarizes that distance matrix using Principal Coordinates Analysis (PCoA), and then generates PCoA plots. You can run this workflow as follows::

	beta_diversity_through_plots.py -o bdiv_even258/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt -t $reference_tree -e 258

The parameters here are mostly straight-forward: we're passing our OTU table as ``-i``, our metadata mapping file as ``-m``, our phylogenetic tree as ``-t``, and the output directory as ``-o``. The one additional parameter we're passing here is ``-e``, which is the even sampling depth that we want to apply in this analysis. This is extremely important: in order to accurately compare our microbial communities with UniFrac, each sample must have the same number of sequences: otherwise we may see samples cluster by their depth of sequencing coverage, which is not representative of the biology of the samples, but rather a technical artifact. ``-e 258`` tells QIIME to randomly subsample each of the samples in the OTU table to exactly 258 sequences per sample, without replacement. The importance of this step is discussed `here <>`_.


Generating taxonomic summaries of microbial communities
-------------------------------------------------------

We can additionally generate taxonomic summaries of these samples using the ``summarize_taxa_through_plots.py`` script. These can be run on a per-sample basis as follows::

	summarize_taxa_through_plots.py -o taxa_summaries/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt

After this command completes, there will be two ``html`` files in the new ``taxa_summaries`` directory: ``bar_charts.html`` and ``area_charts.html``. For categorical data the bar charts are generally more informative, and for continuous data the area charts are generally more informative. If working with continuous data you will likely want to call `sort_otu_table.py <>`_ first, sorting by the continuous variable in your metadata file. The `filter_samples_from_otu_table.py <>`_ script may also be useful here, to filter out samples that you may not want in your taxa summary plot (e.g., control samples, or human skin samples if you're trying to generate a plot illustrating the change in your human gut communities over time).

You may alternatively be interesting in a taxonomic summary of your samples collapsed by some metadata category. For example, in this data set collapsing by the sample type (left palm, right palm, tongue, and gut) is a useful way to see the differences across community types. We can achieve this by adding the ``-c`` parameter to our call to ``summarize_taxa_through_plots.py``. Here we collapse by ``SampleType``, which is a column header in our mapping file::

	summarize_taxa_through_plots.py -o taxa_summaries_by_SampleType/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt -c "SampleType"

As before, be can view either bar charts or area charts by opening the corresponding file.

Generating alpha rarefaction plots
----------------------------------

Alpha rarefaction plots are a useful way to compare the relative alpha diversities of our samples, and also to determine if we are approaching complete coverage of our microbial communities. We can generate alpha rarefaction plots with QIIME as follows::

	alpha_rarefaction.py -o arare_max258/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt -t $reference_tree -e 258

Notice that we again pass ``-e 258`` here. In this case, this specifies the maximum rarefaction depth that ... **pick up here, too tired to keep going now...**

Next steps
----------

Link to Procrustes tutorial











