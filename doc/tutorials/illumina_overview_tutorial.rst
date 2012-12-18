.. _illumina_overview_tutorial:

===============================================================================
Illumina Overview Tutorial: closed reference OTU picking and diversity analyses
===============================================================================

This tutorial covers a variety of QIIME using Illumina sequencing data. It expands on the :ref:`Processing Illumina Data <processing_illumina_data>` tutorial, which focused primarily on demultiplexing Illumina data. This tutorial is is intended to be quick to run, and as such, uses only a subset of a full Illumina Genome Analyzer II (GAIIx) run. 

To further reduce runtime, we take a big shortcut: in this tutorial, we perform *closed-reference* OTU picking, where reads are searched against a reference collection (Greengenes in this case), and reads that do not hit the reference collection are discarded. This works well for well-characterized environments such as human-associated habitats, but not as well for environments that harbor a lot of novel diversity with respect to the reference collection (e.g., the `Guerrero Negro microbial mats <http://www.ncbi.nlm.nih.gov/pubmed/22832344>`_ are an example of where this does not work well). In general, we *highly* recommend that you use our `subsampled open-reference OTU picking protocol <open_reference_illumina_processing.html#option-2-subsampled-open-reference-otu-picking>`_ for processing Illumina data.

The data used in this analysis is derived from the `Moving Pictures of the Human Microbiome <http://www.ncbi.nlm.nih.gov/pubmed/21624126>`_ study, where two human subjects collected daily samples from four body sites: the tongue, the palm of the left hand, the palm of the right hand, and the gut (via fecal samples obtained by swapping used toilet paper). This data was sequenced across six lanes of an Illumina GAIIx, using the barcoding amplicon sequencing protocol described in `Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample <http://www.ncbi.nlm.nih.gov/pubmed/20534432>`_. A more recent version of this protocol that can be used with the Illumina HiSeq 2000 and MiSeq can be found `here <http://www.ncbi.nlm.nih.gov/pubmed/22402401>`_. 

IPython Notebook
----------------

The steps in this tutorial are also presented in an IPython Notebook which you can view `here <http://nbviewer.ipython.org/urls/raw.github.com/qiime/qiime/master/examples/ipynb/illumina_overview_tutorial.ipynb>`_. 

We'd like to thank the IPython developers for their help with using their tools with QIIME. For more information on using QIIME with IPython, see `our recent paper here <http://www.nature.com/ismej/journal/vaop/ncurrent/full/ismej2012123a.html>`_. You can find more information on the `IPython Notebook here <http://ipython.org/ipython-doc/stable/interactive/htmlnotebook.html>`_, and the `nbviewer tool (which we use to display the notebook) here <http://nbviewer.ipython.org/>`_.

Obtaining the data
------------------

The data used in this tutorial is available for download: `moving_pictures_tutorial.tgz <https://s3.amazonaws.com/qiime-tutorial/moving_pictures_tutorial.tgz>`_. You can pull this data to your system by running the following command::

	wget https://s3.amazonaws.com/qiime-tutorial/moving_pictures_tutorial.tgz

.. note:: If you are using an Apple computer, we recommend downloading the dataset directly using the above link as OS X does not come with ``wget`` pre-installed. Alternatively, for the adventurous, ``wget`` can be compiled from `source <ftp://ftp.gnu.org/gnu/wget/>`_ or you can use ``curl``, which is similar to ``wget`` although the usage is different.

Once you have this data on your system you can unzip it as follows::

	tar xzf moving_pictures_tutorial.tgz

Defining reference filepaths with environment variables
-------------------------------------------------------

Through-out this tutorial we make use of a reference sequence collection, tree, and taxonomy derived from the Greengenes database. As these files may be store in different locations on your system, we'll define them as environment variables using the paths as they would be if you're running in a QIIME virtual machine (e.g., on AWS or with the Virtual Box). We'll then reference the environment variables through-out this tutorial when they are used. If you're not working on either of these systems, you'll have to modify these paths. Run the following::

	export reference_seqs=/home/ubuntu/qiime_software/gg_otus-12_10-release/rep_set/97_otus.fasta
	export reference_tree=/home/ubuntu/qiime_software/gg_otus-12_10-release/trees/97_otus.tree
	export reference_tax=/home/ubuntu/qiime_software/gg_otus-12_10-release/taxonomy/97_otu_taxonomy.txt


Demultiplexing and quality filtering
------------------------------------

We start by demultiplexing our sequences (i.e. assigning barcoded reads to the samples they are derived from). In general, you should get seperate fastq files for your sequence and barcode reads. On the multiple-lane Illumina platforms, we typically reuse barcodes across lanes, so we must demultiplex each lane independently. To do that, run the following command (*will run for a few minutes*)::

	split_libraries_fastq.py -o slout/ -i moving_pictures_tutorial/subsampled_fastq/subsampled_s_1_sequence.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_2_sequence.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_3_sequence.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_4_sequence.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_5_sequence.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_6_sequence.fastq -b moving_pictures_tutorial/subsampled_fastq/subsampled_s_1_sequence_barcodes.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_2_sequence_barcodes.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_3_sequence_barcodes.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_4_sequence_barcodes.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_5_sequence_barcodes.fastq,moving_pictures_tutorial/subsampled_fastq/subsampled_s_6_sequence_barcodes.fastq -m moving_pictures_tutorial/filtered_mapping_l1.txt,moving_pictures_tutorial/filtered_mapping_l2.txt,moving_pictures_tutorial/filtered_mapping_l3.txt,moving_pictures_tutorial/filtered_mapping_l4.txt,moving_pictures_tutorial/filtered_mapping_l5.txt,moving_pictures_tutorial/filtered_mapping_l6.txt

This is a big command, but it's relatively straight-forward. We're telling QIIME that we have six lanes of sequence data (specified as a comma-separated list of files passed as ``-i``), six lanes of barcode data (specified as a comma-separated list of files passed as ``-b``), and a metadata mapping file corresponding to each lane (specified as a comma-separated list of files passed as ``-m``). The metadata mapping file contains the sample-to-barcode mapping that we need for demultiplexing. **Important**: The order of files passed for ``-m``, ``-b``, and ``-i`` must be consistent, so if you pass the lane 1 sequence data first for ``-i``, you must pass the lane 1 barcode data first for ``-b``, and the lane 1 metadata mapping file first as ``-m``. The only other parameter here is the output directory, which we'll call ``slout``, for *split libraries output*.

When this completes you can view the files that were created by calling ``ls`` on the output directory::
	
	ls slout/

You'll see that we now have three files. The ones we care the most about are ``split_library_log.txt``, and ``seqs.fna``. The former provides a summary of what was filtered during quality filtering, and the latter is the demultiplexed sequence data, combined across all lanes.

Review the ``split_library_log.txt`` file as follows::

	less slout/split_library_log.txt

The command that we're using here, ``less``, allows us to view files in read-only mode. Use the up and down arrows to scroll through the file. When you're done, hit the ``q`` key to exit from ``less``. 

If you want to see how many sequences remain after demultiplexing and quality filtering you can call the following QIIME command::
	
	count_seqs.py -i slout/seqs.fna

You can use ``count_seqs.py`` to count the sequences and summarize the sequence lengths in any ``fasta`` or ``fastq`` file.

OTU picking
-----------

Now that we have demultiplexed sequences, we're ready to cluster these sequences into OTUs. As mentioned above, in the interest of providing a tutorial that can be run quickly for educational purposes, we're using a closed-reference OTU picking protocol here, although typically you'll want to use open-reference OTU picking, as discussed `here <open_reference_illumina_processing.html>`_). For closed-reference OTU picking we use `pick_reference_otus_through_otu_table.py` (*will run for a few minutes*)::

	pick_reference_otus_through_otu_table.py -o ucrC_fast/ -i slout/seqs.fna -r $reference_seqs -t $reference_tax -p moving_pictures_tutorial/ucrC_fast_params.txt

Note that this command takes the ``seqs.fna`` file that was generated in the previous step, as well as the reference fasta file (``$reference_seqs`` here) and the taxonomies associated with the reference sequences (``$reference_tax`` here). We're also taking on an additional shortcut here for the sake of reduced run time: we're using the *fast uclust* parameters. To allow this to run in a just a couple of minutes, we're using parameters that are optimized for reduced runtime at the expense of accuracy. These correspond to ``uclust``'s default parameters. QIIME uses slightly more stringent parameter settings by default. These parameters are specified the the *parameters file* which is passes as ``-p``. You can find information on defining parameters files `here <../documentation/file_formats.html#qiime-parameters>`_.

The primary output that we can about from this command is the *OTU table*, or the number of times each operational taxonomic unit (OTU) is observed in each sample. QIIME uses the Genomics Standards Consortium *candidate standard* Biological Observation Matrix (BIOM) format for representing these files. You can find additional information on the `BIOM format here <http://www.biom-format.org>`_, and information on converting this files to tab-separated text that can be view in spreadsheet programs `here <http://biom-format.org/documentation/biom_conversion.html>`_. 

To see some summary statistics of the OTU table we can run the following command::

	per_library_stats.py -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom

We started with six lanes of data but have now summarized these in a single OTU table. However, we still need to merge the per-lane mapping files into a single *combined* mapping file that represents all six lanes, and therefore all of our data. Note that we will have duplicated barcodes in our mapping file, but that's OK as we've already demultiplexed our reads. We don't use the barcodes again. We can merge the six mapping files as follows::

	merge_mapping_files.py -o combined_mapping_file.txt -m moving_pictures_tutorial/filtered_mapping_l1.txt,moving_pictures_tutorial/filtered_mapping_l2.txt,moving_pictures_tutorial/filtered_mapping_l3.txt,moving_pictures_tutorial/filtered_mapping_l4.txt,moving_pictures_tutorial/filtered_mapping_l5.txt,moving_pictures_tutorial/filtered_mapping_l6.txt

From this point on, we'll work with ``combined_mapping_file.txt``.

The OTU table is a key piece of data, and essentially all of the additional analyses that you'll want to do with QIIME use that as input. We'll now explore some of the additional analyses. As these all branch from the OTU table, it's not necessary to run this in order. 

Comparing microbial communities: beta diversity
-----------------------------------------------

Now that we have an OTU table, and we're working with a reference phylogenetic tree for our analysis, we can compute UniFrac distances between our samples. To do this, we will use one of QIIME's *workflow* scripts, which computes a beta diversity distance matrix containing distances between all samples, summarizes that distance matrix using Principal Coordinates Analysis (PCoA), and then generates PCoA plots. You can run this workflow as follows (*will run for a few minutes*)::

	beta_diversity_through_plots.py -o bdiv_even258/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt -t $reference_tree -e 258

The parameters used are described as follows: we're passing our OTU table as ``-i``, our metadata mapping file as ``-m``, our phylogenetic tree as ``-t``, the output directory as ``-o`` and last, ``-e`` to specify an even sampling depth that we want to apply in this analysis. The sampling depth is extremely important: in order to accurately compare our microbial communities with UniFrac, each sample must have the same number of sequences otherwise we may see samples cluster by their depth of sequencing coverage, which is not representative of the biology of the samples, but rather a technical artifact. ``-e 258`` tells QIIME to randomly subsample each of the samples in the OTU table to exactly 258 sequences per sample, without replacement.

 .. warning:: If you're working on a remote system (e.g., EC2) and want to download the results of this analysis for viewing, you'll need to download the whole directory for the plots to be viewable. You can zip this directory (``tar -czf bdiv_even258.tgz bdiv_even258``) and then `download it using Cyberduck <./working_with_aws.html#working-with-cyberduck>`_ or `via the command line <http://qiime.org/tutorials/working_with_aws.html#working-with-command-line-tools>`_. 


Generating taxonomic summaries of microbial communities
-------------------------------------------------------

We can additionally generate taxonomic summaries of these samples using the ``summarize_taxa_through_plots.py`` script. These can be run on a per-sample basis as followsi (*will run for a few minutes*)::

	summarize_taxa_through_plots.py -o taxa_summaries/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt

After this command completes, there will be two ``html`` files in the new ``taxa_summaries`` directory: ``bar_charts.html`` and ``area_charts.html``. For categorical data the bar charts are generally more informative, and for continuous data the area charts are generally more informative. If working with continuous data you will likely want to call `sort_otu_table.py <../scripts/sort_otu_table.html>`_ first, sorting by the continuous variable in your metadata file. The `filter_samples_from_otu_table.py <../scripts/filter_samples_from_otu_table.html>`_ script may also be useful here to filter out samples that you may not want in your taxa summary plot (e.g., control samples, or human skin samples if you're trying to generate a plot illustrating the change in your human gut communities over time).

You may alternatively be interesting in a taxonomic summary of your samples collapsed by some metadata category. For example, in this data set collapsing by the sample type (left palm, right palm, tongue, and gut) is a useful way to see the differences across community types. We can achieve this by adding the ``-c`` parameter to our call to ``summarize_taxa_through_plots.py``. Here we collapse by ``SampleType``, which is a column header in our mapping file::

	summarize_taxa_through_plots.py -o taxa_summaries_by_SampleType/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt -c "SampleType"

As before, be can view either bar charts or area charts by opening the corresponding file.

 .. warning:: If you're working on a remote system (e.g., EC2) and want to download the results of this analysis for viewing, you'll need to download the whole directory for the plots to be viewable. You can zip this directory (``tar -czf taxa_summaries.tgz taxa_summaries``) and then `download it using Cyberduck <./working_with_aws.html#working-with-cyberduck>`_ or `via the command line <http://qiime.org/tutorials/working_with_aws.html#working-with-command-line-tools>`_. 

Comparing microbial communities: alpha diversity
-------------------------------------------------

**WARNING: This step can be require approximately 20 minutes to run.**

Alpha rarefaction plots are a useful way to compare the relative alpha diversities across samples, and also to determine if we are approaching complete coverage of our microbial communities. We can generate alpha rarefaction plots with QIIME as follows (*will run for over 10 minutes*)::

	alpha_rarefaction.py -o arare_max258/ -i ucrC_fast/uclust_ref_picked_otus/otu_table.biom -m combined_mapping_file.txt -t $reference_tree -e 258

Notice that we again pass ``-e 258`` here. In this case, this specifies the maximum rarefaction depth: in general you want to choose the same value as specified for the even sampling depth to `beta_diversity_through_plots.py` if you are interested in looking at alpha diversity and rarefaction by metadata category.

 .. warning:: If you're working on a remote system (e.g., EC2) and want to download the results of this analysis for viewing, you'll need to download the whole directory for the plots to be viewable. You can zip this directory (``tar -czf arare_max258.tgz arare_max258``) and then `download it using Cyberduck <./working_with_aws.html#working-with-cyberduck>`_ or `via the command line <http://qiime.org/tutorials/working_with_aws.html#working-with-command-line-tools>`_. 

Next steps
----------

This illustrates some of the basic features of QIIME, and there are a lot of places to go from here. If you're interested in seeing additional visualizations, you should check out the `QIIME overview tutorial <tutorial.html>`_. We also highly recommend reviewing how to perform open-reference OTU picking on Illumina data, which you can find `here <open_reference_illumina_processing.html#option-2-subsampled-open-reference-otu-picking>`_. The `Procrustes analysis tutorial <procrustes_analysis.html>`_ illustrates a really cool analysis, allowing you to continue with the same data used here, comparing against the samples sequenced on 454 (rather than Illumina, as in this analysis). If you're interested in some possibilities for statistical analyses you can try our `supervised learning <running_supervised_learning.html>`_ or `distance matrix comparison <distance_matrix_comparison.html>`_ tutorials, both of which can be adapted to use data generated in this tutorial.

Modified Procrustes Analysis Steps (temporary)
----------------------------------------------

We're in the process of modifying the `Procrustes analysis tutorial <procrustes_analysis.html>`_ to more directly follow from this one. In the meantime, these commands will allow you to continue::

	pick_reference_otus_through_otu_table.py -i moving_pictures_tutorial/subsampled_454_seqs.fna -o 454_ucrC_fast/ -r $reference_seqs -t $reference_tax -p moving_pictures_tutorial/ucrC_fast_params.txt
	per_library_stats.py -i 454_ucrC_fast/uclust_ref_picked_otus/otu_table.biom
	beta_diversity_through_plots.py -o bdiv_even135/ -i 454_ucrC_fast/uclust_ref_picked_otus/otu_table.biom -e 135 -t $reference_tree -m moving_pictures_tutorial/454_map.txt
	transform_coordinate_matrices.py -o 454_v_illumina/ -i bdiv_even258/unweighted_unifrac_pc.txt,bdiv_even135/unweighted_unifrac_pc.txt -s moving_pictures_tutorial/procrustes_sid_map.txt -r 100
	compare_3d_plots.py -o 454_v_illumina/plots/ -i 454_v_illumina/pc1_transformed.txt,454_v_illumina/pc2_transformed.txt -m moving_pictures_tutorial/procrustes_metadata_map.txt --custom_axes days_since_epoch



.. warning:: If you're working on a remote system (e.g., EC2) and want to download the results of this analysis for viewing, you'll need to download the whole directory for the plots to be viewable. You can zip this directory (``tar -czf 454_v_illumina.tgz 454_v_illumina``) and then `download it using Cyberduck <./working_with_aws.html#working-with-cyberduck>`_ or `via the command line <http://qiime.org/tutorials/working_with_aws.html#working-with-command-line-tools>`_. 














