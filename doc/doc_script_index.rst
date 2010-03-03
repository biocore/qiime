.. _scripts: Script Index

==================================
Scripts - Analyses and Parameters
==================================

General Notes
-------------

All QIIME analyses are performed using python (.py) scripts, which are located in the Qiime/scripts directory. To access QIIME python scripts, it may be useful to set an environment variable to the location of the innermost QIIME directory (the one containing `check_id_map.py <./scripts/check_id_map.html>`_, for example)::

	qdir=/path/to/QIIME/

Further commands of the form :file:`python QIIME_script.py -o option` can be invoked as :file:`python $qdir/QIIME_script.py -o option`. For all path description throughout the Documentation, the "/path/to/" refers to the physical location of each program on your local computer. For instance, the "/path/to/QIIME/" on my computer refers to "/Users/Jesse/Qiime/" on Mac OS X version 10.5.

The user can obtain help about the arguments which can be passed to each script, as well as examples and usage notes, by typing the following in the bash shell::

	python $qdir/script_of_interest.py –h

Workflow Scripts
----------------

.. toctree::
   :maxdepth: 1

   scripts/alpha_rarefaction 
   scripts/beta_diversity_through_3d_plots 
   scripts/jackknifed_upgma 
   scripts/pick_otus_through_otu_table 

User-Generated Mapping File
---------------------------

.. toctree::
   :maxdepth: 1

   scripts/check_id_map 
   scripts/merge_mapping_files 

Quality Checking
----------------

.. toctree::
   :maxdepth: 1

   scripts/denoise
   scripts/identify_chimeric_seqs 

Demultiplexing
--------------

.. toctree::
   :maxdepth: 1

   scripts/split_libraries 

Read OTU File (OTU Mapping)
---------------------------

.. toctree::
   :maxdepth: 1

   scripts/add_taxa 
   scripts/filter_otus_by_sample 
   scripts/merge_otu_maps
   scripts/parallel_pick_otus_blast 
   scripts/pick_otus 
   scripts/summarize_taxa 

Process Sequences
--------------------

.. toctree::
   :maxdepth: 1

   scripts/adjust_seq_orientation 
   scripts/exclude_seqs_by_blast 
   scripts/extract_seqs_by_sample_id
   scripts/pick_rep_set 

Sequence Alignment
------------------

.. toctree::
   :maxdepth: 1

   scripts/align_seqs 
   scripts/filter_alignment 
   scripts/parallel_align_seqs_pynast 

Tree-Building
-------------

.. toctree::
   :maxdepth: 1

   scripts/make_phylogeny 

Taxonomy Assignment
-------------------

.. toctree::
   :maxdepth: 1

   scripts/assign_taxonomy 
   scripts/parallel_assign_taxonomy_blast 
   scripts/parallel_assign_taxonomy_rdp

OTU Table Processing
--------------------

.. toctree::
   :maxdepth: 1

   scripts/filter_otu_table 
   scripts/filter_by_metadata
   scripts/make_otu_table 
   scripts/multiple_rarefactions
   scripts/otu_category_significance
   scripts/parallel_multiple_rarefactions 
   scripts/single_rarefaction 
   scripts/summarize_otu_by_cat 

Alpha-Diversity
---------------

.. toctree::
   :maxdepth: 1


   scripts/alpha_diversity 
   scripts/alpha_diversity_metrics
   scripts/collate_alpha 
   scripts/make_rarefaction_averages
   scripts/parallel_alpha_diversity 

Beta-Diversity
---------------

.. toctree::
   :maxdepth: 1

   scripts/beta_diversity_metrics
   scripts/beta_diversity 
   scripts/dissimilarity_mtx_stats 
   scripts/parallel_beta_diversity

Principal Coordinates Analysis (PCoA)
-------------------------------------

.. toctree::
   :maxdepth: 1

   scripts/principal_coordinates 
   scripts/transform_coordinate_matrices

Visualization
-------------

.. toctree::
   :maxdepth: 1

   scripts/cytoscape_usage
   scripts/compare_3d_plots
   scripts/make_2d_plots 
   scripts/make_3d_plots 
   scripts/make_3d_plot_prefs_file
   scripts/make_distance_histograms 
   scripts/make_pie_charts 
   scripts/make_otu_heatmap_html 
   scripts/make_otu_network
   scripts/make_rarefaction_plots 
   scripts/preferences_file

Jackknifing
-----------

.. toctree::
   :maxdepth: 1

   scripts/make_bootstrapped_tree 
   scripts/tree_compare 
   scripts/upgma_cluster 

SRA Submission
--------------

.. toctree::
   :maxdepth: 1

   scripts/make_fastq 
   scripts/make_library_id_lists 
   scripts/make_per_library_sff 
   scripts/make_sra_submission 
   scripts/per_library_stats 
   scripts/process_sff 
   scripts/sra_spreadsheet_to_map_files
   scripts/trim_sff_primers 

Utilities
---------

.. toctree::
   :maxdepth: 1

   scripts/blast_wrapper 
   scripts/cluster_jobs
   scripts/convert_unifrac_sample_mapping_to_otu_table 
   scripts/fix_arb_fasta 
   scripts/make_qiime_py_file
   scripts/make_qiime_rst_file
   scripts/parallel_blast 
   scripts/poller
   scripts/poller_example
   scripts/print_qiime_config
