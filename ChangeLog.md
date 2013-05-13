QIIME 1.6.0-dev (changes since QIIME 1.6.0 go here)
===================================================
* Added section for using USEARCH 6.1 chimera checking with ``identify_chimeric_seqs.py`` to "Chimera checking sequences with QIIME" tutorial.
* ``compare_alpha_diversity.py`` is changed so that the output now includes average alpha diversity values as well as the comparison p and t vals. 
* ``compare_distance_matrices.py`` has a new option ``--variable_size_distance_classes`` for running Mantel correlogram over distance classes that vary in size (i.e. width) but contain the same number of pairwise distances in each class.
* ``qiime.filter.sample_ids_from_category_state_coverage`` now supports splitting on a category.
* Modified add_qiime_labels script to use normal metadata mapping file with a column specified for fasta file names to make more consistent with other scripts.
* otu_category_significance.py now makes better use of the BIOM Table API, addressing a performance issue when using CSMat as the sparse backend.
* Required biom-format version is now 1.1.2.
* Added qiime.group.get_adjacent_distances, which is useful for plotting distances between "adjacent" sample ids in a list provided by the user. This is useful, for example, in plotting distances between adjacent temporal samples in a time series.
* core_qiime_analyses.py has been replaced with core_diversity_analyses.py. This follows a re-factoring to support only "downstream" analyses (i.e., starting with a BIOM table). This makes the script more widely applicable as it's now general to any BIOM data and/or different OTU picking strategies.
* Fixed a bug in make_3d_plots.py related to biplot calculations. This bug would change the placement of taxonomic groups based on how many taxa were included in the biplot analysis. Examples and additional details can be found here: [#677](https://github.com/qiime/qiime/issues/677).
* Major refactoring of workflow tests and organization of workflow code. The workflow library code and tests have now been split apart into separate files. This makes it a lot more manageable, which will support a more general refactoring of the workflow code in the future to make it easier to develop new workflows. The workflow tests have also been updated to use the new test data described in [#582](https://github.com/qiime/qiime/issues/582), which is now accessible through ``qiime.test. get_test_data()`` and ``qiime.test.get_test_data_fps()``. This provides improved testing of boundary cases in each workflow, as well as more consistent tests across the workflows.
* otu_category_significance.py now supports an input directory of biom tables, and can write out either a single collated results file or an individual file for every input table in the directory. The -o output_fp is now a required parameter rather than an optional parameter.
* simsam.py now has a -m/--mapping_fp option and writes output to a directory instead of a single file. -n/--num and -d/--dissim now accept a single number or comma-separated list of values.
* supervised_learning.py can now handle input directorys of otu tables, can write a single collated results file if the input directory is of rarefied otu tables, and the -o output fp option is now a required parameter.
* The qiime_test_data repository has been merged into the main qiime repository, which will facilitate development by not requiring users to time pull requests against two repositories. Users will no longer have to specify qiime_test_data_dir in their qiime_config files to include the script usage tests in runs of all_tests.py. all_tests.py will now know how to find qiime_test_data, and will run all of the script usage tests by default.
* pick_reference_otus_through_otu_table.py now outputs otu_table.biom in top-level output directory rather than nested in the otu picking output directory.
* pick_reference_otus_through_otu_table.py has been renamed pick_closed_reference_otus.py (issue #708).
* pick_subsampled_reference_otus_through_otu_table.py has been renamed pick_open_reference_otus.py (issue #708).
* pick_otus_through_otu_table.py has been renamed pick_de_novo_otus.py (issue #708).
* make_distance_comparison_plots.py now supports auto-sizing of distribution plots via --distribution_width (which is the new default) and better handles numeric label types with very large or small ranges (e.g. elevation) by scaling x-axis units to [1, (number of data points)]. --group_spacing has been removed in favor of the new auto-sizing feature.
* per_library_stats.py removed in favor of biom-format's print_biom_table_summary.py.
* Add SourceTracker tutorial, and changed QIIME to depend on SourceTracker 0.9.5 (which is modified to facilitate use with QIIME).
* Moran's I (in compare_categories.py) now supports identical samples (i.e. zeros in the distance matrix that aren't on the diagonal).
* summarize_taxa.py now outputs taxa summary tables in both classic (TSV) and BIOM formats by default. This will allow taxa summary tables to be used with other QIIME scripts that expect BIOM files as input. This change is the first step towards adding full support for BIOM taxon tables in QIIME. summarize_taxa.py also has two new options: --suppress_classic_table_output and  --supress_biom_table_output.
* make_distance_boxplots.py and make_distance_comparison_plots.py now explicitly state the alternative hypothesis used in the t-tests.
* parallel_blast.py now has a different option for providing a blast db (--blast_db). This implies that the current --refseqs_path should be used only for providing a fasta file of reference sequences. The --suppress_format_blastdb option has been removed since it is no longer needed.

QIIME 1.6.0 (18 Dec 2012)
=========================
* Added ``filter_taxa_from_otu_table.py`` to support filtering OTUs with (or without) specific taxonomy assignments from an OTU table.
* Added parameters to ``pick_subsampled_reference_otus_through_otu_table.py`` to suppress taxonomy assignment (``--suppress_taxonomy_assignment``), and alignment and tree building steps (``--suppress_align_and_tree``). These are useful for cases where a taxonomy may not exist for the reference collection (not too common) or when the region doesn't work well for phylogenetic reconstruction (e.g., fungal ITS). Additionally fixed a bug where alternate ``assign_taxonomy`` parameters provided in the parameters file would be ignored when running in parallel.
* Detrending of quadratic curvature in ordination coordinates now a feature of QIIME. This approach was used in [Harris JK, et al. "Phylogenetic stratigraphy in the Guerrero Negro hypersaline microbial mat."](http://www.nature.com/ismej/journal/v7/n1/full/ismej201279a.html).
* Supervised learning mislabeling output now includes binary "mislabeled" columns at 5%, 10%, ..., 95%, 99%.
* Added tutorial on Fungal ITS analysis.
* Added tutorial on predicting mislabeled samples.
* Modified the parameters (de novo chimera detection, reference chimera detection, and size filtering) for USEARCH options with ``pick_otus.py`` to ``suppress_X`` and ``False`` by default, rather than ``True`` and turned off by calling, to make them more intuitive to use and work better with the workflow scripts.
* Added a ``simpson_reciprocal`` measure of alpha diversity, which is ``1/D``, following the [definition here](http://www.countrysideinfo.co.uk/simpsons.htm) among other places. Note the measure ``reciprocal_simpson`` is ``1/simpson``, not ``1/D``. It was removed for clarity.
* Added new script, ``compute_core_microbiome.py``, which identifies the core OTUs (i.e., those defined in some user-defined percentage of the samples).
* Major refactoring of parallel QIIME. Repetitive code was consolidated into the ParallelWrapper class, which may ultimately move to PyCogent. The only script interface changes are that the ``-Y/--python_exe_fp``, ``-N (serial script filepath)``, and ``-P/--poller_fp`` parameters are no longer available to the user. These were very infrequently (if ever) modified from defaults, so it doesn't make sense to continue to support these. These changes will allow for easier development of new parallel wrappers and facilitate changes to the underlying parallel functionality. 
* Added new script, ``compare_taxa_summaries.py``, and supporting library and test code (``qiime/compare_taxa_summaries.py`` and ``tests/test_compare_taxa_summaries.py``) to allow for the comparison of taxa summary files, including sorting and filling, expected, and paired comparisons using pearson or spearman correlation. Added accompanying tutorial (``doc/tutorials/taxa_summary_comparison.rst``).
* New script for parallel trie otu picker.
* Made ``loaddata.r`` more robust when making mapping files, distance matrices, etc. compatible with each other. There were rare cases that caused some R functions (e.g. ``betadisper``) to fail if empty levels were left in the parsed mapping file.
* Fixed issue in ``ParallelWrapper`` class that could have caused a deadlock if run from within a subprocess with pipes.
* ``make_distance_boxplots.py`` and ``make_distance_comparison_plots.py`` can now perform Student's two-sample t-tests to determine whether a pair of boxplots/distributions are significantly different (using both parametric and nonparametric Monte Carlo-based tests of significance). These changes include three new options to the two scripts (``--tail_type``, ``--num_permutations``, and ``--suppress_significance_tests``), as well as a new function ``all_pairs_t_test`` in ``qiime.stats``. The accompanying tutorial has also been updated to cover the new statistical tests.
* Checks are now in place to prevent asymmetric and non-hollow distance matrices from being used in ``make_distance_boxplots.py``, ``make_distance_comparison_plots.py``, ``make_distance_histograms.py``, ``compare_categories.py``, and ``compare_distance_matrices.py``. The relevant script help and underlying library code has been documented to warn against their use, and the symmetry checks can be easily disabled if performance becomes an issue in the future.
* ``qiime.util.DistanceMatrix`` has new method ``is_symmetric_and_hollow``.
* Added the new Illumina Overview Tutorial which was developed for the ISME 14 Bioinformatics Workshop and added the IPython notebook files that were used in the ISME 14 workshop under the new ``examples/ipynb`` directory. These can be used by changing to the ``ipynb`` directory and running ``ipython notebook`` on a system with IPython and the IPython Notebook dependencies installed. Also moved the ``qiime_tutorial`` directory to the new ``examples`` folder.
* Added support for translated database mapping through ``map_reads_to_reference.py`` and ``parallel_map_reads_to_reference.py`` and related library code, parallel code, etc. This is analogous to closed-reference OTU picking, but can translate queries so is useful for mapping metagenomic or metatranscriptomic data against databases of functional genes (e.g., KEGG). Currently BLAT and usearch are supported for translated searching.
* ``qiime.util.qiime_system_call`` now has an optional shell parameter that is passed through to ``subprocess.Popen``.
* Changed ``compare_categories.py`` script interface such that ``--method rda`` is no longer supported and must now be ``--method dbrda`` as the method we provide is db-RDA (capscale), not traditional RDA; added the ability to pass the number of permutations (``-n``) for PERMDISP and db-RDA (these were previously not supported); updated script documentation, statistical method descriptions, and accompanying tutorial to be of overall better quality and clarity; output filename when method is PERMDISP is now ``permdisp_results.txt`` instead of ``betadisper_results.txt``, which is consistent with the rest of the methods; significant refactor of underlying code to be better tested and maintained easier; added better error checking and handling for the types of categories that are accepted by the statistical methods (e.g. checking that categories are numeric if they need to be, making sure categories do not contain all unique values, or a single value); fixed output format for BEST method to be easier to read and consistent with the other methods; ``qiime.util.MetadataMap`` class has a few new utility methods to suppport some of these changes.
* ``compare_alpha_diversity.py`` now supports both parametric and nonparametric two sample t-tests (nonparametric is the default) with the new optional options ``-t/--test_type`` and ``-n/--num_permutations``. Also fixed a bug that used the wrong degrees of freedom in the t-tests, yielding incorrect t statistics and p-values, and added correction for multiple comparisons.
* Removed tree method ``raxml`` from ``make_phylogeny.py``'s choices for ``-t/--tree_method``. Tree method ``raxml_v730`` should now be used instead. RAxML v703 is no longer supported.
* Minimum PyNAST version requirement upgraded to PyNAST 1.2.
* ``make_distance_boxplots.py``, ``make_distance_comparison_plots.py``, and ``make_distance_histograms.py`` now correctly output TSV data files with ``.txt`` extension instead of ``.xls`` (this allows them to be opened easier in programs such as Excel).
* ``make_distance_boxplots.py`` has a new option ``--color_individual_within_by_field`` that allows the "individual within" boxplots to be optionally colored to indicate their membership in another mapping file field. A legend is also included.
* Added ``sample_ids_from_category_state_coverage`` function to ``qiime/filter.py`` to support filtering of samples based on a subject's category coverage. For example, this function is useful for filtering individuals out of a time series study that do not meet some sort of timepoint coverage criteria.
* ``assign_taxonomy.py`` now supports assignment with tax2tree version 1.0 and mothur version 1.25.0.
* Added new script ``load_remote_mapping_file.py`` and accompanying tutorial to allow exporting and downloading of mapping files stored as Google Spreadsheets.
* Fixed bug in ``parallel_assign_taxonomy_blast.py`` which would cause the script to hang if a relative path was passed for ``-o``.
* Added the [``qiime_test_data``](https://github.com/qiime/qiime_test_data) repository which contains example input and output for most QIIME scripts. The individual script documentation was completely refactored so that usage examples correspond to the example input and output files. The *basic script testing* functionality was removed from ``all_tests.py`` and replaced with more detailed testing of the scripts based on their usage examples.
* ``add_taxa.py`` was removed in favor of ``add_metadata.py`` (a ``biom-format`` project script). See the new [tutorial on adding metadata to BIOM files](biom-format.org/documentation/adding_metadata.html).
* Updated ``qiime.util.get_qiime_library_version`` to return git commit hash rather than svn revision number (as we're using git for revision control now).
* Added java version in output of ``print_qiime_config.py`` to assist with debugging.
* Changed ``plot_rank_abundance_graph.py`` so ``-o`` specifies the filename of the figure, not the output directory anymore.
* Added new script ``add_alpha_to_mapping_file.py`` which adds alpha diversity data to a mapping file for incorporation in plots, etc.
* Moved the QIIME website files from ``Qiime/web`` to their own GitHub repository: [qiime.github.com](https://github.com/qiime/qiime.github.com).
* Fixed bug in installation of QIIME Denoiser with setup.py.
* ``supervised_learning.py`` now produces mislabeling.txt and cv_probabilities.txt that look like QIIME mapping files, allowing them to be used for coloring points in PCoA plots, etc.
* Updated RDP Classifier training code to allow any number of ranks in training files, as long as number of ranks is uniform. This removes the need for special RDP training files in reference OTU collections.
* Added table density and metadata listings to ``per_library_stats.py``.
* Updates to several dependencies. New dependencies (for those that changed in this release) are: Python 2.7.3; PyCogent 1.5.3; biom-format 1.1.1; PyNAST 1.2; usearch 5.2.236; rtax 0.983; AmpliconNoise 1.27; Greengenes OTUs 12_10; and RDP Classifier 2.2.

QIIME 1.5.0 (8 May 2012)
==================================
* OTU tables are now stored on disk in the BIOM file format (see http://biom-format.org). The BIOM format webpage describes the motivation for the switch, but briefly it will support interoperability of related tools (e.g., QIIME/MG-RAST/mothur/VAMPS), and is a more efficient representation of data/metadata. The biom-format projects DenseTable and SparseTable objects are now used to represent OTU tables in memory. See the convert_biom.py script in the biom-format project for converting between 'classic' and BIOM formatted OTU tables.
* Added a script, add_qiime_labels, that allows users to specify a directory of fasta files, along with a mapping file of SampleID<tab>fasta file name, and combines the fasta files into a single combined fasta file with QIIME compatible labels.  This is to handle situations where sequencing centers perform their own proprietary demultiplexing into separate fasta files per sample, instead of supplying raw data, but users would like to use QIIME to analyze their data.
* Added new compare_categories.py script to perform significance testing of categories/sample grouping. Added accompanying tutorial and new RExecutor class to util.py. Methods supported by compare_categories.py are Adonis, Anosim, BEST, Moran's I, MRPP, PERMANOVA, PERMDISP, and RDA. See doc/tutorials/category_comparison.rst for details. 
* compare_distance_matrices.py can now perform partial Mantel and Mantel correlogram tests in addition to the traditional Mantel test. Additionally, the script has several new options. Added new supporting tutorial and generic statistical method library code (doc/tutorials/distance_matrix_comparison.rst, qiime/stats.py, qiime/compare_distance_matrices.py), and two new classes (DistanceMatrix and MetadataMap) to util.py.
* make_3d_plots.py added a new option "-s" which by default only outputs the unscaled points, whereas user can choose to show scaled, unscaled or both.
* split_libraries_fastq.py default parameters updated based on evaluation of parameter settings on real and mock community data sets. A manuscript describing these results is currently in preparation. Briefly, the -p/--min_per_read_length parameter was modified to take a fraction of the full read length that is acceptable as the minimum, rather than an absolute (integer) length. Additionally the --max_bad_run_length default was changed from 1 to 3.
* check_id_map.py code was completely refactored to increase readability and ease of modification.  Now also creates html output to display locations of errors and warnings in the mapping file.
* Altered default value of min_length in align_seqs.py and parallel_align_seqs_pynast.py. This was previously set to 150 based on 454 FLX data, but it is now computed as 75% of the median input sequence length. This will scale better across platforms and read length, and allow for more consistent handling in of data from different sources. The user can still pass --min_length with a specific value to override the default.
* Altered the way split_libraries.py handles errors/warnings from the mapping file, and fixed a bug where suppression of warnings about variable length barcodes was not being properly passed.  Now warnings will not cause split_libraries.py to halt execution, although more serious problems (errors) will.  These includes problems with headers, SampleIDs, and invalid characters in DNA sequence fields.
* Increased allowed ambiguous bases in split_libraries.py default values from 0 to 6.  This is to accommodate the FLX+ long read technology which will often make ambiguous base calls but still have quality sequences following the ambiguous bases.  Also added an option to truncate at the first "N" character option (-x) to allow users to retain these sequences but remove ambiguous bases if desired.
* Updated merge_mapping_files.py to support merging of mapping files with overlapping sample ids.
* Added support for CASAVA 1.8.0 quality scores in split_libraries_fastq.py. This involved deprecating the --last_bad_quality_char parameter in favor of --phred_quality_threshold. The latter is now computed from the former on the basis of detecting which version of CASAVA is being used from the fastq headers (unfortunately they don't include this information in the file, but it is possible to detect).
* Added the possibility of printing the function of the curve that was fit to the points in plot_semivariogram.py 
* Replaced filter_otu_table.py with filter_otus_from_otu_table.py. The interface was redesigned, and the script was renamed for clarity.
* Replaced filter_by_metadata.py with filter_samples_from_otu_table.py. The interface was redesigned, and the script was renamed for clarity.
* Add new script to compute the coverage of a  sample (or its inverse - the conditional uncovered probability) in the script conditional_uncovered_probability.py. Current estimators include lladser_pe, lladser_ci, esty_ci and robbins.
* Updated usearch application wrapper, unit test, and documentation to handle usearch v5.2.32 as earlier version supported has bugs regarding consensus sequence generation (--consout parameter).
* Added support for the RTAX taxonomy assignment. RTAX is designed for assigning taxonomy to paired-end reads, but additionally works for single end reads. QIIME currently supports RTAX 0.981.
* Added the pick_subsampled_reference_otus_through_otu_tables.py, a more efficient open reference OTU picking workflow script for processing very large Illumina (or other) data sets. This is being used to process the Earth Microbiome Project data, so is designed to scale to tens of HiSeq runs. A new tutorial has been added that describes this process (doc/tutorials/open_reference_illumina_processing.rst).
* Added new script convert_fastqual_to_fastq.py to convert fasta/qual files to fastq. 
* Added ability to output demultiplexed fastq from split_libraries_fastq.py. 
* Added a new sort option to summarize_taxa_through_plots.py which is very useful for web-interface. By default, sorting is turned off.
* Added ability to output OTUs per sample instead of sequences per sample to per_library_stats.py.
* Updates and expansions to existing tutorials, including the using AWS and procrustes analysis tutorials. 
* Added insert_seqs_into_tree.py to insert reads into an existing tree. This script wraps RAxML, ParsInsert, and PPlacer. 
* Updated split_libraries_fastq.py to handle look only at the first n bases of the barcode reads, where n is automatically determined as the length of the barcodes in the mapping file. This feature is only use if all of the barcodes are the same length. It allows qiime to easily handle ignoring of a 13th base call in the barcode files - this is a technical artifact that sometimes arises.
* Added new stats.py module that provides an API for running biogeographical statistical methods, as well as a framework for creating new method implementations in the future (this code was moved over from qiimeutils/microbiogeo). Also added two new classes to the util module (DistanceMatrix and MetadataMap) that are used by the stats module.
* Updated Mothur OTU picker support from 1.6.0 to the latest (1.25.0) version.
* Added start_parallel_jobs_sg.py to support parallel jobs on SGE queueing systems.
* Modified split_libraries_fastq.py and format.py to show SampleIDs with zero sequence count and to show the total sum of sequences written in the log file.

QIIME 1.4.0 (13 Dec 2011)
==================================
* Implemented usearch (ie OTUPIPE) as chimera detection/quality filtering/OTU picking in the pick_otus.py module.
* All workflows now log the md5 sums of all input files (trac #92).
* Testing of QIIME with new dependency versions, updating of warnings and test failures (in print_qiime_config.py). No code changes were required to support new versions.
* split_libraries_fastq.py can now handle gzipped input files.
* Addition of code and tutorial to support plotting of raw distance data in QIIME (scripts/make_distance_comparison_plots.py, scripts/make_distance_boxplots.py, qiime/group.py, doc/tutorials/creating_distance_comparison_plots.rst).
* Updates to many scripts to support PyCogent custom option types (new_filepath, new_dirpath, etc.). 
* Fixes to workflows to fail immediately on certain types of bad inputs (e.g., missing tree when building UniFrac plots) rather than failing only when the script reaches the relevant step in the workflow.
* Added ability to merge otu tables with overlapping sample ids (in merge_otu_tables.py). Values are summed when an OTU shows up in the same sample in different OTU tables.
* Added a new script (filter_distance_matrix.py) to filter samples directly from distance matrices.
* Added script nmds.py Non-Metric Multidimensional Scaling (NMDS).
* Added in the calculation of standard error in rarefaction plots, since only standard deviation was calculated. Also added an optional option choice for this.
* Support for pick_otus_through_otu_table.py to allow for uclust_ref to be run in parallel with creation of new clusters. 
* Added script distance_matrix_from_mapping.py which allows to create a distance matrix from a metadata column.
* assign_taxonomy_reference_seqs_fp and assign_taxonomy_id_to_taxonomy_fp were added to qiime_config, which allows users to set defaults for the dataset they'd like to perform taxonomy assignment against. This works for the serial and parallel versions of assign_taxonomy for both BLAST and RDP.
* Added in make_3d_plots.py the possibility of calculating RMS vectors, using two methods: avg and trajectory, to assess power (movement) of the trajectories. Additionally this feature will return the significance of the difference of the trajectories using ANOVA.
* Added in make_3d_plots.py the possibility of adding vectors or traces of individuals in space; this can be helpful in time series analysis.
* Added additional allowed characters to data fields in mapping files.  These include space and /:,; characters.  All characters allowed now are: alphanumeric, underscore, space, and +-%./:,; characters.
* split_otu_table.py now can keep duplicated rows in the resulting mapping files and can rename sample names (SampleID), both in the resulting mapping files and the otu tables, with other column of the mapping file; this can be helpful for Procrustes analysis.
* plot_semivariogram.py now lets you control colors and axis of the resulting plots, and ignore missing samples, this can be useful when samples are missing after rarefying.
* default num_dimensions for transform_coordinate_matrices.py changed from all dimensions to 3 (trac ticket #119). This more closely corresponds with how we use this test (e.g., to determine if we would draw the same biological conclusions from two different methods of generating a PCoA plot). This was in response to our noticing that monte carlo p-values were lower than we would expect in controls.
* Removed the --suppress_distance_histograms option from beta_diversity_through_plots.py in favor using the -c/--histogram_categories option to determine whether these will be generated. If the user passes -c, distance histograms are generated. If they do not, these are not generated.
* Added support for fastq files in count_seqs.py.
* Several new tutorials including retraining of the RDP classifier, working with Amazon Web Services, basic unix/linux commands, and others.
* Fixed bug in process_qseq.py that would result in only a single input file per lane have it's data stored in the fastq.
* Fixed bug in filter_otu_table where sampleIDs would remain despite all OTU counts being zero.
* Fixed bug in serial pick_reference_otus_through_otu_table.py that was causing uclust to be used rather than uclust_ref as the default method for otu picking.
* Added option to support reverse complements of golay barcodes in the mapping file.
* Modifed beta_diversity_through_plots.py so distance histograms are only generated if the user specifies --histogram_categories on the command line. These are very slow to generate for all mapping categories, so it makes more sense for the user to turn on histogram plotting for the specific categories they're interested in.
* Added option, --reverse_primer_mismatches to split_libraries.py to allow setting of distinct mismatches from forward primer.
* Added option (-e/--max_rare_depth) to the command line of alpha_rarefaction.py. This allows for a convenient way for users to specify the maximum rarefaction depth on the command line, and is useful for when it needs to be set to something other than the median rarefaction depth. Also added option to control minimum rarefaction depth from the alpha_rarefaction.py command line.
* Added support for 5- and 10-fold and leave-one-out cross-validation to supervised learning.
* Added filter_by_metadata.py state string handling to filter_fasta.py for metadata-based fasta filtering.
* Added subsample_fasta module for randomly subsampling fasta files.
* Added script to split a post-split-lib fasta file into per-sample fasta files. This is useful for sharing Illumina data with collaborators or creating per-sample files for DB submission.
* Fixed bug where multiple_rarefactions_even_depth didn't work with --lineages_included.
* Modified pick_otus_through_otu_table.py so filter_alignment.py can be applied when the method is other than PyNAST. This previously wasn't possible because we only filtered with the lanemask, but we now allow entropy filtering, so this is relevant.
* Fixed two serious bugs in make_distance_histograms.py related to p-value calculations (both Monte Carlo and parametric p-values were affected).
* Removed several obsolete scripts (make_pie_charts.py and several denoiser-related scripts).
* Added muscle_max_memory option to align_seqs script.
* Changed default num_dimensions to 3 in transform_coordinate_matrices.py. This more closely corresponds with how we use this test (e.g., to determine if we would draw the same biological conclusions from two different methods of generating a PCoA plot). This was in response to our noticing that monte carlo p-values were lower than we would expect in controls.

QIIME 1.3.0 (29 June 2011)
==================================
* uclust and uclust_ref OTU pickers now incorporate a pre-filtering step where identical sequences are collapsed before calling uclust and then expanded after calling uclust. This gives a big speed improvement (5-20x) on reasonably sized input sets (>200k sequences) with no effect on the resulting OTUs. This is now the default behavior for pick_otus.py, and can be disabled by passing --suppress_uclust_prefilter_exact_match to pick_otus.py.
* Added ability to pass a file to sort_otu_table.py that contains a sorted list of sample ids, and use that information rather than the mapping file for sorting the OTU table. This allows users to, e.g., pass sorted mapping files as input. 
* Added core_analyses.py script and workflow function. This plugs together many components of QIIME (split libraries, pick_otus_through_otu_table.py, beta_diversity_through_3d_plots.py, alpha_rarefaction.py) into a single command and parameters file.
* Added script (split_otu_table_by_taxonomy.py) which will create taxon-specific OTU tables from a master OTU table for taxon-specific analyses of alpha/beta diversity, etc.
* Changed default behavior of single_rarefaction.py. Now lineage information is included by default, but can be turned off with --suppress_include_lineages
* Added script (compare_distance_matrices.py) for computing mantel correlations between a set of distance matrices.
* Interface changes to summarize_otu_by_cat.py. This allows the user to pass the output file name, rather than a directory where the output file should be written.
* Parameter -r reassignment in parallel_assign_taxonomy_rdp.py. Now -r is used for reference_seqs_fp as before was for rdp_classifier_fp.
* Added script inflate_denoiser_output.py to expand clusters to fasta representing all sequences. This allows denoiser results to be passed directly to the OTU pickers (and OTU picking workflows) which should greatly reduce the complexity of denoiser runs. The "Denoising 454 Data" tutorial has been updated to reflect how the pipeline should now be run. The denoising functionality was removed from the pick_otus_through_otu_table.py workflow script as that could only be used in very special circumstances - this allows us to focus our attention on supporting the new pipeline described in the updated tutorial.
* Reorganized output from pick_otus_through_otu_table.py to get rid of the confusing output directory structure. 
* Added script plot_semivariogram.py to plot semivariograms using two distance matrices. This script also plots a fitting curve of the data values.
* Changed beta diversity scripts to do unweighted_unifrac,weighted_unifrac by default.
* Changed output of summarize_taxa.py to a directory instead of filepath. This allows for multiple levels to be processed simultaneously.
* The beta_diversity_through_3d_plots.py now contains some additional functionality -- 2d plots and distance histograms. It has therefore been renamed beta_diversity_through_plots.py. Any of the plots can be disabled by passing the options  --suppress_distance_histograms, --suppress_2d_plots, and --suppress_3d_plots.
* Updated required version of FastTree to 2.1.3 as this version contains some bug fixes over version 2.1.0.
* Modified single_rarefaction.py so default is to include lineages (previously did not include these by default).
* Added split_otu_table.py script which splits a single OTU table into several OTU tables based on the values in a specified column of the mapping file. This is useful, for example, when a single OTU table is generated that covers multiple studies.
* Fixed bug in mouseovers in taxa area and bar charts. These were misaligned when a lot of samples were included.
* Added support for RDP classifier 2.2. Versions 2.0 and 2.2 are both supported.
* Added support for AmpliconNoise with the ampliconnoise.py script.
* Added new page to the documentation to cover upgrades between versions of QIIME.
* Updated the make_distance_histograms.py output filepaths and HTML layout to be more consistent with other plotting scripts.
* Added a new taxonomy summary workflow (summarize_taxa_through_plots.py).
* Modified workflow scripts so stdout and stderr are written to the log file. This is very useful for debugging.
* Added new script (simsam.py) to simulate samples using a phylogentic tree.
* Complete overhaul of Illumina data processing code. QIIME now treats fastq format as the default for Illumina data, and various other formats can be converted to fastq using process_qseq.py and process_iseq.py. The "Processing Illumina Data tutorial" has also been completely overhauled and describes these changes. The primary script for demultiplexing Illumina data is now split_libraries_fastq.py.
* Dropped support for PyroNoise in favor of AmpliconNoise (the successor to PyroNoise) and the QIIME denoiser.
* Added inflate_denoiser_output.py script to simply the integration of denoiser results into the QIIME pipeline. See the "Denoising 454 Data" tutorial, which has been overhauled in this release. To reduce the possible pathways through QIIME with denoising, support for denoising was removed from pick_otus_through_otu_table.py in favor of working with the pipeline presented in the tutorial.
* Changed default behavior of split_libraries.py so unassigned reads are not stored by default. There is now a --retain_unassigned_reads option to achieve the previous behavior.
* Many clean-ups to the script documentations through-out QIIME.
* Adding scripts to plot semivariograms.
* Modified all workflow scripts so parameter files are now optional. This will simplify working with 'default' analyses in these scripts.
* Added more thorough support for floating point values in OTU tables. This was previously supported only in specific cases.
* Added support for users to pass jobs_to_start on the command line for all of the workflow scripts. This overrides this value in the parameters file and qiime_config, and is a more convenient way of controlling this.
* Added entropy filtering option to filter_alignment.py. This can be useful for position-filtering de novo alignments, or other alignments where no lanemask is available.
* Added new script (count_seqs.py) which will count the number of sequences in one or more fasta file, as well as the mean/stddev sequence lengths, and print the results to stdout or file.
* Added the plot_taxa_summary.py workflow script, which includes summarizing the OTU table by category.
* Overhauled the QIIME overview tutorial. 
* Added new script (start_parallel_jobs_torque.py) which can be used for running parallel QIIME on clusters using torque for the queueing system. A new qiime_config value, torque_queue, can be specified to define the default queue.
* Integrated the QIIME Denoiser (Reeder and Knight, 2011) into Qiime.
* Added script (compare_alpha_diversity.py) for comparing rarefied alpha diversities across different mapping file categories.
* Fixed bug in pick_otus.py where reverse strand matching did not work for uclust/uclust_ref.
* Modified location where temp files are written for more consistency through-out QIIME. Temp files are now written the temp_dir (from qiime_config) or /tmp/ if temp_dir is not defined. There may still be a few temp files being written to other locations, but the goal is that all will write to the same user-defined (or default) directory.
* Added split_otu_table.py script which splits a single OTU table into several OTU tables based on the values in a specified column of the mapping file. This is useful, for example, when a single OTU table is generated that covers multiple studies.
* Added script (make_tep.py) that makes TopiaryExplorer project file (.tep) from an otu table, sample metadata table and tree file.
* Removed the rdp_classifier_fp from qiime_config. This was used inconsistently through-out QIIME, so was somewhat buggy, and with the switch to RDP 2.2 in QIIME 1.3.0 I think it will save a lot of support headaches to just get rid of it.
* Added tutorial for processing 18S data, along with a small 3 domain sample sequence file in the qiime_tutorial/18S_tutorial_files/ folder.
* Added filter_tree.py script, which functions similarly to filter_fasta.py. Moved some functions from filter_fasta.py to filter_tree.py that were generally useful.

QIIME 1.2.1 (22 Feb 2011)
==================================
* Added submit_to_mgrast.py script which takes a post-split-libraries fasta file and submits it to the MG-RAST database.
* Added sort_otu_table.py script which allows for sorting samples in an OTU table based on their associated values in a mapping file.
* Remove DOTUR OTU picker. This was requested by Pat Schloss as Mothur has replaced DOTUR.
* Removed support of SRA submission and processing scripts along with related documentation and tutorial. This included the following scripts: make_sra_submission, sra_spreadsheets_to_map_files, process_sra_submission (starting revision 1786).
* Added categorized_dist_scatterplot.py script.
* Added OTU gain as a new beta diversity metric to compute non-phylogenetic gain (G).
* Added features to split_libraries to allow truncation or removal of sequences with quality score windows, and increased information deposited in log file about sliding window quality score tests.  Added unit test for quality score truncation/removal.
* Added reference-based OTU picking workflow script. This can be applied for database OTU picking, as well as for applying Shotgun UniFrac (Caporaso et al. 2011, PLoS One, accepted).
* Added a new list of distinct colors to the colors.py module
* Added Area and Bar taxa summary plots to a new script plot_taxa_summary.py.  This script allows for writing of Pie Charts as well, thereby deprecating the make_pie_charts.py script.
* Added support for output of biplot coords to make_3d_plots script (SF feature req. 3124713).
* --stable_sort option enabled by default for uclust OTU pickers.
* Changed defaults for uclust and uclust_ref OTU pickers. The new parameters make both OTU pickers about 2-3x slower, but the resulting clusters are significantly better in terms of making the best choice of OTU for a given sequence, and ensuring that cluster seeds are less than 97% identical to one another. The default rep seq picking method was also changed to "first" from "most_abundant" which ensures that the seed sequence is chosen as the representative for a cluster. Abundance is instead taken into account at the otu picking stage (as it has been for a while) by pre-sorting the sequences by abundance so most abundant sequences are more likely to be seeds. In practice, with presorting by abundance, the same sequence is usually chosen as the representative when passing first or most_abundant as the OTU picking method.
* Added support for generating inVUE plots in make_3d_plots.py.
* Changed tree type default for upgma comparisons, to consensus tree rather than the upgma tree based on the full otu table.
* Disabled the check that jobs_to_start > 1 in a user's qiime_config before allowing them to start parallel jobs. This is inconvenient in several places (e.g., EC2 images when used with n3phele), and after some discussion we decided that it should be up to the user to have understood how parallel qiime should be configured before using it.
* Added ability to pool primers for mapping files passed to check_id_map and split_libraries.py.  Primers are separated by commas, and autodetected.
* Added sort_otu_table.txt for sorting the sample IDs in an OTU table based on their value in a mapping file.
* Changed the method for p-value calculation in Procrustes analysis Monte Carlo in response to SF bug # 3189200.

QIIME 1.2.0 (10 Nov 2010)
==================================
* When computing jackknife support for sample clustering (e.g.: UPGMA sample trees), Qiime can now compute a consensus tree from the jackknife replicates, in addition to the existing functionality of using the full dataset as the master tree, and annotating that tree with jackknife support values. See jackknifed_beta_diversity.py --master_tree and consensus_tree.py .
* Added the ability to write out the flowgram file in process_sff.py, ability to define an output directory and convert Titanium reads to FLX length.
* SRA submission protocol updated to perform human screening with uclust_ref against 16S reference sequences, rather than cdhit/blast against reference sequences. This can be a lot faster, and reduces the complexity of the code by requiring users to have uclust installed for the human screen rather than cdhit and blast.
* Updated SRA protocol to allow users to skip the human screening step as this takes about 2/3 or more of the total analysis time, and is not relevant for non-human-derived samples (e.g., soil samples).
* Added ability to pass --max_accepts, --max_rejects, and --stable_sort through the uclust otu pickers. 
* Added a -r parameters to pick_rep_set.py to allow users to pass "preferred" representative sequences in a fasta file. This is useful, for example, if users have picked OTUs with uclust_ref, and would like to use the reference sequences as their representatives, rather than sequences from their sequencing run. 
* Renamed Qiime/scripts/jackknifed_upgma.py to Qiime/scripts/jackknifed_beta_diversity.py to reflect the addition of generating jackknifed 2d and 3d plots to this workflow script.
* Updated parallel_multiple_rarefactions.py, parallel_alpha_diversity.py, and parallel_beta_diversity.py to use the jobs_to_start value for better control over the number of parallel runs.
* uclust_ref otu picker now outputs an additional failures file listing the sequences which failed to cluster if the user passed --suppress_new_clusters. This is done for ease of parsing in downstream applications which want to do something special with these sequences. The failures list is no longer written to the log file (although the failures count is still written to the log file).
* Added the filter_fasta.py script which allows users to build a fasta file from an existing fasta where specified sequences are either included or excluded from the new file. The sequences to keep or exclude can be specified by a variety of different inputs, for example as a list of sequence identifiers in a text file.
* Added parallel version of uclust_ref OTU picker.
* Added negative screen option to process_sra_submission.py -- this allows users to screen by discarding all sequences that match a reference set, while the (default) positive screen allows users to screen by retaining only sequences that match a reference set.
* Added options to split_libraries.py to enable the detection and removal of reverse primers from input sequences, and an option to record a filtered quality score output file that matches the bases found in the output seqs.fna file.
* Added the trflp_file_to_otu_table.py script that allows users to create an OTU table simile from a Terminal restriction fragment length polymorphism (T-RFLP) text file.
* Added min_aligned_length parameter to the BLAST OTU picker. By default, BLAST alignments now must cover at least 50% of the input sequence for OTU assignment to occur.
* Changed default randomization strategy in Procrustes monte carlo from shuffling within coordinate vectors to shuffing the labels on the vectors themselves. This doesn't appear to affect clearly significant cases at all, but is more conservative and therefore favors non-significance of results in borderline cases.
* Added ability to run beta diversity calculations in parallel at the single OTU table level to improve performance when computing diversity on very large collections of samples. This functionality is now hooked up to the beta_diversity_through_3d_plots.py workflow script, and includes the new -r parameter to beta_diversity.py which allows users to specify samples to compute diversity vectors for (rather than requiring that the full all-against-all diversity matrix is created).
* uclust-based analyses now retain the .uc files as these contain a lot of useful information that was previously being discarded.
* Improved handling of blank lines in parse_otu_table -- these are now ignored. Other improvements were made to the parse_otu_table format to better support these files coming from sources other than QIIME (such as MG-RAST).
* Allow the -R option to be passed to ChimeraSlayer. Closes feature request 3007445.
* Added capability for pairwise sample/sample, monte carlo significance tests. These are frequently done via the unifrac web interface. Users hitting max size limitations on the web can now thrash their own hardware.
* Fixed a bug in make_rarefaction_plots where the table below the plots had column labels sorted by natsort, while the values in the table were sorted arbitrarily by dict keys. The plots themselves were fine.
* Added a Procrustes analysis/plotting tutorial.
* Added code to exclude OTU ids from an OTU table when building the OTU table. This allows users to discard OTUs that were identified as chimeric. Accessible by passing --exclude_otus_fp to make_otu_table.py. 
* Modified identify_chimeric_sequences.py to no longer require the ref db in unaligned format when using chimeraSlayer.
* Added a tutorial document on applying chimera checking in QIIME.
* Added ability to pass -F T/F to parallel_blast to allow disabling of the low-complexity filtering in BLAST.
* Added new script (shared_phylotypes.py) for computing shared OTUs between pairs of samples. Batch mode can be used in combination with dissimilarity_mtx_stats.py to calculate stats for a set or rarefied OTU tables.
* Added min_aligned_percent parameter to BLAST OTU picker workflow, with default set at 50%. This will now require that an alignment must cover at least 50% of a sequence OTU assignment to occur.
* Add script to draw rank abundance graphs (plot_rank_abundance_graph.py).
* Modified interface of make_distance_histograms so --html_output is now the default. A new parameter, --suppress_html_output, was added to produce the old behavior.
* Added script (quality_scores_plot.py) to plot quality score by position given a .qual file. This is useful with another new script (truncate_fasta_qual_files.py) to truncate fasta/qual files at the point where quality begins to decrease, and has been useful in controlling for quality issues on 454 Ti runs.
* Added binary SFF parsing module from PyCogent, removed sfftools dependency from workflow test, process_sff, and other areas of QIIME.
* Added ACE calculation to alpha_diversity.py.
* Updated documentation on file formats used by Qiime.
* Added more extensive error checking in parse_mapping_file to handle some cryptic error messages that were arising from scripts that were passed bad mapping files.
* Added capability to perform supervised classification of metadata categories using the Random Forests classifier. Outputs include a ranking of OTUs by discriminatory power, and the estimated probability of each metadata category for each sample. The latter may be useful for detecting potentially mislabeled samples.


QIIME 1.1.0 (14 May 2010)
=========================
* Additional field added to BLAST assign taxonomy output to indicate the best BLAST hit of the query sequence -- this is in response to Sourceforge feature request 2988407.
* Added presorting by abundance to uclust OTU picker. The idea here is that sequences which are more abundant are better representatives when clustering, so they should come first in the file. Also added ability to pass the optimal flag to uclust, which should also improve uclust-picked OTUs, which comes with a performance hit.
* Added Confidence interval display (jackknifed pcoa) in make_2d_plots and make_3d_plots. After performing multiple_rarefactions, beta_diversity and principal_coordinates on an OTU table, the user can supply the resulting directory to both of these scripts.  Currently the user has the option of performing InterQuartile Range (IQR) or standard-deviation (sdev) on the principal coordinate files and ellipses are drawn around each point to represent the confidence interval in each P.C.  Along with this option, the user can manipulate the opacity of the ellipses as well.
* Updated the display for rarefaction plots, so the legend does not overlap with the plots and fixed the display of the rarefaction average table in the webpage.  Now the user can switch between plots with different metrics and categories by using the drop down menus.  The user can also display the samples that contribute to the average for that group.  Below the plots, a table is displayed to show the rarefaction average data with all the distance metric values.
* Merged the make_rarefaction_averages into the make_rarefaction_plots script.  Also removed the inputs (--rarefaction_ave and --ymax) options, since they are determined by the script.  Also, restructured the output directory format and combined all metric data into one html.
* Added the uclust_ref OTU picker, which uses uclust to pick OTUs against a reference collection. Sequences which are within the similarity threshold to a reference sequnece will cluster to an OTU defined by that reference sequence, and sequences which are outside of the similarity threshold to any reference sequence will form new OTUs.
* The interface for exclude_seqs_by_blast.py has changed.  -M and -W options are now lowercase to avoid conflicts with parallel scripts.  Users can avoid formatting the database by passing --no_format_db.  By default the files created by formatdb are now cleaned up. Users can choose not to  clean up these files  using the --no_clean option.  Output file extensions have changed from ".excluded" to ".matching" and from ".screened" to ".non-matching" to be clear regardless of whether the sequences matching the database, or not matching the database, are to be excluded. A check was added for user-supplied BLAST databases in exclude_seqs_by_blast.py when run with --no_format_db: if the required files do not exist a parser error is thrown
* Added ability to chimera check sequences with ChimeraSlayer. See identify_chimeric_seqs.py for details. 
* Added workflow script for second-stage SRA submission, process_sra_submission.py. The SRA submission tutorial has been extensively updated to reflect the use of this new script.
* Added the ability to supply a tree and sort the heatmap based on the supplied tree.
* Added the ability to handle variable length barcodes, variable length primers, and no primers with split_libraries.py. Error-correction is not supported for barcode types other than golay_12 and hamming_8. split_libraries.py also now throws an error if the barcode length passed on the commands line does not match the barcode length in the mapping file.
* Updated the print_qiime_config.py script to print useful debugging information about the QIIME environment. 
* Added high-level logging functionality to the workflow scripts.
* Added RUN_ALIAS field to SRA experiment.txt spreadsheet in make_sra_submission.xml.



QIIME 1.0.0 - (8 Apr 2010)
===========================
* uclust made default OTU picker (instead of cdhit).
* uclust made default pairwise aligner for PyNAST (instead of BLAST).
* Minimum PyNAST version requirement upgraded to PyNAST 1.1.
* Minimum PyCogent version requirement upgraded to PyCogent 1.4.1.
* tree_compare now can compare trees where some tips aren't present in all trees.
* --small_included option removed from rarefaction scripts.
* Added "remove outliers" functionality to filter_alignment.py.  After removing lanemasked columns and gap columns, -r will remove outlying sequences, preventing odd spikes in phylo trees when some seqs are poorly aligned.
* Absent samples are now included in the output of unifrac like metrics - 0 dist between two samples that aren't there, 1 dist between an absent and a present sample.
* make phylogeny now does good midpoint rooting (still off by default).
* Consolidated parsing functionality to qiime.parse.
* Removed dependence on several qiime_config values - users should run Qiime/scripts/print_qiime_config.py -t to get information on parameter settings which are outdated.
* Added an example 'cluster_jobs' -- start_parallel_jobs.py -- script which will give users in multi-core or multi-proc environments very easy access to parallel QIIME. This also adds parallel support to the QIIME virtual box.
* Modified the default value of jobs_to_start to be 1 -- because of the addition of the example cluster_jobs script, the default value of 24 no longer makes sense (if it ever really did...). Because the new script is built for multi-core/multi-proc environments, 24 is too high for most cases. Users will need to modify this value from 1 (corresponding to no parallelization) to a value that makes sense for their environment (e.g., 2 for dual core, or 24 to get the previous default).
* Added colors module and tests to consolidate and standardize coloring code in QIIME - also updated the graphics scripts to use the colors module.
* Added ability for user to specify the background colors of plots in prefs files or on the command line.
* Tweaked SRA submission routines in accordance with accepted format from JCVI's 
survey of multiple body sites.
* Fixed SF bug #2971581, which was an issue with the path to qiime's scripts directory not being determined correctly when qiime was installed using setup.py. qiime_config now contains a key (empty by defualt) for the qiime_scripts_dir. If this is not specified by the user, it is determined from the qiime project dir.
* Renamed scripts/make_3d_prefs_file.py as scripts/make_prefs_file.py to reflect that the prefs files are now used by other scripts.
* Changed behavior of color-by option to make_3d_plots, make_2d_plots, and make_rarefaction_plots, so if no -b option or prefs files is provided, scripts default to coloring by all values. Consequently, mapping files are also now required for these scripts.
* Added a split_libraries_illumina.py script to handle processing of Illumina GAIIx data.
* Added an additional rarefaction script for clarity. There are now 3 scripts to handle rarefaction: single_rarefaction takes one input otu table into one output table, allows manual naming,  multiple_rarefactions makes auto-named rarefied otu tables at a range of depths, and multiple_rarefactios_even_depth.py makes auto-named tables all at the same depth.
* Added workflow unit tests (with timeout functionality). 
* Added default alpha and beta diversity metrics to qiime_parameters.txt.
* Integrated Denoiser (Jens Reeder's 454 denoiser) wrappers, and tied this into the workflow scripts.
* Added biplot functionality.  make_3d_plots now takes the -t option (off by default) to include taxa on the pcoa plot.
* Updated the QIIME tutorial to use the workflow scripts where possible. Additionally added the tutorial data set in the svn repository.
* Reorganization and expansion of the documentation through-out.
* Added sanity checks to print_qiime_config.py. This will now allow users to evaluate their environment, and should help with debugging.
* Added new field to qiime_config (temp_dir) which will be used to specify where temp files should be written. Currently this is only used by the workflow tests, and is intended to allow users to specify something other than /tmp for cases when /tmp is not shared between all nodes that might be working on a job. This will eventually be used for all temp dir creation.
* Added ability to make summary plots for a directory of coordinate files in make_3d_plots and make_2d_plots. The summary plot adds ellipsoidal confidence intervals around each point in the plot.
	


QIIME 0.92 - (3 Mar 2010)
=======================
* Removed outdated documentation PDFs, along with references to those PDFs in the README and INSTALL documents.

QIIME 0.91 - (3 Mar 2010)
=======================

* Addition of a uclust-based OTU picker.
* Transfer of all command line interfaces from Qiime/qiime to Qiime/scripts -- this was an important change as it allowed us to get away from the previously one-to-one relationship between files in our library code (in Qiime/qiime) and the command line interfaces.
* Standardized command line interfaces for all code in Qiime/scripts by using a new function, Qiime.qiime.util.parse_command_line_parameters to handle the command line interfaces.
* Moved to Sphinx for documentation, and developed a framework for extracting script documentation directly from the scripts to populate the web documentation. 
* Bug fixes through-out the code base, including but not limited to fixes for Sourceforge tickets: 2957503, 2953765, 2945548, 2942443, 2941925, 2941926, 2941717, 2941396, 2939588, 2939575, 2935939.
* Updated the all_tests.py script to perform a minimal test of the scripts (getting help text works as expected), and to alert users if unit tests may be failing due to missing external applications, in which case they may not be critical.
* Created a directory for pycogent_backports, where we can temporarily store new code that has been added to PyCogent, but which has not been added to a PyCogent release yet. This will allow us to keep QIIME's dependencies on the latest PyCogent version despite rapid and frequently related changes in both packages.
* Added code for performing Procrustes analyses of coordinate matrices, and graphing the results of those analyses in 3d plots (see transform_coordinate_matrices.py and compare_3d_plots.py).
* Performance enhancements related to golay barcode decoding.
* Added setup.py to help with installation of QIIME - this will put the library code in site-packages, and the scripts in /usr/local/bin (both locations can be changed via command line options to setup.py).
* Created a support_files directory to hold jar, js, png, and other required files.
* Added Pearson correlation to list of options in otu_category_significance.py.
* Workflow scripts added for running large repetitive processes with a single command rather than multiple commands -- in scripts, see beta_diversity_through_3d_plots.py, pick_otus_through_otu_table.py, alpha_rarefaction.py, jackknifed_upgma.py.



QIIME 0.9 - (25 Jan 2010)
=======================

* Initial release
