.. _taxa_summary_comparison:

========================
Comparing Taxa Summaries
========================

Introduction
------------
This tutorial explains how to use a few different taxa summary comparison techniques that are available in `compare_taxa_summaries.py <../scripts/compare_taxa_summaries.html>`_. Taxa summary files are compared using either Pearson or Spearman correlation coefficients, which are computed between pairs of samples. The script provides two different comparison modes. `Paired` comparison mode compares samples that match between the two taxa summary files, while `expected` comparison mode compares all samples in one taxa summary file to an "expected" sample in the other file.

One potential use for this script is to compare two taxa summary files that have been constructed using different taxonomic classification techniques. For example, you might want to see whether the taxonomic composition of samples assigned using the RDP classifier are correlated with the taxonomic composition of samples derived from BLAST or RTAX. Another example is to see whether using different confidence levels or E-values (for RDP and BLAST, respectively) yield taxonomic compositions that are similar/correlated.

Another example use-case is to compare the taxonomic composition of several mock community replicate samples to a single expected, or known, sample community.

This script is also useful for sorting and filling taxa summary files so that each sample has the same taxa listed in the same order (with missing taxa reporting an abundance of zero). The sorted and filled taxa summary files can then be passed to a script, such as `plot_taxa_summary.py <../scripts/plot_taxa_summary.html>`_, to visually compare the differences using the same taxa coloring scheme.

Please note that this tutorial does not attempt to cover every possible option that can be used in the taxa summary comparison script. Instead, it attempts to provide useful examples to give you an idea of how to use this script in your own study, as well as customize some of the output to your liking. For a complete listing of the available options, please refer to the `compare_taxa_summaries.py <../scripts/compare_taxa_summaries.html>`_ script documentation.

Input Files
-----------
You can obtain the files used in this tutorial `here <https://s3.amazonaws.com/s3-qiime_tutorial_files/taxa_summary_comparison_tutorial_data.zip>`_. The files are derived from the `overview tutorial <./tutorial.html>`_. The unzipped directory will contain two taxa summary files named :file:`ts_rdp_0.60.txt` and :file:`ts_rdp_0.80.txt`. These taxa summaries were constructed from the overview tutorial representative sequence set by assigning taxonomy to it using the RDP classifier with two different confidence levels (0.60 and 0.80). The gg_otus_4feb2011 97% representative set was used as training input to the RDP classifier.

Output Files
------------
Depending on the options that you invoke on the script, either three or four output files will be generated. They will all be placed in the directory specified by the required -o option. Most of the output files will be tab-separated text files containing information about the test that was performed and its results. These can easily be viewed in a spreadsheet program such as Excel.

Two of the output files will be the sorted and filled versions of the input taxa summary files, which can then be used in `plot_taxa_summary.py <../scripts/plot_taxa_summary.html>`_ to visualize the differences in taxonomic composition. These files will be named based on the basename of the input taxa summary files. If the input files' basenames are the same, the output files will have '0' and '1' appended to their names to keep the filenames unique. The first input taxa summary file will have '0' in its filename and the second input taxa summary file will have '1' in its filename.

The third output file will contain the results of the overall comparison of the input taxa summary files using the specified comparison mode. The file will be named :file:`overall_comparison.txt` and will be in tab-delimited format for easy viewing in a program such as Excel.

If --perform_detailed_comparisons is specified, a fourth output file named :file:`detailed_comparisons.txt` will be created. This file will also be in tab-delimited format and will contain the results of comparing individual pairs of samples (whereas :file:`overall_comparison.txt` contains results from comparing all pairs of samples at once).

Paired Comparison
-----------------
In order to see whether the results of using the RDP classifier at a confidence
level of 0.60 versus 0.80 are correlated, we need to compute the correlation
between each of the paired samples in the two taxa summary files. To do this,
run the following command: ::

    compare_taxa_summaries.py -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp

This command will create a new output directory named :file:`taxa_comp`, which will contain three files named :file:`ts_rdp_0.60_sorted_and_filled.txt`, :file:`ts_rdp_0.80_sorted_and_filled.txt`, and :file:`overall_comparison.txt`. The first two files are sorted and filled versions of the input taxa summary files, meaning that the taxa (rows) have been sorted and any missing taxa have been added with zero abundance. This results in taxa summary files that have the same taxa listed in the same order. Open up :file:`overall_comparison.txt` to see the results of the test:

.. note::

    * # Correlation coefficient: pearson. Performed a two-tailed test of significance using a t-distribution.
    * # Number of samples that matched between the taxa summary files: 9
    * Correlation coefficient p-value
    * 0.9993  0.0000

The comment at the top of the file tells us that Pearson's correlation coefficient was used to compute the correlation between samples. The value of 0.9993 is very close to +1, which indicates strong positive correlation between the paired samples in the two input taxa summary files. The p-value indicates that the correlation is statistically significant. Thus, it appears that using the RDP classifier at a confidence level of 0.60 or 0.80 yields highly correlated results (i.e. you might draw the same biological conclusions from using either taxa summary file).

Note the second line of the file, which tells us that 9 samples matched between the two taxa summary files. All 9 samples matched between the two files because both files have the same sample IDs. The script determines which samples to compute the correlation between based on matching sample IDs between the two files (the order of sample IDs does not matter between the two files, but the names must match exactly). Any sample IDs that do not match will not be considered when calculating the results.

If you need to compare taxa summaries that do not have matching sample IDs, you
can use a sample ID map to provide a mapping between sample IDs. To illustrate
this, run the following command: ::

    compare_taxa_summaries.py -i ts_rdp_0.80.txt,ts_rdp_0.60_renamed.txt -m paired -o taxa_comp_using_sample_id_map -s sample_id_map.txt -c spearman

The second input taxa summary file (:file:`ts_rdp_0.60_renamed.txt`) is simply the original :file:`ts_rdp_0.60.txt` file with all sample IDs starting with 'PC.' renamed to 'S.'. The sample ID map (:file:`sample_id_map.txt`) specifies how to pair up samples between the two taxa summary files now that the sample IDs do not exactly match.

This command will create a new output directory named :file:`taxa_comp_using_sample_id_map`, which will contain the three files like before. Open up :file:`overall_comparison.txt` to see the results of the test:

.. note::

    * # Correlation coefficient: spearman. Performed a two-tailed test of significance using a t-distribution.
    * # Number of samples that matched between the taxa summary files: 9
    * Correlation coefficient p-value
    * 0.9687  0.0000

Notice that all 9 samples were still included because we provided a sample ID map. If we had not provided a sample ID map, the script would have raised an error stating that the two taxa summary files were incompatible because no matches could be found.

Also notice that we used the -c option to use the Spearman correlation coefficient instead of the Pearson correlation coefficient. This is why we obtained a different value for the correlation coefficient.

If you'd like to see the correlation tests applied to each pair of samples (instead of to all pairs of samples at one time), use the --perform_detailed_comparisons option: ::

    compare_taxa_summaries.py -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp_detailed --perform_detailed_comparisons

In addition to the three previous file, you will find the file :file:`detailed_comparisons.txt` in the newly created directory :file:`taxa_comp_detailed`: ::

    # Correlation coefficient: pearson. Performed a two-tailed test of significance using a t-distribution.
    # Number of samples that matched between the taxa summary files: 9
    Sample ID	Sample ID	Correlation coefficient	p-value	p-value (Bonferroni-corrected)
    PC.354	PC.354	0.9997	0.0000	0.0000
    PC.355	PC.355	0.9998	0.0000	0.0000
    PC.356	PC.356	0.9999	0.0000	0.0000
    PC.481	PC.481	0.9998	0.0000	0.0000
    PC.593	PC.593	0.9999	0.0000	0.0000
    PC.607	PC.607	0.9990	0.0000	0.0000
    PC.634	PC.634	1.0000	0.0000	0.0000
    PC.635	PC.635	0.9963	0.0000	0.0000
    PC.636	PC.636	0.9999	0.0000	0.0000

We see that all of the comparisons indicate strong, statistically-significant correlation. The final column provides p-values that have been corrected for multiple comparisons (though still zero in this case because the uncorrected p-values are zero).

Expected Comparison
-------------------
By specifying `-m expected` to the script, all samples in the first input taxa summary file will be compared to a single "expected" sample in the second input taxa summary file. This comparison mode is especially useful if you need to compare the taxonomic composition of various samples to a sample that has a known taxonomic composition. The second input taxa summary file must contain only a single sample or you must tell the script which one to use (if there are multiple samples) using the --expected_sample_id option. The output files will be in the same format as those seen previously when using paired comparison mode.
