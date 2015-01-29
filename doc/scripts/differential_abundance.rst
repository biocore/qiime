.. _differential_abundance:

.. index:: differential_abundance.py

*differential_abundance.py* -- Identify OTUs that are differentially abundance across two sample categories
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

OTU differential abundance testing is commonly used to identify OTUs that
differ between two mapping file sample categories (i.e. Palm and Tongue body
sites).  We would recommend having at least 5 samples in each category.  These 
methods can be used in comparison to `group_significance.py <./group_significance.html>`_ on a rarefied matrix, 
and we would always recommend comparing the results of these approaches to the
rarefied/`group_significance.py <./group_significance.html>`_ approaches.  We would also recommend treating the
differentially abundant OTUs identified by these (metagenomeSeq zero-inflated Gaussian,
or ZIG, and DESeq2 negative binomial Wald test) techniques with caution, as they assume a
distribution and are therefore parametric tests.  Parametric tests can do poorly if
the assumptions about the data are not met.  These tests are also newer techniques
that are less well tested compared to rarefying with a `group_signficance.py <./group_signficance.html>`_ test.  

The input is a raw (not normalized, not rarefied) matrix having uneven column sums.
With these techniques, we would still recommend removing low depth samples (e.g. below
1000 sequences per sample), and low abundance/rare OTUs from the data set.  The DESeq2 method
should NOT be used if the fit line on the dispersion plot (one of the diagnostic plots
output by the -d, or --DESeq2_diagnostic_plots option) does not look smooth, there are big
gaps in the point spacings, and the fitted line does not look appropriate for the data.

DESeq2 is stronger at very small/smaller data sets, but the run-time beyond 100 
total samples becomes very long.  MetagenomeSeq's fitZIG is a better algorithm for larger
library sizes and over 50 samples per category (e.g. 50 Palm samples), the more the better. 
In simulation, these techniques have higher sensitivity, but sometimes higher false positive
rate compared to the non-parametric tests (e.g. Wilcoxon rank sum) in `group_significance.py <./group_significance.html>`_, 
especially with very uneven library sizes (starting at 2-3 fold difference).  In practice 
and with real data, we do not observe much of a difference between these results and 
the tests in `group_significance.py <./group_significance.html>`_.  

For more on these techniques please see McMurdie, P. and Holmes, S. 'Waste not want not 
why rarefying microbiome data is inadmissible.' PLoS Comp. Bio. 2014.  For more on 
metagenomeSeq and fitZIG, please read Paulson, JN, et al. 'Differential abundance analysis
for microbial marker-gene surveys.'  Nature Methods 2013.  For DESeq2/DESeq please read
Love, MI et al. 'Moderated estimation of fold change and dispersion for RNA-Seq data 
with DESeq2,' Genome Biology 2014.  Anders S, Huber W. 'Differential expression analysis 
for sequence count data.' Genome Biology 2010.  Additionally, you can also read the
vignettes for each of the techniques on the Bioconductor/R websites.  Also, if you use
these methods, please CITE the proper sources above (metagenomeSeq and DESeq2) as well as QIIME.


**Usage:** :file:`differential_abundance.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input_path
		Path to the input BIOM file (e.g., the output from OTU picking) or directory containing input BIOM files for batch processing [REQUIRED if not passing -l]
	-o, `-`-out_path
		Output filename for single file operation, or output directory for batch processing [REQUIRED if not passing -l]
	-a, `-`-algorithm
		Differential abundance algorithm to apply to input BIOM table(s) [default: metagenomeSeq_fitZIG] Available options are: metagenomeSeq_fitZIG, DESeq2_nbinom
	-m, `-`-mapping_file_path
		Path to mapping file [REQUIRED if not passing -l]
	-c, `-`-mapping_file_category
		Mapping file category [REQUIRED if not passing -l]
	-x, `-`-mapping_file_subcategory_1
		Mapping file subcategory [REQUIRED if not passing -l]
	-y, `-`-mapping_file_subcategory_2
		Mapping file subcategory [REQUIRED if not passing -l]
	-l, `-`-list_algorithms
		Show available differential abundance algorithms and exit [default: False]
	-d, `-`-DESeq2_diagnostic_plots
		Show a MA plot - y axis: log2 fold change, x axis: average size factor normalized OTU value. Also show a Dispersion Estimate plot - visualize the fitted dispersion vs. mean relationship [default: False]


**Output:**

The resulting output OTU txt file contains a list of all the OTUs in the input matrix, along with their associated statistics and FDR p-values.


**Single File OTU Differential Abundance Testing with metagenomeSeq_fitZIG:**

Apply metagenomeSeq_fitZIG differential OTU abundance testing to a raw (NOT normalized) BIOM table to test for differences in OTU abundance between samples in the Treatment:Control and Treatment:Fast groups.

::

	differential_abundance.py -i otu_table.biom -o diff_otus.txt -m map.txt -a metagenomeSeq_fitZIG -c Treatment -x Control -y Fast

**Single File OTU Differential Abundance Testing with DESeq2_nbinom:**

Apply DESeq2_nbinom differential OTU abundance testing to a raw (NOT normalized) BIOM table to test for differences in OTU abundance between samples in the Treatment:Control and Treatment:Fast groups, including output of plots.

::

	differential_abundance.py -i otu_table.biom -o diff_otus.txt -m map.txt -a DESeq2_nbinom -c Treatment -x Control -y Fast -d

**Multiple File OTU Differential Abundance Testing with metagenomeSeq_fitZIG:**

Apply metagenomeSeq_fitZIG differential OTU abundance testing to a folder of raw (NOT normalized) BIOM tables to test for differences in OTU abundance between samples in the Treatment:Control and Treatment:Fast groups.

::

	differential_abundance.py -i otu_tables/ -o diff_otus/ -m map.txt -a metagenomeSeq_fitZIG -c Treatment -x Control -y Fast


