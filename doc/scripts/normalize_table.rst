.. _normalize_table:

.. index:: normalize_table.py

*normalize_table.py* -- Matrix normalization alternatives to rarefaction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

To perform many downstream analyses after OTU picking (besides
metagenomeSeq's fitZIG and DESeq OTU differential abundance testing), the OTU
matrix must be normalized to account for uneven column (sample) sums that are a
result of most modern sequencing techniques.  These methods attempt to correct
for compositionality too.  Rarefying throws away some data by rarefying to a
constant sum and throwing away extremely low depth samples.

Even with these new normalization techniques, we would recommend throwing away
low depth samples (e.g. less that 1000 sequences/sample).  DESeq/DESeq2 outputs
negative values for lower abundant OTUs as a result of its log transformation.
For most ecologically useful metrics (e.g. UniFrac/Bray Curtis) this presents
problems. No good solution exists at the moment for this issue.  Note that one
is added to the matrix to avoid log(0).  It has been shown that clustering
results can be highly dependent upon the choice of the pseudocount (e.g. should
it be 0.01 instead of 1?), for more information read Costea, P. et al. (2014)
"A fair comparison", Nature Methods.

DESeq/DESeq2 can also have a very slow runtime, especially for larger datasets.
In this script, we implement DESeq2's variance stabilization technique. If you do use these
alternatives to rarefying, we would recommend metagenomeSeq's CSS (cumulative sum
scaling) transformation for those metrics that are abundance-based.  It is not 
recommended to use these new methods with presence/absence metrics, for example
binary Jaccard or unweighted UniFrac. 

For more on metagenomeSeq's CSS, please see Paulson, JN, et al. 'Differential
abundance analysis for microbial marker-gene surveys' Nature Methods 2013.  For DESeq
please see Anders S, Huber W. 'Differential expression analysis for sequence
count data.' Genome Biology 2010.  For DESeq2 please read Love, MI et al. 
'Moderated estimation of fold change and dispersion for RNA-Seq data 
with DESeq2,' Genome Biology 2014.  If you use these methods, please CITE the
appropriate reference as well as QIIME.  For any of these methods, clustering by
sequence depth MUST BE CHECKED FOR as a confounding variable, e.g. by coloring
by sequences/sample on a PCoA plot and testing for correlations between
taxa abundances and sequencing depth with e.g. adonis in `compare_categories.py <./compare_categories.html>`_, 
or `observation_metadata_correlation.py <./observation_metadata_correlation.html>`_.



**Usage:** :file:`normalize_table.py [options]`

**Input Arguments:**

.. note::

	
	**[OPTIONAL]**
		
	-i, `-`-input_path
		Path to the input BIOM file (e.g., the output from OTU picking) or directory containing input BIOM files for batch processing [REQUIRED if not passing -l]
	-o, `-`-out_path
		Output filename for single file operation, or output directory for batch processing [REQUIRED if not passing -l]
	-s, `-`-output_CSS_statistics
		Output CSS statistics file. This will be a directory for batch processing, and a filename for single file operation [default: False]
	-z, `-`-DESeq_negatives_to_zero
		Replace negative numbers produced by the DESeq normalization technique with zeros [default: False]
	-a, `-`-algorithm
		Normalization algorithm to apply to input BIOM table(s). [default: CSS] Available options are: CSS, DESeq2
	-l, `-`-list_algorithms
		Show available normalization algorithms and exit [default: False]


**Output:**

BIOM table with normalized counts.


**Single File CSS Matrix Normalization:**

Normalize a raw (non-normalized/non-rarefied) otu_table.biom using CSS:

::

	normalize_table.py -i otu_table.biom -a CSS -o CSS_normalized_otu_table.biom

**Single File DESeq2 Matrix Normalization:**

Normalize a raw (non-normalized/non-rarefied) otu_table.biom using DESeq2:

::

	normalize_table.py -i otu_table.biom -a DESeq2 -o DESeq2_normalized_otu_table.biom

**Multiple File Matrix Normalization:**

Normalize a folder of raw (non-normalized/non-rarefied) otu tables using e.g. DESeq2:

::

	normalize_table.py -i otu_tables/ -a DESeq2 -o normalized_tables/


