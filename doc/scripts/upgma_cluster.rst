.. _upgma_cluster:

.. index:: upgma_cluster.py

*upgma_cluster.py* -- Build a UPGMA tree comparing samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

In addition to using PCoA, it can be useful to cluster samples using UPGMA (Unweighted Pair Group Method with Arithmetic mean, also known as average linkage). As with PCoA, the input to this step is a distance matrix (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_).


**Usage:** :file:`upgma_cluster.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path.  directory for batch processing, filename for single file operation
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation


**Output:**

The output is a newick formatted tree compatible with most standard tree viewing programs. Batch processing is also available, allowing the analysis of an entire directory of distance matrices.


**UPGMA Cluster (Single File):**

To perform UPGMA clustering on a single distance matrix (e.g.: beta_div.txt, a result file from `beta_diversity.py <./beta_diversity.html>`_) use the following idiom:

::

	upgma_cluster.py -i beta_div.txt -o beta_div_cluster.tre

**UPGMA Cluster (Multiple Files):**

The script also functions in batch mode if a folder is supplied as input. This script operates on every file in the input directory and creates a corresponding upgma tree file in the output directory, e.g.:

::

	upgma_cluster.py -i beta_div_weighted_unifrac/ -o beta_div_weighted_clusters/


