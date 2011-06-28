.. _neighbor_joining:

.. index:: neighbor_joining.py

*neighbor_joining.py* -- Build a neighbor joining tree comparing samples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Description:**

The input to this step is a distance matrix (i.e. resulting file from `beta_diversity.py <./beta_diversity.html>`_).


**Usage:** :file:`neighbor_joining.py [options]`

**Input Arguments:**

.. note::

	
	**[REQUIRED]**
		
	-i, `-`-input_path
		Input path.  directory for batch processing, filename for single file operation
	-o, `-`-output_path
		Output path. directory for batch processing, filename for single file operation


**Output:**

The output is a newick formatted tree compatible with most standard tree viewing programs. Batch processing is also available, allowing the analysis of an entire directory of distance matrices.


**neighbor joining (nj) cluster (Single File):**

To perform nj clustering on a single distance matrix (e.g.: beta_div.txt, a result file from `beta_diversity.py <./beta_diversity.html>`_) use the following idiom:

::

	neighbor_joining.py -i beta_div.txt -o beta_div_cluster.tre

**neighbor joining (Multiple Files):**

The script also functions in batch mode if a folder is supplied as input. This script operates on every file in the input directory and creates a corresponding neighbor joining tree file in the output directory, e.g.:

::

	neighbor_joining.py -i beta_div_weighted_unifrac/ -o beta_div_weighted_clusters/


