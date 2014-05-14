.. _feature_selection:

==============================================
Selecting Differentiating Features Using Fizzy
==============================================

Extremely high dimensional data sets are quite common in metagenomic data sets. The biologists collecting and analyzing data need efficient methods to determine relationships between classes in a data set and the variables that are capable of differentiating between multiple groups in a study. Feature selection is one such method to reduce the dimensionality of a data set, such that the remaining features provide the highest level of differentiation between the multiple groups in the study. 

The Fizzy tool (available via fizzy.py) provides access to several information theoretic feature selection algorithms described by Brown et al. (2012) by wrapping the PyFeast_ package into the QIIME suite of tools.

.. _PyFeast: https://github.com/EESI/PyFeast

References

- G Brown, A Pocock, M-J Zhao, and M Lujan, "Conditional likelihood maximization: A unifying framework for information theoretic feature selection," Journal of Machine Learning Research, vol. 13, pp. 27–66, 2012.
- G Ditzler, R Polikar, and G Rosen, "Information theoretic feature selection for high dimensional metagenomic data," in International Workshop on Genomic Signal Processing and Statistics, Washington, DC, 2012, pp. 143--146.






What is Required
----------------

* a biom file containing the abundance profiles
* a mapping file containing the meta-data
* the column in the mapping file that contains the class information (e.g., unhealthy or control).
* output file location containing the results from feature selection


Using Feature Selection on Your Data
------------------------------------

Fizzy requires four flags to be specified: a biom-format file, mapping file, a column specifying the labels in the map file, and a file to save the results. The biom-file must be in the standard biom format. The map file must be tab-delimited and contain the class labels in one of the columns. The labels need not be specified with integer, rather a labeling scheme like 'healthy' or 'unhealthy' will work just fine. Given that you have these files properly formatted you can call the Fizzy module from the command line as follows: ::

	fizzy.py -i data.biom -m map.txt -c Class -o output.txt

Given this expression, Fizzy will fix the number of features being selected to 15 and the mutual information maximization algorithm will be used. The class labels are determined by the column 'Class' in the mapping file `map.txt`. The results will be written to `output.txt` in the directory from which `fizzy.py` was called. The user has control over the number of features being selected and the feature selection algorithm by setting the `-k` and `-f` flags, respectively. As an example, ::

	fizzy.py -i data.biom -m map.txt -c Class -k 15 -f JMI -o ~/results/output.txt

will run the joint mutual information feature selection algorithm on `data.biom` and `map.txt`. There will be 15 features selected and saved in `~/results/output.txt`.
