.. _feature_selection:

========================
Fizzy - Information Theoretic Feature Selection for Metagenomics
========================




Extremely high dimensional data sets are quite common in
metagenomic data sets. The biologists collecting and analyzing data need efﬁcient methods to determine relationships between classes in a data set and the variables that are capable of differentiating between multiple groups in a study. Feature selection is one such method to reduce the dimensionality of a data set, such that the remaining features provide the highest level differentation between the multiple groups in the study. 

References

* G. Ditzler, R. Polikar, and G. Rosen, "Information theoretic feature selection for high dimensional metagenomic data," in International Workshop on Genomic Signal Processing and Statistics, Washington, DC, 2012, pp. 143--146.
* G. Brown, A. Pocock, M.-J. Zhao, and M. Lujan, "Conditional likelihood maximisation: A unifying framework for information theoretic feature selection," Journal of Machine Learning Research, vol. 13, pp. 27–66, 2012.







What is Required
------------------

* a biom file containing the abundance profiles
* a mapping file containing the meta-data
* the column in the mapping file that contains the class information (e.g., unhealthy or control).


Using Feature Selection on Your Data
------------------
First, make sure you have installed the build of QIIME, which has the fizzy.py scripts located in the repo. Fizzy requires two inputs to be specified: a biom-format file, and mapping file. The biom-file must be in the standard biom format. There is no restriction as to whether or not the file is in sparse or dense format. Either file format will work.  The mapping file must be tab-delimited and contain the class labels in one of the columns. The labels need not be specified with integer, rather a labeling scheme like 'healthy' or 'unhealthy' will work just fine. Given that you have these files properly formatted you can call the Fizzy module from the commandline as follows: ::

	fizzy.py -i data.biom -m map.txt -c Class

Given this call, Fizzy will fix the number of features being selected to 15 and the mutual information maximization algorithm will be used. The class labels are determined by the column 'Class' in the mapping file `map.txt`. The results will be written to `output.txt` in the directory from which `fizzy.py` was called. The user has control over the number of features being selected and the feature selection algorithm by setting the `-k` and `-f` flags, respectively. Furthermore, the user can also control the output file location by setting the setting te `-o` flag. As an example, ::

	fizzy.py -i data.biom -m map.txt -c Class -k 15 -f JMI -o ~/results/output.txt

will run the joint mutual information feature selection algorithm on `data.biom` and `map.txt`. There will be 15 features selected and saved in `~/results/output.txt`.
