.. _feature_selection:

========================
Fizzy - Information Theoretic Feature Selection for Metagenomics
========================



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
