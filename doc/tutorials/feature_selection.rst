.. _feature_selection:

========================
Fizzy - Information Theoretic Feature Selection for Metagenomics
========================



What is Required
------------------
* a biom file containing the abundance profiles
* a map file containing the meta-data
* the column in the map file that contains the class information (e.g., unhealthy or control).


Using Feature Selection on Your Data
------------------
There are a couple of ways that you can use the Fizzy tool: commandline and calling Fizzy from Python. Let us first examine how to call Fizzy from the commandline. First, make sure you have installed the build of QIIME in this repo. Fizzy requires two inputs to be specified: an input file, and label file. The input file must be in biom format, and more specificially in sparse format. The label file must be tab-delimited and contain the numerical class labels. For example, 1 could correspond to healthy and 2 could correspond to unhealthy. Given that you have these files properly formatted you can call the Fizzy module from the commandline as follows: ::

	fizzy.py -i data.biom -m labels.map -c Class

Given this call, Fizzy will fix the number of features being selected to 15 and the mutual information maximization algorithm will be used. The results will be written to `output.txt` in the directory from which `fizzy.py` was called. The user has control over the number of features being selected and the feature selection algorithm by setting the `-k` and `-f` flags, respectively. Furthermore, the user can also control the output file location by setting the setting te `-o` flag. As an example, ::

	fizzy.py -i data.biom -m labels.map -c Class -k 15 -f JMI -o ~/results/output.txt

will run the joint mutual information feature selection algorithm on `data.biom` and `label.map`. There will be 15 features selected and saved in `~/results/output.txt`.
