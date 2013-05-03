# QIIME: Quantitative Insights Into Microbial Ecology


[![Build Status](http://ci.qiime.org/job/QIIME/badge/icon)](http://ci.qiime.org/job/QIIME/)

The official QIIME source code repository. For details on QIIME, see www.qiime.org. 

See the [QIIME GitHub organization](https://github.com/qiime) for related software projects and data.

For questions on QIIME, head to the [QIIME Forum](https://groups.google.com/forum/#!forum/qiime-forum)



# What is different with this build of QIIME?
This build of QIIME extends the base build to include feature selection methods for metagenomic data. The QIIME build with feature selection was developed in the Ecological & Evolutionary Signal Processing and Informatics Lab ([EESI](http://www.ece.drexel.edu/gailr/EESI/)) at Drexel University. The existing QIIME scripts have not been modified, rather, we have added in a few new scripts to implement the desired feature selection routines in the QIIME pipeline. 

Direct any questions to <gregory.ditzler@gmail.com>

## A few notes 
* [PyFeast](https://github.com/mutantturkey/PyFeast) must be installed
* [Biom file](http://biom-format.org/) must be sparse format
* The label file must be separated by tabs and be integers. The user is required to track a mapping between the integers and physical meaning. An example file could look something like: 1\t2\t2\t1\t2\t

## To do 
* The current implementation of Fizzy processes biom files that are formatted in the sparse format. A future build of the file should include dense biom files and the traditional OTU tables. 
* I have not include all of the FEAST feature selection methods; however, implementing this should be quite easy. All of the feature selection algorithms should be included in the unit test (i.e., `./tests/test_fizzy.py`).
* Offer the ability to extract the class labels from the the map file used in the QIIME analysis. This bullet is not very urgent.
* Allow for the label file to contain labels other than integers. 


## Using feature selection
There are a couple of ways that you can use the Fizzy tool: commandline and calling Fizzy from Python. Let us first examine how to call Fizzy from the commandline. First, make sure you have installed the build of QIIME in this repo. Fizzy requires two inputs to be specified: an input file, and label file. The input file must be in [biom format](http://biom-format.org/), and more specificially in sparse format. The label file must be tab-delimited and contain the numerical class labels. For example, 1 could correspond to healthy and 2 could correspond to unhealthy. Given that you have these files properly formatted you can call the Fizzy module from the commandline as follows: 
```bash
fizzy.py -i data.biom -c label.tsv
```
Given this call, Fizzy will fix the number of features being selected to 15 and the mutual information maximization algorithm will be used. The results will be written to `output.txt` in the directory from which `fizzy.py` was called. The user has control over the number of features being selected and the feature selection algorithm by setting the `-k` and `-f` flags, respectively. Furthermore, the user can also control the output file location by setting the setting te `-o` flag. As an example,
```bash
fizzy.py -i data.biom -c label.tsv -k 15 -f jmi -o ~/results/output.txt
```
will run the joint mutual information feature selection algorithm on `data.biom` and `label.tsv`. There will be 15 features selected and saved in `~/results/output.txt`.

Fizzy can also be called from within Python as long as QIIME is installed. For example,  from the Python interpreter you can use,
```python
import qiime.fizzy as fizzy
data_file = open('./test_fizzy/data.biom','U')    # path to the biom file
label_file = open('./test_fizzy/labels.tsv','U')  # path to the label file 
n_select = 25      # number of features to select from data.biom
method = 'mim'     # select the mutual information maximization algorithm
outfile = './test_fizzy/output.txt'               # path to the output file
fizzy.run_feature_selection(data_file, label_file, outfile, method, n_select) # run Fizzy 
f = open('./test_fizzy/output.txt','U')           # open the result file
print f.read()     # print the contents of the result file
```

## Installing PyFeast
You must have [PyFeast](https://github.com/mutantturkey/PyFeast) installed on your machine to be able to take advantage of the Fizzy feature selection tool for QIIME. PyFeast is a interface for the FEAST feature selection toolbox, which was originally written in C with a interface to Matlab. Note that since QIIME already requires [Numpy](http://www.numpy.org/) to be installed, there are no other dependancies for PyFeast. 

To install PyFeast, make a temporary folder and clone the repo from Github
```bash
mkdir tmp
cd tmp
git clone git://github.com/mutantturkey/PyFeast.git
cd PyFeast
```
Build the mutual information toolbox
```bash
cd FEAST/MIToolbox
make
sudo make install
cd ../../
```
Build the FEAST toolbox
```bash
cd FEAST/FSToolbox
make
sudo make install
cd ../../
sudo ldconfig
```
Install the PyFeast module
```bash
python setup.py build
sudo python setup.py install # requires root access
```
