.. _running_supervised_learning:

==========================
Running Supervised Learning
==========================

This document describes how to run supervised classification with QIIME. The goal of supervised classification is to classify new, unlabeled communities based on a set of labeled training communities. See (Knights et al. 2010) for a general discussion of the application of supervised classification to microbiota. Supervised classification using the Random Forests (Breiman, 2001) classifier is implemented in the QIIME script `supervised_learning.py <../scripts/supervised_learning.html>`_. When you run this script you will get several output files:

* `summary.txt`. This gives the estimated generalization error (an estimate of how much error the classifier would have on a novel data set), a description of the method used to estimate the generalization error, and the total number of features used by the model (for Random Forests this is not relevant, as it will usually be close to the total number of OTUs). The estimated generalization error should be compared to the "baseline" error of an uninformed classifier that merely guesses the most common class. For example, if 90% of your samples belong to a single class, then your baseline error is 10%, and you should not be impressed if a classifier achieves, say, 8% error. In contrast, if your data is evenly divided into 10 classes, then your baseline error for random guessing is 90%, and even 50% classifier error is an indication that the classifier is doing quite a bit better than random.
* `cv_probabilities.txt`. Cross-validation estimates of class probabilities for the samples. For each sample, gives the classifier's estimated probability that the sample belonged to each of the possible classes or categories. To avoid overfitting, the estimates for a given sample are always predicted by models that did not contain that sample in their training sets.
* `mislabeling.txt`. For each sample, gives the estimated probability of the alleged class (the class shown in the mapping file), and the probability for the most likely other class. This is just a convenient summary of the information contained in `cv_probabilities.txt`.
* `feature_importance_scores.txt`. The classifier's estimates of the importance (i.e., discriminative power) of each of the features (e.g. OTUs). For Random Forests, the importance reported is the expected mean decrease in accuracy when the feature is ignored. There is no accepted threshold for determining which features are significantly discriminative, as this may depend on the overall error rate and other factors. For example, if the overall error rate is only 5%, then a feature whose removal increases error by 1% (i.e. importance score 0.01) might be considered highly discriminative.
* `params.txt`. Shows any non-default parameters passed to the Random Forests (this should eventually be replaced by a log file once Random Forests is ported to QIIME).


Input files
------------------
This script requires a QIIME OTU table (or equivalent), and a QIIME metadata mapping file. 

Running the script
--------------------------------------------------------------------------------

To run supervised classification on the `QIIME tutorial <./tutorial.html>`_ data set, where the "Treatment" metadata column gives the class labels::

	supervised_learning.py -i otu_table.txt -m Fasting_Map.txt -c Treatment -o ml -v
	
All of the result files described above will be contained in `ml/random_forest`. The `-v` flag causes verbose output, including a trace of the classifier's progress while it is running. This runs Random Forests with the default setting of 500 trees. For larger data sets, the expected generalization error may decrease slightly if more trees are used. To run Random Forests with 1,000 trees, create a `params` file:

.. note::

   * # file params.txt
   * # example params file for random forests classifier
   * # commented lines (starting with "#") are ignored
   * 
   * # number of trees in forest; more improves generalization error estimate
   * params$ntree = 1000
   * 
   * # seed integer for the random number generator (between 0 and 1e9)
   * # can be used to replicate results for stochastic processes
   * params$seed = 0

The parameter `params$ntree` can be used to set the number of trees for Random Forests. The parameter `params$seed` can be used to set a predetermined seed for the random number generator, in order to replicate results with stochastic models. To run supervised classification on the `QIIME tutorial <./tutorial.html>`_ data set with this custom params file (should eventually replaced by command-line parameters once Random Forests is ported to QIIME)::

	supervised_learning.py -i otu_table.txt -m Fasting_Map.txt -c Treatment -o ml -v -p params.txt
	



Cautions
---------
Supervised classification is most useful for larger data sets. When data sets are too small, the estimates of the generalization error, feature importance, and class probabilities may be quite variable. How large a data set needs to be depends on, among other things, how subtle are the differences between classes, and how many noisy features (e.g. OTUs) there are.

Note: we recommend running `single_rarefaction.py` on your OTU table before using it as input to `supervised_learning.py`, to control for variation in sequencing effort.

References
------------
Knights D, Costello EK, Knight R (2010). "Supervised Classification of Human Microbiota". FEMS Microbiology Reviews 35, 343-359

Breiman L (2001). "Random forests". Maching Learning 45: 5–32.





