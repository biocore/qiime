#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division
from cogent.app.util import ApplicationNotFoundError
from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file_to_dict
import numpy 


__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison", "Gail Rosen"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"


def get_fs_methods():
    """
        get_fs_methods()
        return the feature selection methods that are 
        available for use in a list. note that the options
        are case sensitive. 
    """
    return ['CIFE','CMIM','CondMI','Condred','ICAP','JMI','MIM','MIFS','mRMR']

def parse_biom(biom_f): 
    """
        data, variable_names, sample_ids = parse_biom(biom_f)
        @biom_f - this is the file handle.
        @data_matrix (return) - dense matrix for feature selection
        @variable_names (return) - feature names in a list
        @sample_ids (return) - names of the samples in the 
            database found in biom_f. 
    """
    biom_table = parse_biom_table(biom_f)
    sample_ids = list(biom_table.SampleIds)
    variable_names = list(biom_table.ObservationIds)
    data_matrix = []
    for data in biom_table.iterObservationData():
        data_matrix.append(data)
    data_matrix = numpy.array(data_matrix)

    return data_matrix.transpose(), variable_names, sample_ids


def parse_map_file(map_f, column_name, sample_ids):
    """
        parse_map_file(map_f, column_name, sample_ids)
        @map_f - file handle
        @column_name - name of the column that contains the class 
            labels
        sample_ids - names of the samples in the order
            of which they appear in the data set. 
        @labels - numpy array of class labels
    """
    obj, comm = parse_mapping_file_to_dict(map_f)
    label_full = []
    labels = []

    # grab the class labels which are likely to be in string 
    # format. 

    for id_set in sample_ids:
        if(id_set not in obj):
            raise ValueError("Unknown sample ID supplied (" + str(id_set) + "). Make sure that the ID is in map file you specified")
        if(column_name not in obj[id_set]):
            raise ValueError("Unknown class name supplied ("+ str(column_name) + "). Make sure that the column name is in map file you specified")
        label_full.append(obj[id_set][column_name])
    
    # now that we have the full class labels we need to determine
    # the number of unique classes in the data. if the number of
    # classes is equal to the number of observations, throw an error. 
    # its likely the user does not know what they are doing. 
    unique_classes = numpy.unique(label_full)
    if len(unique_classes) == len(sample_ids):
        raise ValueError("Number of classes is equal to the number of IDs.  The number of classes must be less than the number of IDs in map file that was specified.")
    if len(unique_classes) == 1: 
        raise ValueError("There must be multiple classes specified. Only one unique class was detected.")

    for str_lab in label_full:
        for n,uclass in enumerate(unique_classes):
            if str_lab == uclass:
                labels.append(float(n))
                break
    return numpy.array(labels)


def run_pyfeast(data, labels, features, method='MIM', n_select=15):
    """
        run_pyfeast(data, labels, method)
        @data - numpy data (dense)
        @labels - vector of class labels (discrete)
        @features - list of feature names
        @method - feature selection method
        @n_select - number of features to select

        The feature selection method is based off of the FEAST 
        C variable selection toolbox. 

        Reference:
        Gavin Brown, Adam Pocock, Ming-Jie Zhao, and Mikel Lujan, 
            "Conditional Likelihood Maximisation: A Unifying Framework 
            for Information Theoretic Feature Selection," Journal of 
            Machine Learning Research, vol. 13, pp. 27--66, 2012.
            (http://jmlr.csail.mit.edu/papers/v13/brown12a.html)
    """
    
    try:
        import feast
    except ImportError:
        raise ApplicationNotFoundError("Error loading the PyFeast module. It is likely that you do not have PyFeast installed locally.")

    try:
        fs_method = getattr(feast, method)
    except AttributeError:
        raise AttributeError("Unknown feature selection method is being specified for PyFeast. Make sure the feature selection method being selected is a valid one. ")

    if len(data.transpose()) < n_select:
      raise ValueError("n_select must be less than the number of observations.")

    sf = fs_method(data, labels, n_select)
    reduced_set = []
    for k in range(len(sf)):
        reduced_set.append(features[int(sf[k])])
    return reduced_set


def run_feature_selection(biom_f, map_f, column_name, method='MIM', n_select=15):
    """
        run_feature_selection(biom_f, map_f, column_name, method, n_select)
        @biom_f - handle of the biom file
        @map_f - handle of the csv file
        @column_name - column name containing the class labels found 
            in the map file. 
        @method - feature selection method [see PyFeast docs]
        @n_select - number of features to selection (integer)
    """
    data_matrix, variable_names, sample_ids = parse_biom(biom_f)
    label_vector = parse_map_file(map_f, column_name, sample_ids)
    reduced_set = run_pyfeast(data_matrix, label_vector, variable_names, method, n_select)
    return reduced_set 
