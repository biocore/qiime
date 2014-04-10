#!/usr/bin/env python
import sys
import time
from sklearn import manifold
from numpy import asarray
from biom.parse import parse_biom_table
from biom.table import DenseTable
from qiime.format import format_coords

__author__ = "Joshua Haas"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Joshua Haas,","Antonio Gonzalez Pena","Gregory Ditzler"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Joshua Haas"
__email__ = "laptopdude2@gmail.com"

def compute_manifold(file_name,params):

    """compute the specified manifold on the specified file"""

    otu_table = parse_biom_table(open(file_name,"U"))

    samples = otu_table.SampleIds

    if isinstance(otu_table, DenseTable):
        otumtx = otu_table._data.T
    else:
        otumtx = asarray([v for v in otu_table.iterSampleData()])

    alg = opts["algorithm"]

    if alg=="isomap":
        defaults = {"n_neighbors":5,"n_components":3,"eigen_solver":"auto",
            "tol":0,"max_iter":None,"path_method":"auto","neighbors_algorithm":"auto"}
        params = fill_args(defaults,params)
        mapper = manifold.Isomap(
            n_neighbors=params["n_neighbors"],
            n_components=params["n_components"],
            eigen_solver=params["eigen_solver"],
            tol=params["tol"],
            max_iter=params["max_iter"],
            path_method=params["path_method"],
            neighbors_algorithm=params["neighbors_algorithm"])
    elif alg=="lle":
        defaults = {"n_neighbors":5,"n_components":3,"reg"=0.001,"eigen_solver":"auto",
            "tol":1e-06,"max_iter":100,"method":"standard","hessian_tol":0.0001,
            "modified_tol":1e-12,"neighbors_algorithm":"auto","random_state":None}
        params = fill_args(defaults,params)
        mapper = manifold.LocallyLinearEmbedding(
            n_neighbors=params["n_neighbors"],
            n_components=params["n_components"],
            reg=params["reg"],
            eigen_solver=params["eigen_solver"],
            tol=params["tol"],
            max_iter=params["max_iter"],
            method=params["method"],
            hessian_tol=params["hessian_tol"],
            modified_tol=params["modified_tol"],
            neighbors_algorithm=params["neighbors_algorithm"],
            random_state=params["random_state"])
    elif alg=="spectral-embedding":
        defaults = {"n_components":3,"affinity":"nearest_neighbors","gamma":None,
            "random_state":None,"eigen_solver":None,"n_neighbors":None}
        params = fill_args(defaults,params)
        mapper = manifold.SpectralEmbedding(
            n_components=params["n_components"],
            affinity=params["affinity"],
            gamma=params["gamma"],
            random_state=params["random_state"],
            eigen_solver=params["eigen_solver"],
            n_neighbors=params["n_neighbors"])
    elif alg=="ltsa":
        defaults = {"n_neighbors":5,"n_components":3,"reg"=0.001,"eigen_solver":"auto",
            "tol":1e-06,"max_iter":100,"method":"ltsa","hessian_tol":0.0001,
            "modified_tol":1e-12,"neighbors_algorithm":"auto","random_state":None}
        params = fill_args(defaults,params)
        mapper = manifold.LocallyLinearEmbedding(
            n_neighbors=params["n_neighbors"],
            n_components=params["n_components"],
            reg=params["reg"],
            eigen_solver=params["eigen_solver"],
            tol=params["tol"],
            max_iter=params["max_iter"],
            method=params["method"],
            hessian_tol=params["hessian_tol"],
            modified_tol=params["modified_tol"],
            neighbors_algorithm=params["neighbors_algorithm"],
            random_state=params["random_state"])
    elif alg=="mds":
        defaults = {"n_components":3,"metric":True,"n_init":4,"max_iter":300,
            "verbose":0,"eps":0.001,"n_jobs":1,"random_state":None,
            "dissimilarity":"euclidean"}
        params = fill_args(defaults,params)
        mapper = manifold.Isomap(
            n_components=params["n_components"],
            metric=params["metric"],
            n_init=params["n_init"],
            max_iter=params["max_iter"],
            verbose=params["verbose"],
            eps=params["eps"],
            n_jobs=params["n_jobs"],
            random_state=params["random_state"],
            dissimilarity=params["dissimilarity"])
    else:
        print("arg in error, unknown algorithm '"+alg+"'")
        exit(1)

    fit = mapper.fit_transform(otumtx)
    fit /= abs(fit).max()

    eigvals = [3.0,2.0,1.0]
    pcnts = [30.0,20.0,10.0]
    
    return format_coords(samples, fit, eigvals, pcnts)

def multiple_file_manifold(input_dir, output_dir, algorithm):
    
    """perform manifolds on all distance matrices in the input_dir"""
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_names = os.listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')
                                                        or os.path.isdir(fname))]

    for fname in file_names:
        base_fname, ext = os.path.splitext(fname)
        infile = os.path.join(input_dir, fname)
        manifold_res_string = compute_manifold(infile,algorithm)
        outfile = os.path.join(output_dir, algorithm + '_' + base_fname + '.txt')
        outfile = open(outfile, 'w')
        outfile.write(manifold_res_string)
        outfile.close()

def fill_args(defaults,params):
    result = {}
    for key in defaults:
        result[key] = defaults[key]
    if params is not None:
        for key in params:
            result[key] = params[key]
    return result
