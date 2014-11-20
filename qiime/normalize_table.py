#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

import rpy2.robjects as robjects
from os.path import exists, splitext, join, isdir
from os import makedirs, listdir, remove
from biom import load_table

__author__ = "Sophie Weiss"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Sophie Weiss", "Joseph Paulson"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Sophie Weiss"
__email__ = "sophie.sjw@gmail.com"




robjects.r('''
        CSS <- function(input_path, out_path, output_CSS_stats=NULL) {
            library(metagenomeSeq)
            library(biom)
            x = read_biom(input_path)
            obj = biom2MRexperiment(x)
            p = cumNormStatFast(obj)
            data = cumNorm(obj, p = p)
            mat = MRcounts(data, norm = TRUE, log = TRUE)
            if (!is.null(output_CSS_stats)) {
                exportStats(obj, p=p, file = file.path(output_CSS_stats))
            }
            write_biom(make_biom(mat), out_path)
        }
        ''')


CSS = robjects.r('CSS')

def normalize_CSS(input_path, out_path, output_CSS_stats):
    """performs metagenomeSeq's CSS normalization on a single raw abundance OTU matrix
    """
    base_fname, ext = splitext(input_path)
    json_infile = base_fname+'_json.biom'
    open(str(json_infile),'w').write(load_table(input_path).to_json('forR'))

    if output_CSS_stats==False:
        CSS(json_infile, out_path)
        remove(json_infile)
    else:
        base_fname, ext = splitext(out_path)
        stats_outfile = base_fname+'_CSS_stats.txt'
        CSS(json_infile, out_path, output_CSS_stats=stats_outfile)
        remove(json_infile)


def multiple_file_normalize_CSS(input_dir, output_dir, output_CSS_stats):
    """performs metagenomeSeq's CSS normalization on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        original_fname = base_fname+'.biom'
        json_fname = base_fname+'_json.biom'
        hdf5_infile = join(input_dir, original_fname)
        json_infile = join(input_dir, json_fname)
        open(str(json_infile),'w').write(load_table(hdf5_infile).to_json('forR'))
        outfile = join(output_dir, 'CSS_'+base_fname+'.biom')
        if output_CSS_stats==False:
            CSS(json_infile, outfile)
            remove(json_infile)
        else:
            stats_outfile = join(output_dir, 'CSS_stats_'+base_fname+'.txt')
            CSS(json_infile, outfile, output_CSS_stats=stats_outfile)
            remove(json_infile)



robjects.r('''
        DESeq <- function(input_path, out_path, DESeq_negatives_to_zero=NULL, sampleConditions = rep("A", ncol(foo)), 
                method = "blind", sharingMode = "maximum", fitType = "local", locfit_extra_args = list(maxk=300)) {
                library("DESeq")
                library("biom")
                foo = read_biom(input_path)
                x = as(biom_data(foo), "matrix")
                # avoid zeros
                x = x + 1
                # Create annotated data.frame with the taxonomy table
                taxADF = as(data.frame(as(observation_metadata(foo), "matrix"), 
                   stringsAsFactors = FALSE), "AnnotatedDataFrame")
                cds = newCountDataSet(x, sampleConditions, featureData = taxADF)
                # First estimate library size factors
                cds = estimateSizeFactors(cds)
                # Variance estimation, passing along additional options
                cds = estimateDispersions(cds, method, sharingMode, fitType, locfit_extra_args)
                # Determine which column(s) have the dispersion estimates
                dispcol = grep("disp_", colnames(fData(cds)))
                # Enforce that there are no infinite values in the dispersion estimates
                if (any(!is.finite(fData(cds)[, dispcol]))) {
                    fData(cds)[which(!is.finite(fData(cds)[, dispcol])), dispcol] <- 0
                }
                vsmat = exprs(varianceStabilizingTransformation(cds))
                if (!is.null(DESeq_negatives_to_zero)) {
                    vsmat[vsmat<0]=0
                }
                DESeq_otu_table <- make_biom(vsmat)
                write_biom(DESeq_otu_table, out_path)
            } 
        ''')


DESeq = robjects.r('DESeq')

def normalize_DESeq(input_path, out_path, DESeq_negatives_to_zero):
    """performs DESeqVS normalization on a single raw abundance OTU matrix
    """
    base_fname, ext = splitext(out_path)
    json_infile = base_fname+'_json.biom'
    open(str(json_infile),'w').write(load_table(input_path).to_json('forR'))

    if DESeq_negatives_to_zero==False:
        DESeq(json_infile, out_path)
        remove(json_infile)
    else:
        base_fname, ext = splitext(out_path)
        DESeq(input_path, out_path, DESeq_negatives_to_zero=True)
        remove(json_infile)


def multiple_file_normalize_DESeq(input_dir, output_dir, DESeq_negatives_to_zero):
    """performs DESeqVS normalization on a directory of raw abundance OTU matrices
    """
    if not exists(output_dir):
        makedirs(output_dir)
    file_names = listdir(input_dir)
    file_names = [fname for fname in file_names if not (fname.startswith('.')\
        or isdir(fname))]

    for fname in file_names:
        base_fname, ext = splitext(fname)
        original_fname = base_fname+'.biom'
        json_fname = base_fname+'_json.biom'
        hdf5_infile = join(input_dir, original_fname)
        json_infile = join(input_dir, json_fname)
        open(str(json_infile),'w').write(load_table(hdf5_infile).to_json('forR'))
        outfile = join(output_dir, 'DESeqVS_'+base_fname+'.biom')
        if DESeq_negatives_to_zero==False:
            DESeq(json_infile, outfile)
            remove(json_infile)
        else:
            DESeq(json_infile, outfile, DESeq_negatives_to_zero=True)
            remove(json_infile)


def algorithm_list():
    """ returns list of normalization algorithms from qiime.normalize_table
    """
    return ['CSS', 'DESeq']

if __name__ == "__main__":
    main()
