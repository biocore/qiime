#!/usr/bin/env python
# File created on 14 Jul 2014
from __future__ import division

import rpy2.robjects as robjects
from qiime.parse import parse_mapping_file_to_dict
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
        fitZIG <- function(input_path, out_path, mapping_category, subcategory_1, subcategory_2) {
            library(metagenomeSeq)
            library(biom)
            foo = read_biom(input_path)
            MGS = biom2MRexperiment(foo)
            MGS = cumNorm(MGS,p = cumNormStat(MGS))
            samplesToKeep = which(pData(MGS)[,mapping_category]%in%c(subcategory_1,subcategory_2))
            MGS = MGS[samplesToKeep,]
            MGS_category = pData(MGS)[,mapping_category]
            ###if continuous or a factor, then error message??
            mod = model.matrix(~MGS_category)
            settings = zigControl(maxit=1, verbose=FALSE)
            fit = fitZig(obj=MGS, mod=mod, control=settings)
            MRfulltable(fit, number = nrow(assayData(MGS)$counts), file=out_path, group=3)
        }
        ''')


fitZIG = robjects.r('fitZIG')

def DA_fitZIG(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2):
   """perform metagenomeSeq's Zero Inflated Gaussian (ZIG) OTU differential abundance testing"""
   base_fname, ext = splitext(input_path)
   json_infile = base_fname+'_json.biom'
   tmp_bt = load_table(input_path)
   tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
   tmp_bt.add_metadata(tmp_pmf, 'sample')
   open(str(json_infile),'w').write(tmp_bt.to_json('forR'))

   fitZIG(json_infile, out_path, mapping_category, subcategory_1, subcategory_2)
   remove(json_infile)

def multiple_file_DA_fitZIG(input_dir, output_dir, mapping_fp, mapping_category, subcategory_1, subcategory_2):
    """perform metagenomeSeq's Zero Inflated Gaussian (ZIG) OTU differential abundance test on a directory of raw abundance OTU matrices
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
        tmp_bt = load_table(hdf5_infile)
        tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
        tmp_bt.add_metadata(tmp_pmf, 'sample')
        open(str(json_infile),'w').write(tmp_bt.to_json('forR'))
        outfile = join(output_dir, 'fitZIG_DA_'+base_fname+'.txt')

        fitZIG(json_infile, outfile, mapping_category, subcategory_1, subcategory_2)
        remove(json_infile)


robjects.r('''
        DESeq2 <- function(input_path, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots = NULL) {
            library(DESeq2)
            library(biom)
            foo = read_biom(input_path)
            x = as(biom_data(foo), "matrix")
            # avoid zeros
            x = x + 1
            colData <- data.frame(sample_metadata(foo))
            dds <- DESeqDataSetFromMatrix(x, colData, design = as.formula(paste("~",mapping_category)))
            suppressWarnings(dds <- try(DESeq(dds, quiet = TRUE), silent = TRUE))
            if (inherits(dds, "try-error")) {
                # If the parametric fit failed, try the local.
                suppressWarnings(dds <- try(DESeq(dds, fitType = "local", quiet = TRUE), 
                silent = TRUE))
            if (inherits(dds, "try-error")) {
                # If local fails, try the mean
                suppressWarnings(dds <- try(DESeq(dds, fitType = "mean", quiet = TRUE), 
                    silent = TRUE))
                }
            if (inherits(dds, "try-error")) {
                # If still bad, quit with error.
                return(NULL)
                }
            }
            res <- results(dds, contrast = c(mapping_category, subcategory_1, subcategory_2))
            resOrdered <- res[order(res$padj),]
            df1 <- data.frame(resOrdered)
            df1 <- cbind(OTU = rownames(df1), df1)
            write.table(df1, out_path, sep="\t", quote=F, row.names=F)
            # # #add independent filtering?
            if (!is.null(DESeq2_diagnostic_plots)) {
               plotMA(res, ylim = c(-2,2))
               readline("Press <return to continue")
               plotDispEsts(dds, ylim = c(1e-6, 1e1))
            }
        }
        ''')


DESeq2 = robjects.r('DESeq2')

def DA_DESeq2(input_path, out_path, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots):
    """perform DESeq2 negative binomial differential abundance test on a raw abundance OTU matrix
    """
    base_fname, ext = splitext(input_path)
    json_infile = base_fname+'_json.biom'
    tmp_bt = load_table(input_path)
    tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
    tmp_bt.add_metadata(tmp_pmf, 'sample')
    open(str(json_infile),'w').write(tmp_bt.to_json('forR'))

    if DESeq2_diagnostic_plots==False:
        DESeq2(json_infile, out_path, mapping_category, subcategory_1, subcategory_2)
        remove(json_infile)
    else:
        DESeq2(json_infile, out_path, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots=True)
        remove(json_infile)

def multiple_file_DA_DESeq2(input_dir, output_dir, mapping_fp, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots):
    """perform DESeq2 negative binomial differential abundance test on a directory of raw abundance OTU matrices
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
        tmp_bt = load_table(hdf5_infile)
        tmp_pmf, _ = parse_mapping_file_to_dict(mapping_fp)
        tmp_bt.add_metadata(tmp_pmf, 'sample')
        open(str(json_infile),'w').write(tmp_bt.to_json('forR'))
        outfile = join(output_dir, 'DESeq2_DA_'+base_fname+'.txt') 

        if DESeq2_diagnostic_plots==False:
            DESeq2(json_infile, outfile, mapping_category, subcategory_1, subcategory_2)
            remove(json_infile)
        else:
            print "now showing diagnostic plots for %s " %(original_fname)
            DESeq2(json_infile, outfile, mapping_category, subcategory_1, subcategory_2, DESeq2_diagnostic_plots=True)
            remove(json_infile)

def algorithm_list():
    """ returns list of differential abundance detection algorithms from qiime.differential_abundance
    """
    return ['metagenomeSeq_fitZIG', 'DESeq2_nbinom']

if __name__ == "__main__":
    main()