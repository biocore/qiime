args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--source_dir', args)){
    stop("\n\nPlease use '--source_dir' to specify the R source code directory.\n\n")
}
sourcedir <- args[which(args == '--source_dir') + 1]
source(sprintf('%s/util.r',sourcedir))

load.library('optparse')
load.library('DESeq')
load.library('biom')

# make option list and parse command line
option_list <- list(
    make_option(c("--source_dir"), type="character",
        help="Path to R source directory [required]."),
    make_option(c("-i", "--input_path"), type="character",
        help="Input otu table [required]."),
    make_option(c("-o", "--out_path"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("-z", "--DESeq_negatives_to_zero"), type="character", default=NULL,
        help="set the negatives that result from DESeq transformation to zero")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)

# Error checking
if(is.null(opts$input_path)) stop('Please supply an otu table.')


"DESeq" <- function(input_path, out_path, DESeq_negatives_to_zero=NULL, sampleConditions = rep("A", ncol(foo)), 
                method = "blind", sharingMode = "maximum", fitType = "local", locfit_extra_args = list(maxk=300)) {
                foo = read_biom(input_path)
                x = as(biom_data(foo), "matrix")
                # avoid zeros
                x = x + 1
                #Create annotated data.frame with the taxonomy table
                if (!is.null(observation_metadata(foo))) {
                	taxADF = as(data.frame(as(observation_metadata(foo), "matrix"), 
                   		stringsAsFactors = FALSE), "AnnotatedDataFrame")
                	cds = newCountDataSet(x, sampleConditions, featureData = taxADF)
                } else {
	               	cds = newCountDataSet(x, sampleConditions)
	            }
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

DESeq(opts$input_path, opts$out_path, opts$DESeq_negatives_to_zero)
