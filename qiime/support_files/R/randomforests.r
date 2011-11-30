# Runs random forests on QIIME otu table
# usage:
# R --slave --args -i otus.txt -m map.txt -c Treatment -o rf < randomforests.r
# 
# print help string:
# R --slave --args -h < randomforests.r
#
# Requires environment variable QIIME_DIR pointing to  top-level QIIME directory.

# load libraries and source files
"load.libraries" <- function(sourcedir, quietly=TRUE){
    sourcefiles <- c('loaddata.r', 'util.r', 'randomforests_util.r')
    for(sourcefile in sourcefiles) source(sprintf('%s/%s',sourcedir, sourcefile))

    if(quietly) hide.warnings()
    library('randomForest',warn.conflicts=FALSE,quietly=quietly)
    library('e1071',warn.conflicts=FALSE,quietly=quietly)
}


# make option list and parse command line
library('optparse',warn.conflicts=FALSE,quietly=TRUE)
option_list <- list(
    make_option(c("--sourcedir"), type="character",
        help="Path to QIIME R source directory [required]."),
    make_option(c("-i", "--otutable"), type="character",
        help="Input otu table [required]."),
    make_option(c("-m", "--mapfile"), type="character",
        help="Input metadata mapping file [required]."),
    make_option(c("-c", "--category"), type="character",
        help="Metadata column header giving cluster IDs [required]"),
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
        help="Print warnings and other additional information"),
    make_option(c("--ntree"), type="integer", default=1000,
        help="Number of trees in forest [default %default]"),
    make_option(c("-e", "--errortype"), type="character", default='oob',
        help="Type of error estimation: oob (out-of-bag, fastest), 
              cv5 (5-fold cross validation, provides mean and standard deviation of error),
              cv10 (10-fold cross validation, provides mean and standard deviation of error),
              loo (leave-one-out cross validation, useful for small data sets) [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default='.',
    help="Output directory [default %default]")
)
opts <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly=TRUE))
load.libraries(opts$sourcedir, quietly=!opts$verbose)
hide.warnings(!opts$verbose)

# File requirements
if(is.null(opts$mapfile)) stop('Please supply a mapping file.')
if(is.null(opts$category)) stop('Please supply a mapping file header.')
if(is.null(opts$otutable)) stop('Please supply an otu table.')

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir,showWarnings=FALSE, recursive=TRUE)

# load qiime data
map <- load.qiime.mapping.file(opts$mapfile)
otus <- load.qiime.otu.table(opts$otutable)
data.list <- remove.nonoverlapping.samples(map=map, otus=otus)

# run random forests
x <- data.list$otus
y <- factor(data.list$map[[opts$category]])
if(opts$errortype == 'oob'){
    result <- rf.out.of.bag(x, y, verbose=opts$verbose, ntree=opts$ntree)
    result$error.type <- 'out-of-bag'
} else {
    if(opts$errortype == 'loo') nfolds=-1
    if(opts$errortype == 'cv5') nfolds=5
    if(opts$errortype == 'cv10') nfolds=10
    nfolds <- min(nfolds, length(y))
    result <- rf.cross.validation(x,y,nfolds=nfolds,verbose=opts$verbose,ntree=opts$ntree)
    result$error.type <- sprintf('%d-fold cross validation',nfolds)
}
print.rf.results(result, opts, colnames(x))
