# Runs random forests on QIIME otu table
# usage:
# R --slave --args -i otus.txt -m map.txt -c Treatment -o rf --qiime_dir $QIIME_HOME < randomforests.r
# 
# print help string:
# R --slave --args -h --qiime_dir $QIIME_HOME < randomforests.r
#
# Requires command-line param --qiime_dir pointing to top-level QIIME dir

# load libraries and source files
args <- commandArgs(trailingOnly=TRUE)
if(!is.element('--qiime_dir', args)){
    stop("\n\nPlease use '--qiime_dir' to specify the QIIME home directory.\n\n")
}
qiimedir <- args[which(args == '--qiime_dir') + 1]
source(sprintf('%s/qiime/support_files/R/loaddata.r',qiimedir))
source(sprintf('%s/qiime/support_files/R/util.r',qiimedir))
source(sprintf('%s/qiime/support_files/R/randomforests_util.r',qiimedir))
load.library('optparse')
load.library('randomForest')

# make option list and parse command line
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
              cv (k-fold cross validation, provides mean and standard deviation of error),
              loo (leave-one-out cross validation, useful for small data sets) [default %default]"),
    make_option(c("-o", "--outdir"), type="character", default='.',
        help="Output directory [default %default]"),
    make_option(c("--nfolds"), type="integer", default=10,
        help="Number of folds in cross-validation (ignored if --errortype is not 'cv') [default %default]"),

    make_option(c("--qiime_dir"), type="character", default=NULL,
        help="QIIME home directory [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
set.warning.visibility(opts$verbose) # hide extraneous library output

# Error checking
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

    if(opts$errortype == 'loo' || opts$nfolds >= length(y)) {
        opts$nfolds <- -1
        error.type <- 'leave-one-out cross validation'
    } else {
        error.type <- sprintf('%d-fold cross validation', opts$nfolds)
    }
    result <- rf.cross.validation(x,y,nfolds=opts$nfolds,
            verbose=opts$verbose,ntree=opts$ntree)
    result$error.type <- error.type
}

save.rf.results(result, opts, colnames(x))
