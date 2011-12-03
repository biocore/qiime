# Runs random forests on QIIME otu table
# usage:
# R --slave --args -i otus.txt -m map.txt -c Treatment -o rf < randomforests.r
# 
# print help string:
# R --slave --args -h < randomforests.r
#
# Requires environment variable QIIME_DIR pointing to  top-level QIIME directory.

# Attempts to load a given library. If does not exists, fails gracefully
# and prints instructions for installing the library
"load.library" <- function(lib.name, quietly=TRUE){

    # if R_LIBRARY_PATH environment variable is set, add it to the library paths
    lib.loc <- .libPaths()
    envvars <- as.list(Sys.getenv())
    if(is.element('R_LIBRARY_PATH', names(envvars))){
        lib.loc <- c(envvars[['R_LIBRARY_PATH']], lib.loc)
    }
    
    # attempt to load the library, suppress warnings
    options(warn=-1)
    has.library <- library(lib.name,character.only=TRUE,logical.return=TRUE,
                         verbose=F,warn.conflicts=FALSE,lib.loc=lib.loc, quietly=quietly)
    options(warn=0)
    
    # if does not exists, fail gracefully
    if(!has.library){
        help_string1 <- sprintf(
            'To install: open R and run the command "install.packages("%s")".', 
            lib.name)
        cat(sprintf('\n\nLibrary %s not found.\n\n',lib.name),file=stderr())
        cat(help_string1,'\n\n',sep='',file=stderr())
        
        help_string2 <- sprintf(
"If you already have the %s package installed in a local directory,
please store the path to that directory in an environment variable
called \"R_LIBRARY_PATH\". This may be necessary if you are running
QIIME on a cluster, and the cluster instances of R don't know about
your local R libraries. If you don't know your R library paths, you
can list them by opening R and running with the command, \".libPaths()\".
The current R instance knows about these paths:
[%s]", lib.name, paste(.libPaths(),collapse=', '))

        cat(help_string2,'\n\n',file=stderr())
        q(save='no',status=2,runLast=FALSE);
    }
}

# make option list and parse command line
load.library('optparse', quietly=TRUE)
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
sourcefiles <- c('loaddata.r', 'util.r', 'randomforests_util.r')
for(sourcefile in sourcefiles) source(sprintf('%s/%s',opts$sourcedir, sourcefile))
hide.warnings(!opts$verbose)
load.library('randomForest',quietly=TRUE)

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

    if(opts$errortype == 'loo') {
        nfolds=-1
        error.type <- sprintf('leave-one-out cross validation')
    } else {
        if(opts$errortype == 'cv5') nfolds=5
        if(opts$errortype == 'cv10') nfolds=10
        error.type <- sprintf('%d-fold cross validation', nfolds)
    }
    nfolds <- min(nfolds, length(y))
    result <- rf.cross.validation(x,y,nfolds=nfolds,verbose=opts$verbose,ntree=opts$ntree)
    result$error.type <- error.type
}
print.rf.results(result, opts, colnames(x))
