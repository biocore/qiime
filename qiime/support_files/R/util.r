# File util.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2011, The QIIME Project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.9.1-dev"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
#

# Attempts to load a given library. If does not exists, fails gracefully
# and prints instructions for installing the library
"load.library" <- function(lib.name, quietly=TRUE, bioconductor=FALSE){

    include.custom.library.path() # ensure that custom library is included

    # attempt to load the library, suppress warnings if needed
    warnings.visible <- get.warning.visibility()
    if(quietly && warnings.visible) set.warning.visibility(FALSE)
    library.result <- try (
            library(lib.name,character.only=TRUE,logical.return=TRUE,
                     verbose=F,warn.conflicts=FALSE),
            silent=FALSE
        )
    if(quietly && warnings.visible) set.warning.visibility(TRUE)

    if(bioconductor){
        install.command <- sprintf('source("http://bioconductor.org/biocLite.R"); biocLite("%s")', lib.name)
    } else {
        install.command <- sprintf('install.packages("%s")', lib.name)
    }
    # if library does not exist or failed, fail gracefully
    if(class(library.result)=='try-error'){
        cat(sprintf('\n\nError encounted loading library %s:\n\n',lib.name),sep='',file=stderr())
        cat(library.result[1],'\n\n',sep='',file=stderr())
    } else if(!library.result){
        help_string1 <- sprintf(
            'To install, open R and run the command:\n  %s', install.command)
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

# Ensure that custom library path is included in search path
# default environment variable name is 'R_LIBRARY_PATH'
"include.custom.library.path" <- function(env.var.name = 'R_LIBRARY_PATH'){
    # if R_LIBRARY_PATH environment variable is set, add it to the library paths
    envvars <- as.list(Sys.getenv())
    if(is.element(env.var.name, names(envvars))){
        custom.path <- envvars[[env.var.name]]
        if(!is.element(custom.path, .libPaths())){
            # if custom path is not already in the path list, add it
            .libPaths(c(custom.path, .libPaths()))
        }
    }
}

# Show warning
"get.warning.visibility" <- function(){
    if(sink.number('message') > 2)
        return(FALSE)
    else
        return(TRUE)
}

# Show warning
"set.warning.visibility" <- function(show=FALSE){
    # drop previous sink stack
    while(sink.number('message') > 2)
        sink(NULL,type='message')
    if(!show) sink(file('/dev/null',open='w'),type='message')
}

# Get probability of mislabeling by several measures
# returns matrix of p(alleged), max(p(others)), p(alleged) - max(p(others))
# also includes binary "mislabeled" columns at 5%, 10%, ..., 95%, 99%
# 95% means 95% probability of being mislabeled.
"get.mislabel.scores" <- function(y,y.prob){
    result <- matrix(0,nrow=length(y),ncol=3)
    # get matrices containing only p(other classes), and containing only p(class)
    mm <- model.matrix(~0 + y)
    y.prob.other.max <- apply(y.prob * (1-mm),1,max)
    y.prob.alleged <- apply(y.prob * mm, 1, max)
    result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
    rownames(result) <- rownames(y.prob)
    colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
    for(threshold in c(seq(0.05, 0.95, .05), .99)){
         result <- cbind(result, as.character(y.prob.alleged < (1-threshold)))
         colnames(result)[ncol(result)] <- sprintf('mislabeled_probability_above_%.02f',threshold)
    }
    return(result)
}

# Get balanced folds where each fold has close to overall class ratio
"balanced.folds" <- function(y, nfolds=10){
    folds = rep(0, length(y))
    classes = levels(y)
    # size of each class
    Nk = table(y)
    # -1 or nfolds = len(y) means leave-one-out
    if (nfolds == -1 || nfolds == length(y)){
        invisible(1:length(y))
    }
    else{
    # Can't have more folds than there are items per class
    nfolds = min(nfolds, max(Nk))
    # Assign folds evenly within each class, then shuffle within each class
        for (k in 1:length(classes)){
            ixs <- which(y==classes[k])
            folds_k <- rep(1:nfolds, ceiling(length(ixs) / nfolds))
            folds_k <- folds_k[1:length(ixs)]
            folds_k <- sample(folds_k)
            folds[ixs] = folds_k
        }
        invisible(folds)
    }
}
