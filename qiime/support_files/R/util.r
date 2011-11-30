# File util.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2011, The QIIME Project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.3.0-dev"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Development"

# Attempts to load a given library. If does not exists, fails gracefully
# and prints instructions for installing the library
"load.libraries" <- function(model.name){
    # map user's model name to library name
    lib <- libraries[[model.name]]
    
    # if R_LIBRARY_PATH environment variable is set, add it to the library paths
    lib.loc <- .libPaths()
    envvars <- as.list(Sys.getenv())
    if(is.element('R_LIBRARY_PATH', names(envvars))){
        lib.loc <- c(envvars[['R_LIBRARY_PATH']], lib.loc)
    }
    
    # attempt to load the library, suppress warnings
    options(warn=-1)
    has.library <- library(lib,character.only=TRUE,logical.return=TRUE,
                         verbose=F,warn.conflicts=FALSE,lib.loc=lib.loc)
    options(warn=0)
    
    # if does not exists, fail gracefully
    if(!has.library){
        help_string1 <- sprintf(
            'To install: open R and run the command "install.packages("%s")".', 
            lib)
        cat(sprintf('\n\nLibrary %s not found.\n\n',lib),file=stderr())
        cat(help_string1,'\n\n',sep='',file=stderr())
        
        help_string2 <- sprintf(
"If you already have the %s package installed in a local directory,
please store the path to that directory in an environment variable
called \"R_LIBRARY_PATH\". This may be necessary if you are running
QIIME on a cluster, and the cluster instances of R don't know about
your local R libraries. If you don't know your R library paths, you
can list them by opening R and running with the command, \".libPaths()\".
The current R instance knows about these paths:
[%s]", lib, paste(.libPaths(),collapse=', '))

        cat(help_string2,'\n\n',file=stderr())
        q(save='no',status=2,runLast=FALSE);
    }
}

# Hide warnings
"hide.warnings" <- function(hide=TRUE){
    if(hide) {
        sink(file('/dev/null',open='w'),type='message')
    } else {
        sink(NULL,type='message')
    }
}

# Get probability of mislabeling by several measures
# returns matrix of p(alleged), max(p(others)), p(alleged) - max(p(others))
"get.mislabel.scores" <- function(y,y.prob){
    result <- matrix(0,nrow=length(y),ncol=3)
    # get matrices containing only p(other classes), and containing only p(class)
    mm <- model.matrix(~0 + y)
    y.prob.other.max <- apply(y.prob * (1-mm),1,max)
    y.prob.alleged <- apply(y.prob * mm, 1, max)
    result <- cbind(y.prob.alleged, y.prob.other.max, y.prob.alleged - y.prob.other.max)
    rownames(result) <- rownames(y.prob)
    colnames(result) <- c('P(alleged label)','P(second best)','P(alleged label)-P(second best)')
    return(result)
}

"balanced.folds" <- function(y, nfolds=10){
    # Get balanced folds where each fold has close to overall class ratio
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


