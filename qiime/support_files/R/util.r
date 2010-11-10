# File util.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2010, The QIIME project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.2.0"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Release"

# Global mapping user's model names to library names
libraries <- list('random_forest'='randomForest','elastic_net'='glmnet')

# Attempts to load a given library. If does not exists, fails gracefully
# and prints instructions for installing the library
"load.libraries" <- function(model.name){
    # map user's model name to library name
    lib <- libraries[[model.name]]
    
    # attempt to load the library
    has.library <- library(lib,character.only=TRUE,logical.return=TRUE,
                            verbose=F,warn.conflicts=FALSE)
    
    # if does not exists, fail gracefully
    if(!has.library){
        help_string <- sprintf(
            'To install: open R and run the command "install.packages("%s")"',
            lib)
        cat(sprintf('Library %s not found.\n\n',lib),file=stderr())
        cat(help_string,'\n\n',file=stderr())
        q(save='no',status=2,runLast=FALSE);
    }
}

# Saves the results from training a classifier into separate files:
#    predictions.txt
#    probabilities.txt
#    features.txt
#    summary.txt
#    params.txt
"save.results" <- function(res, output.dir, model.name, seed=NULL){
    # save the current working directory, change into the output directory
    currwd <- getwd()
    setwd(output.dir)
    
    # write predictions for training data
    sink('predictions.txt')
    cat('SampleID\tpredicted_class\n')
    write.table(res$predictions,sep='\t',
        quote=FALSE,col.names=FALSE)
    sink(NULL)
    
    # write probabilities for training data if available
    sink('probabilities.txt')
    if(is.null(res$probabilities)){
        cat('This model does not predict posterior probabilities.')
    } else {
        cat('SampleID',sprintf('\t%s',colnames(res$probabilities)),'\n')
        write.table(res$probabilities,sep='\t',
            quote=FALSE,col.names=FALSE)
    }
    sink(NULL)
    
    # write summary of results
    sink('summary.txt')
    cat(sprintf('Estimated generalization error = %f\n',res$err$err))
    cat(sprintf('Error estimation method = %s\n',res$err$err.method))
    cat(sprintf('Number of features used = %d\n',res$nfeatures))
    sink(NULL)
    sink('params.txt')
    cat('# values of all non-default parameters\n')
    cat(sprintf('# method was "%s"\n',model.name))
    for(i in 1:length(res$params)){
        cat(sprintf('%s:\t%s\n',names(res$params)[i],res$params[[i]]))
    }
    if(!is.null(seed)){
        cat(sprintf('random seed:\t%d\n',seed))
    }
    sink(NULL)

    # write feature importance scores
    sink('features.txt')
    cat(sprintf('# Feature importance scores = %s\n',
            res$importance$importance.method))
    cat('feature_id\timportance_score\n')
    keepix <- res$importance$scores > 0
    write.table(res$importance$scores[keepix],sep='\t',
        quote=FALSE,col.names=FALSE)
    sink(NULL)
    
    # return to the original working directory
    setwd(currwd)
}

