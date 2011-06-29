# File util.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2011, The QIIME Project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.3.0"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Release"

# Global mapping user's model names to library names
libraries <- list('random_forest'='randomForest','elastic_net'='glmnet')

# Loads and preprocesses data
# returns list of x, y, params
"load.data" <- function(arglist){
    x.fp <- arglist[['-i']]
    map.fp <- arglist[['-m']]
    categ <- arglist[['-c']]
    output.dir <- arglist[['-o']]
    modelnames <- strsplit(arglist[['--models']],',')[[1]]

    # if there are seven arguments, last one is params file: parse it
    if(!is.null(arglist[['--params']])) source(arglist[['--params']])
    # generate a random seed if one wasn't provided; set seed
    if(!is.element('seed', names(params))) params$seed=floor(runif(1,0,1e9))
    set.seed(params$seed)

    model.fcns <- list('random_forest'=train.rf.wrapper)

    # load data
    x <- read.table(x.fp,sep='\t',row.names=1,header=TRUE,check.names=FALSE)
    # remove lineage if present
    lineages <- NULL
    lineage.column <- grep("Lineage", colnames(x))
    if(length(lineage.column) > 0){
        lineages <- as.character(x[,lineage.column])
        names(lineages) <- rownames(x)
        x <- x[,-lineage.column]
    }
    x <- t(x)
    
    map <- read.table(map.fp,sep='\t',row.names=1,header=TRUE,check.names=FALSE)
    y <- as.factor(map[,categ])
    names(y) <- rownames(map)

    # drop NA values
    y <- y[!is.na(y)]

    # keep only rows shared between map and data, make sure order is consistent
    shared.rows <- intersect(rownames(x), names(y))
    x <- x[shared.rows,]
    y <- factor(y[shared.rows])

    # Verify that some rows were shared between map and data file
    if(length(shared.rows) == 0){
        cat('Mapping file and OTU table have no sample IDs in common.\n',
             file=stderr())
        q(save='no',status=1,runLast=FALSE);
    }

    # normalize x (skip samples that sum to 0)
    nonzero.rows <- apply(x,1,sum)>0
    x[nonzero.rows,] <- sweep(x[nonzero.rows,], 1, 
                        apply(x[nonzero.rows,], 1, sum), '/')

    # drop singletons
    singletons <- which(apply(x,2,function(x) sum(x>0)) <= 1)
    if(length(singletons) > 0){
        x <- x[,-singletons]
        lineages <- lineages[-singletons]
    }

    return(list(x=x, 
                y=y,params=params, modelnames=modelnames,
                output.dir=output.dir,
                model.fcns=model.fcns, lineages=lineages))
}

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
#    cv_probabilities.txt
#    feature_importance_scores.txt
#    mislabeling.txt
#    summary.txt
#    params.txt
"save.classification.results" <- function(res, output.dir, model.name, seed=NULL){
    # save the current working directory, change into the output directory
    currwd <- getwd()
    setwd(output.dir)

    # write c.v. (or leave-one-out) probabilities for training data if available
    # include p(mislabel) as last column
    sink('cv_probabilities.txt')
    if(is.null(res$cv.probabilities)){
        cat('This model does not predict cross-validated posterior probabilities.')
    } else {
        cat('SampleID',sprintf('\t%s',colnames(res$cv.probabilities)),'\n')
        write.table(res$cv.probabilities,sep='\t',
            quote=FALSE,col.names=FALSE)
    }
    sink(NULL)

    # write c.v. (or leave-one-out) probabilities for training data if available
    # include p(mislabel) as last column
    sink('mislabeling.txt')
    if(is.null(res$cv.probabilities)){
        cat('This model does not predict cross-validated posterior probabilities.')
    } else {
        cat('SampleID\tP(alleged label)',
            '\tP(second best)',
            '\tP(alleged label)-P(second best)\n', sep='')
        p1 <- get.mislabel.scores(res$y,res$cv.probabilities,type='palleged')
        p2 <- get.mislabel.scores(res$y,res$cv.probabilities,type='p2nd')
        p3 <- p1-p2 
        write.table(cbind(p1,p2,p3),sep='\t',
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
    sink(NULL)

    # write feature importance scores
    sink('feature_importance_scores.txt')
    cat(sprintf('# Feature importance scores = %s\n',
            res$importance$importance.method))
    cat('feature_id\timportance_score\n')
    keepix <- res$importance$scores > 0
    for(i in 1:length(res$importance$scores[keepix])){
        cat(sprintf('%s\t%0.8f\n',
            names(res$importance$scores[keepix])[i], 
            res$importance$scores[keepix][i]))
    }
#~     write.table(res$importance$scores[keepix],sep='\t',
#~         quote=FALSE,col.names=FALSE)
    sink(NULL)
    
    # return to the original working directory
    setwd(currwd)
}

# Saves the results from training a classifier with filter into separate files:
#    cv_probabilities.txt
#    feature_importance_scores.txt
#    mislabeling.txt
#    summary.txt
#    params.txt
"save.filter.results" <- function(res, opts, output.dir, 
            model.name, filter.type){
    # save the current working directory, change into the output directory
    currwd <- getwd()
    setwd(output.dir)

    sink('filter_errors.txt')
    cat('Number of features\t')
    write.table(res$error, sep='\t',quote=FALSE)
    sink(NULL)
    
    best.ix <- which.min(res$error[,1])
    best.n <- as.numeric(rownames(res$error))[best.ix]
    sink('filter_features.txt')
    cat('OTU ID\n')
    write.table(res$features[1:best.n], col.names=FALSE, row.names=FALSE, sep='\t',quote=FALSE)
    sink(NULL)
    
    sink('params.txt')
    cat('# values of all non-default parameters\n')
    cat(sprintf('# method was "%s"\n',model.name))
    cat(sprintf('# filter was "%s"\n',filter.type))
    for(i in 1:length(res$params)){
        cat(sprintf('%s:\t%s\n',names(res$params)[i],res$params[[i]]))
    }
    sink(NULL)

    sink('filter_summary.txt')
    cat(sprintf('Estimated generalization error = %f\n',min(res$error[,1])))
    cat(sprintf('Error estimation method = %s\n',res$err.method))
    cat(sprintf('Optimal feature subset size = %d\n', best.n))
    sink(NULL)

    sink('otu_subset_table.txt')
    cat(sprintf('# OTU subset from supervised_learning.py. Classifier = %s, Filter = %s\n',
        model.name, filter.type))
    cat('OTU ID\t')
    otu_subset <- t(opts$x[,res$features[1:best.n]])
    if(!is.null(opts$lineages)){
        lineage_subset <- opts$lineages[res$features[1:best.n]]
        otu_subset <- cbind(otu_subset, lineage_subset)
        colnames(otu_subset)[ncol(otu_subset)] <- "Consensus Lineage"
    }
    write.table(otu_subset,sep='\t',quote=FALSE)
    sink(NULL)
    # return to the original working directory
    setwd(currwd)
}

# Get probability of mislabeling by several measures
"get.mislabel.scores" <- function(y,y.prob,
        type=c("palleged","p2nd","pcross")[1]){

    # get matrix containing only p(other classes), and containing only p(class)
    mm <- model.matrix(~0 + y)
    y.prob.other <- y.prob * (1-mm)
    y.prob.alleged <- apply(y.prob * mm, 1, max)
    if(type=="palleged"){
        return(invisible(y.prob.alleged))
    }
    else if (type=="p2nd"){
        return(invisible(apply(y.prob.other, 1, max)))
    }
    else if (type=="pcross"){
        return(invisible(apply(y.prob.other, 1, max) - y.prob.alleged))
    }
    else{
        stop("Error: mislabel score type not recognized")
        return(NULL)
    }
}


# get the bss/wss ratio for each column of x
"bss.wss" <- function(x,y){
    classes = levels(y)

    if(class(x) != "matrix" && class(x) != "data.frame"){
        return(bss.wss.vector(x,y))
    }
    if((class(x) == "matrix" || class(x) == "data.frame") && ncol(x)==1){
        return(bss.wss.vector(x,y))
    }
    # get the classwise means for each column
    xkj <- matrix(0, nrow=length(classes), ncol=ncol(x))
    for (k in 1:length(classes)){
        xkj[k,] <- apply(x[y==classes[k], ], 2, mean)
    }

    # get the overall means for each column
    xj <- apply(x, 2, mean)

    # get the between-class sum of squares
    BSS <- rep(0, ncol(x))
    for (k in 1:length(classes)){
        BSS <- BSS + sum(y==classes[k]) * (xkj[k,]-xj)^2
    }

    # get the within-class sum of squares
    WSS <- rep(0, ncol(x))
    for (k in 1:length(classes)){
        WSS <- WSS + apply(t(t(x[y==classes[k],])-xkj[k,])^2, 2, sum)
    }

    K <- length(classes)
    N <- length(y)
    bvar <- BSS/(K-1)
    wvar <- BSS/(N-K)
    fstat <- bvar/wvar
    ratio <- BSS/WSS
    ratio[is.na(ratio)] <- 0
    ix <- sort(ratio,index=T,dec=T)$ix
    invisible(list(bss=BSS,wss=WSS,ratio=ratio,
                   bvar=bvar,
                   wvar=wvar,
                   fstat=fstat,
                   pval=df(fstat,K-1,N-K),
                   k=K,n=N,ix=ix))
}

# assumes each x corresponds to a row of y
# returns the range (max,min) of y-values
"my.error.bars" <- function(x, centers, spread, ...){
    width = min(.010,.25/length(x))
    xlim <- range(x)
    barw <- diff(xlim) * width
    upper <- centers + spread
    lower <- centers - spread

    segments(x, upper, x, lower, ...)
    segments(x - barw, upper, x + barw, upper, ...)
    segments(x - barw, lower, x + barw, lower, ...)
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
