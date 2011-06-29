# File error.est.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2011, The QIIME Project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.3.0"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Release"
 
# run with: R --vanilla --slave --args otus.txt map.txt Individual < /bio/../code/r/error.est.r

# assumes that NO args are positional
# allows flags without argument
"parse.args" <- function(allowed.args){
    argv <- commandArgs(trailingOnly=T)
    argpos <- 1
    for(name in names(allowed.args)){
        argpos <- which(argv == name)
        if(length(argpos) > 0){
            # test for flag without argument
            if(argpos == length(argv) || substring(argv[argpos + 1],1,1) == '-')
                allowed.args[[name]] <- TRUE
            else {
                allowed.args[[name]] <- argv[argpos + 1]
            }
        }
    }
    return(allowed.args)
}

# runs filter, returns table of nfeatures, err, stderr
"train.filter" <- function(x,y,filter.min=2, filter.max=200, filter.step=1, nreps=10, type="BSSWSS"){
    filter.max <- min(filter.max, ncol(x))
    n.to.try <- seq(filter.min, filter.max, filter.step)
    errs <- matrix(0,length(n.to.try),nreps)
    if(verbose) cat(sprintf('Applying filter from %d to %d in steps of %d...\n',
                    filter.min, filter.max, filter.step))
    for(n in n.to.try) {
        if(verbose) cat(n,'features','\n')
        for(i in 1:nreps){
            if(verbose) cat('\trep', i, 'of', nreps,'\n')
            folds <- balanced.folds(y)
            err <- 0
            for(fold in 1:max(folds)){
                if(verbose) cat('\t\tfold', fold, 'of', 
                    length(unique(folds)),'\n')
                if(arglist[['--filter']] == 'BSSWSS'){
                    features <- bss.wss(x[folds!=fold,],
                                        y[folds!=fold])$ix
                }
                res <- tune.model(model.fcn,
                                  x[folds!=fold,features[1:n]],
                                  y[folds!=fold],params)
                yhat <- predict(res$model,x[folds==fold,features[1:n]])
                err <- err + sum(yhat != y[folds==fold])
            }
            err <- err/length(y)
            errs[which(n.to.try==n),i] <- res$err$err
        }
    }
    means <- apply(errs,1,mean)
    stderrs <- apply(errs,1,sd)/sqrt(nreps)
    features <- colnames(x)[bss.wss(x,y)$ix]
    
    return(list(error=data.frame("Mean error"=means,
                    "Standard error"=stderrs,
                    row.names=n.to.try,check.names=FALSE),
                err.method=res$err$err.method,
                features=features))
}

# parse arg list
allowed.args <- list('-i'=NULL,'-o'=NULL,'-m'=NULL, '-c'=NULL,
                        '--filter'=NULL, '--params'=NULL,
                        '--models'=NULL, '--sourcedir'=NULL,
                        '--filter_min'='10', '--filter_max'='50',
                        '--filter_step'='10', '--filter_reps'='10',
                        '--verbose'=NULL)
arglist <- parse.args(allowed.args)

# load helper files from qiime source directory
source(sprintf('%s/ml.wrappers.r',arglist[['--sourcedir']]))
source(sprintf('%s/util.r',arglist[['--sourcedir']]))

# load data
params <- NULL
opts <- load.data(arglist)
verbose <- !is.null(arglist[['--verbose']])
if(verbose && is.null(arglist[['--filter']])) params$do.trace=TRUE

# do learning, save results
for(modelname in opts$modelnames){
    subdir <- paste(opts$output.dir,modelname,sep='/')
    load.libraries(modelname)
    model.fcn <- opts$model.fcns[[modelname]]

    if(!is.null(arglist[['--filter']])){
        filter.res <- train.filter(opts$x,opts$y,
                                   as.numeric(arglist[['--filter_min']]),
                                   as.numeric(arglist[['--filter_max']]),
                                   as.numeric(arglist[['--filter_step']]),
                                   as.numeric(arglist[['--filter_reps']]))
        filter.res$params <- params
        save.filter.results(filter.res, opts, subdir,
            modelname, arglist[['--filter']])
    } else {
        res <- tune.model(model.fcn,opts$x,opts$y,params)
        res$params <- params        
        # save results
        save.classification.results(res,subdir,modelname,opts$params$seed)
    }
}

# quit without saving history and workspace
q(runLast=FALSE)
