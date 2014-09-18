# runs cross-validation 
# if predict.fun is NULL, uses S3 predict method
# if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
# ...: additional parameters for train.fun
#
# value:
# y: true values
# predicted: cv predicted values
# probabilities: cv predicted class probabilities (or NULL if unavailable)
# confusion.matrix: confusion matrix (true x predicted)
# nfolds: nfolds
# params: list of additional parameters
# importances: importances of features as predictors
"rf.cross.validation" <- function(x, y, nfolds=10, verbose=FALSE, ...){
    if(nfolds==-1) nfolds <- length(y)
    folds <- balanced.folds(y,nfolds=nfolds)
    result <- list()
    result$y <- as.factor(y)
    result$predicted <- result$y
    result$probabilities <- matrix(0, nrow=length(result$y), ncol=length(levels(result$y)))
    rownames(result$probabilities) <- rownames(x)
    colnames(result$probabilities) <- levels(result$y)
    result$importances <- matrix(0,nrow=ncol(x),ncol=nfolds)
    result$errs <- numeric(length(unique(folds)))

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
        model <- randomForest(x[-foldix,], factor(result$y[-foldix]), importance=TRUE, do.trace=verbose, ...)
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx,nrow=1)
        result$predicted[foldix] <- predict(model, newx)
        probs <- predict(model, newx, type='prob')
        result$probabilities[foldix,colnames(probs)] <- probs
        result$errs[fold] <- mean(result$predicted[foldix] != result$y[foldix])
        result$importances[,fold] <- model$importance[,'MeanDecreaseAccuracy']
    }

    result$nfolds <- nfolds
    result$params <- list(...)
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    return(result)    
}

# Runs standard random forests with out-of-bag error estimation
# This is merely a wrapper that extracts relevant info
# Return values are the same as rf.cross.validation
"rf.out.of.bag" <- function(x,y, verbose=verbose, ...){
    rf.model <- randomForest(x,y,keep.inbag=TRUE,importance=TRUE,do.trace=verbose,...)
    result <- list()
    result$probabilities <- get.oob.probability.from.forest(rf.model,x)
    result$y <- y
    result$predicted <- rf.model$predicted
    result$confusion.matrix <- t(sapply(levels(y), function(level) table(result$predicted[y==level])))
    result$params <- list(ntree=opts$ntree)
    result$errs <- as.numeric(result$predicted != result$y)
    result$importances <- rf.model$importance[,'MeanDecreaseAccuracy']
    return(result)
}

# get probability of each class using only out-of-bag predictions from RF
"get.oob.probability.from.forest" <- function(model,x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- get.oob.votes.from.forest(model,x)
    # convert to probs
    probs <- sweep(votes, 1, apply(votes, 1, sum), '/')
    rownames(probs) <- rownames(x)
    colnames(probs) <- model$classes
    
    return(invisible(probs))
}

# get votes for each class using only out-of-bag predictions from RF
"get.oob.votes.from.forest" <- function(model,x){
    # get aggregated class votes for each sample using only OOB trees
    votes <- matrix(0, nrow=nrow(x), ncol=length(model$classes))
   
    rf.pred <- predict(model, x, type="vote",predict.all=T)
    for(i in 1:nrow(x)){
        # find which trees are not inbag for this sample
        outofbag <- model$inbag[i,]==0
        # get oob predictions for this sample
        votes[i,] <- table(factor(rf.pred$individual[i,][outofbag],levels=model$classes))
    }
    rownames(votes) <- rownames(x)
    colnames(votes) <- model$classes
    
    return(invisible(votes))
}

# prints random forests results file
"save.rf.results" <- function(result, opts, feature.ids){
    save.rf.results.summary(result, opts, outdir=opts$outdir)
    save.rf.results.probabilities(result, outdir=opts$outdir)
    save.rf.results.mislabeling(result, outdir=opts$outdir)
    save.rf.results.importances(result, feature.ids=feature.ids, outdir=opts$outdir)
    save.rf.results.confusion.matrix(result, outdir=opts$outdir)
}

# Print "summary" file
"save.rf.results.summary" <- function(result, opts, filename='summary.txt', outdir='.'){
    err <- mean(result$errs)
    err.sd <- sd(result$errs)
    sink('log.txt');print(names(result));sink(NULL)
    baseline.err <- 1-max(table(y))/length(y)
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat(sprintf('Model\tRandom Forest\n'))
    cat(sprintf('Error type\t%s\n',result$error.type))
    if(opts$errortype == 'oob' || opts$errortype == 'loo'){
        cat(sprintf('Estimated error\t%.5f\n',err))
    } else {
        cat(sprintf('Estimated error (mean +/- s.d.)\t%.5f +/- %.5f\n',err,err.sd))
    }
    cat(sprintf('Baseline error (for random guessing)\t%.5f\n',baseline.err))
    cat(sprintf('Ratio baseline error to observed error\t%.5f\n',baseline.err / err))
    cat(sprintf('Number of trees\t%d\n',result$params$ntree))
    sink(NULL)
}

# Print "probabilities" file
"save.rf.results.probabilities" <- function(result, filename='cv_probabilities.txt', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('#SampleID\t')
    write.table(result$probabilities,sep='\t',quote=F)
    sink(NULL)
}

# Print "mislabeling" file
"save.rf.results.mislabeling" <- function(result, filename='mislabeling.txt', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    sink(filepath)
    cat('#SampleID\t')
    write.table(get.mislabel.scores(result$y,result$probabilities),sep='\t',quote=F)
    sink(NULL)
}

# Print "feature importance scores" file
"save.rf.results.importances" <- function(result, feature.ids, filename='feature_importance_scores.txt', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    if(is.null(dim(result$importances))){
        imp <- result$importances
        imp.sd <- rep(NA,length(imp))
    } else {
        imp <- rowMeans(result$importances)
        imp.sd <- apply(result$importances, 1, sd)
    }
    output.table <- cbind(imp, imp.sd)
    rownames(output.table) <- feature.ids
    output.table <- output.table[sort(imp,dec=T,index=T)$ix,]
    colnames(output.table) <- c('Mean_decrease_in_accuracy','Standard_deviation')

    sink(filepath)
    cat('Feature_id\t')
    write.table(output.table,sep='\t',quote=F)
    sink(NULL)
}

# Print "confusion matrix" file
"save.rf.results.confusion.matrix" <- function(result, filename='confusion_matrix.txt', outdir='.'){
    filepath <- sprintf('%s/%s',outdir,filename)
    
    # add class error column to each row
    x <- result$confusion.matrix
    class.errors <- rowSums(x * (1-diag(nrow(x)))) / rowSums(x)
    output <- cbind(result$confusion.matrix, class.errors)
    colnames(output)[ncol(output)] <- "Class error"
    sink(filepath)
    cat('True\\Predicted\t')
    write.table(output,quote=F,sep='\t')
    sink(NULL)
}