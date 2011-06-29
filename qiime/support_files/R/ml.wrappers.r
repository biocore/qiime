# File ml.wrappers.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2011, The QIIME Project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.3.0"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Release"

# Define new Level 3 methods for ml wrappers
"importance" <- function(object,...) UseMethod("importance")
"get.params" <- function(object,...) object$params
"get.err" <- function(object,...) UseMethod("get.err")
"get.probabilities" <- function(object,newx,...) predict(object,newx,type='prob')
"get.cv.probabilities" <- function(object,newx,...) predict(object,newx,type='cv.prob')
"get.cv.votes" <- function(object,newx,...) predict(object,newx,type='cv.votes')

# Produduces a trained object using wrapper for ML package
#
# wrappers must implement several Level 3 functions:
#   importance, get.params, get.err, predict, and get.probabilities
#
"tune.model" <- function(train.fcn, x, y, params=NULL, ...){
    result <- list()
    result$y <- y
    result$model <- train.fcn(x, y, params=params, ...)
    result$params <- get.params(result$model)
    result$importance <- importance(result$model)
    result$nfeatures <- sum(result$importance$scores > 0)
    result$err <- get.err(result$model)
    result$predictions <- predict(result$model, x)
    result$probabilities <- get.probabilities(result$model, x)
    result$cv.probabilities <- get.cv.probabilities(result$model, x)
    result$cv.votes <- get.cv.votes(result$model, x)
    invisible(result)  
}

######################################################
### Level 3 functions for my Random Forests object ###
######################################################

# Produduces a trained rf.wrapper object, using the default randomForest method
# accepted params: ntree, do.trace
"train.rf.wrapper" <- function(x, y, params, ntree=1000, do.trace=FALSE, ...){
    if(is.element("ntree",names(params))) ntree <- params$ntree
    if(is.element("do.trace",names(params))) do.trace <- params$do.trace
    model <- list("rf.wrapper" = randomForest(x,y,
    	     		       ntree=ntree,do.trace=do.trace,
                           importance=TRUE,
			       keep.inbag=TRUE,...))
    model$params <- list(ntree=ntree, ...)
    class(model) <- "rf.wrapper"
    invisible(model)
}

# Returns mean decrease in Gini for all vars
"importance.rf.wrapper" <- function(model, ...){
    invisible(list('scores' = model$rf.wrapper$importance[,'MeanDecreaseAccuracy'],
                   'importance.method' = 
            "Mean decrease in accuracy when the feature is ignored"))
}

# Returns list of OOB error, and description of error method
"get.err.rf.wrapper" <- function(model, ...){
    list('err'=model$rf.wrapper$err.rate[model$params$ntree,'OOB'],
         'err.method' = "Out-of-bag prediction of training data")
}

# Level 3 predict method for randomforest wrapper class
"predict.rf.wrapper" <- function(model, newx, type='response'){
    if(type == 'prob'){
        yhat <- predict(model$rf.wrapper, newx, type='prob')
    } else if(type == 'cv.prob'){
        yhat <- get.oob.probability.from.forest(model$rf.wrapper, newx)
    } else if(type == 'cv.votes'){
        yhat <- get.oob.votes.from.forest(model$rf.wrapper, newx)
    } else {
        yhat <- predict(model$rf.wrapper, newx)
        names(yhat) <- rownames(newx)
    }
    return(yhat)
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

# get probability of each class using only out-of-bag predictions from RF
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
