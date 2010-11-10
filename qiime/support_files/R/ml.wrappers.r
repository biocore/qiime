# File ml.wrappers.r
#
# __author__ = "Dan Knights"
# __copyright__ = "Copyright 2010, The QIIME project"
# __credits__ = ["Dan Knights"]
# __license__ = "GPL"
# __version__ = "1.2.0"
# __maintainer__ = "Dan Knights"
# __email__ = "daniel.knights@colorado.edu"
# __status__ = "Release"

# Define new Level 3 methods for ml wrappers
"importance" <- function(object,...) UseMethod("importance")
"get.params" <- function(object,...) object$params
"get.err" <- function(object,...) UseMethod("get.err")
"get.probabilities" <- function(object,newx,...) predict(object,newx,type='prob')

# Produduces a trained object using wrapper for ML package
#
# wrappers must implement several Level 3 functions:
#   importance, get.params, get.err, predict, and get.probabilities
#
"tune.model" <- function(train.fcn, x, y, params=NULL, ...){
    result <- list()
    result$model <- train.fcn(x, y, params=params, ...)
    result$params <- get.params(result$model)
    result$importance <- importance(result$model)
    result$nfeatures <- sum(result$importance$scores > 0)
    result$err <- get.err(result$model)
    result$predictions <- predict(result$model, x)
    result$probabilities <- get.probabilities(result$model, x)
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
    model <- list("rf.wrapper" = randomForest(x,y,ntree=ntree,do.trace=do.trace,...))
    model$params <- list(ntree=ntree, ...)
    class(model) <- "rf.wrapper"
    invisible(model)
}

# Returns mean decrease in Gini for all vars
"importance.rf.wrapper" <- function(model, ...){
    invisible(list('scores' = model$rf.wrapper$importance[,'MeanDecreaseGini'],
                   'importance.method' = 
            "Mean decrease in Gini index when the feature is ignored"))
}

# Returns list of OOB error, and description of error method
"get.err.rf.wrapper" <- function(model, ...){
    list('err'=model$rf.wrapper$err.rate[model$params$ntree,'OOB'],
         'err.method' = "Out-of-bag prediction of training data")
}

# Level 3 predict method for randomforest wrapper class
"predict.rf.wrapper" <- function(model, newx, type='response', ...){
    if(type == 'prob'){
        yhat <- predict(model$rf.wrapper, newx, type='prob')
    } else {
        yhat <- predict(model$rf.wrapper, newx, ...)
        names(yhat) <- rownames(newx)
    }
    return(yhat)
}


