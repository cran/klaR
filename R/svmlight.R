svmlight<-function (x, ...) 
    UseMethod("svmlight")


svmlight.default <- function(x, grouping, temp.dir=NULL, pathsvm=NULL, del=TRUE, 
    svm.options=NULL, prior=NULL, out=FALSE, ...)
{
    y <- grouping
    if(is.null(prior)) prior <- table(y) / length(y)

### Construct Dummymatrix for 1-a-a classification: 
### ncol= no. of classes, (i,j)=+1 iff obj i comes from class j, -1 else 
    ys <- as.factor(y)
    tys <- table(ys)
    ymat <- matrix(-1, nrow = nrow(x), ncol = length(tys))
    lev <- levels(ys)
    ymat[cbind(seq(along = ys), sapply(ys, function(x) which(x == lev)))] <- 1
    counts <- as.vector(tys)
    svm.model <- list()
### call svm_learn
    cmd <- if(is.null(pathsvm)) "svm_learn" else file.path(pathsvm, "svm_learn")
    J <- 1:ncol(ymat)
### "paste" file names for training data and model
    train.filename <- paste(temp.dir, "_train_", J, ".dat", sep = "")
    model.filename <- paste(temp.dir, "_model_", J, ".txt", sep = "")
    PWin <- .Platform$OS.type == "windows"
### for all 1-a-a test call svm_learn
    for (j in J){
    ### construct matrix to learn svm in a format, so that svmlight can read it
        train <- svmlight.file(cbind(ymat[,j], x), train = TRUE)    
    ### save to disk
        write.table(train, file = train.filename[j], row.names = FALSE, 
        col.names = FALSE, quote = FALSE)
        if (PWin) 
            system(paste(cmd, svm.options, train.filename[j], model.filename[j]), 
                show.output.on.console = out)
        else 
            system(paste(cmd,svm.options, train.filename[j], model.filename[j]))
    ### store learned model
        svm.model[[j]] <- readLines(model.filename[j])
    }
    if (del) 
       file.remove(c(train.filename, model.filename))

    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    structure(list(prior = prior, counts = counts, lev = lev, 
        temp.dir = temp.dir, pathsvm = pathsvm, del = del, 
        svm.model = svm.model, svm.options = svm.options, call = cl), class = "svmlight")
}



#### svmlight interface for different calls. Copied and adapted from lda.xxx
svmlight.formula <- function(formula, data = NULL, ..., subset, na.action = na.fail) 
{
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) 
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    if (xint > 0) 
        x <- x[, -xint, drop = FALSE]
    res <- svmlight.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}

svmlight.matrix <- function(x, grouping, ..., subset, na.action = na.fail) 
{
    if (!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if (!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x), 
            class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
    res <- NextMethod("svmlight")
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res
}

svmlight.data.frame <- function (x, ...) 
{
   res <- svmlight.matrix(structure(data.matrix(x), class = "matrix"), 
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res
}


### predict method for svmlight
predict.svmlight <- function(object, newdata, scal = TRUE, ...)
{
### copied and adapted form predict.lda
    if (!inherits(object, "svmlight")) 
        stop("object not of class svmlight")
    if (!is.null(Terms <- object$terms)) {
        if (missing(newdata)) 
            newdata <- model.frame(object)
        else {
            newdata <- model.frame(as.formula(delete.response(Terms)), 
                newdata, na.action = function(x) x, xlev = object$xlevels)
        }
        x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch = 0)
        if (xint > 0) 
            x <- x[, -xint, drop = FALSE]
    }
    else {
        if (missing(newdata)) {
            if (!is.null(sub <- object$call$subset)) 
                newdata <- eval.parent(parse(text = paste(deparse(object$call$x, 
                  backtick = TRUE), "[", deparse(sub, backtick = TRUE), 
                  ",]")))
            else newdata <- eval.parent(object$call$x)
            if (!is.null(nas <- object$call$na.action)) 
                newdata <- eval(call(nas, newdata))
        }
        if (is.null(dim(newdata))) 
            dim(newdata) <- c(1, length(newdata))
        x <- as.matrix(newdata)
    }
#######################################

### save test data on disk
    x <- svmlight.file(x, train = FALSE)
    test.filename <- paste(object$temp.dir, "_test_.dat", sep = "")
    write.table(x, file = test.filename, row.names = FALSE, 
        col.names = FALSE, quote = FALSE)
    werte <- NULL
    cmd <- if(is.null(object$pathsvm)) "svm_classify" 
           else file.path(object$pathsvm, "svm_classify")    
    J <- seq(along = object$counts)
    model.filename <- paste(object$temp.dir, "_model_", J, ".txt", sep = "")
    pred.filename <- paste(object$temp.dir, "_pred_", J, ".txt", sep = "")
### for all learned models predict call svm_classify
    for (j in J){
        writeLines(object$svm.model[[j]], model.filename[j])
        system(paste(cmd, test.filename, model.filename[j], pred.filename[j]))
    ### read predicted values
        prognose <- read.table(pred.filename[j], header = FALSE)[ , 1]
        werte <- cbind(werte, prognose)
    }
    if (object$del) 
        file.remove(c(test.filename, pred.filename, model.filename))
### choose class with highest decision value f(x)
    classes <- factor(max.col(werte), levels = seq(along = object$lev), 
        labels = object$lev)
    colnames(werte) <- object$lev
    if (scal) werte <- e.scal(werte)$sv
    return(list(class = classes, posterior = werte))
}

svmlight.file <- function(x, train = FALSE,...)
{
### write data to the following svmlight format. See: http://svmlight.joachims.org/
# <line> .=. <target> <feature>:<value> <feature>:<value> ... <feature>:<value>
# <target> .=. +1 | -1 | 0 | <float> 
# <feature> .=. <integer> | "qid"  ### # qid not supported
# <value> .=. <float>
# The target value and each of the feature/value pairs are separated by a space character. 
# Feature/value pairs MUST be ordered by increasing feature number. Features with value zero can be skipped.
# In classification mode, the target value denotes the class of the example. +1 as the target value marks a positive example, 
# -1 a negative example respectively. So, for example, the line
# -1 1:0.43 3:0.12 9284:0.2
# specifies a negative example for which feature number 1 has the value 0.43, feature number 3 has the value 0.12, 
# feature number 9284 has the value 0.2, and all the other features have value 0. 
# A class label of 0 indicates that this example should be classified using transduction.
# The predictions for the examples classified by transduction are written to the file specified through the -l option. 
# The order of the predictions is the same as in the training data. 
# In regression mode, the <target> contains the real-valued target value.
    if(is.vector(x)) x <- t(x)
    erg <- x
    sn <- 1:nrow(x) 
    if(!train) erg[sn, 1] <- paste("1:", x[sn, 1], sep = "")
    if(ncol(x) > 1){
        j <- 2:ncol(x)
        erg[ , -1] <- matrix(paste(j - train, t(x[,j]), sep = ":"), ncol = ncol(x)-1, byrow = TRUE)
    }
    return(erg)
} 
