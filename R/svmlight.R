svmlight<-function (x, ...) 
    UseMethod("svmlight")


svmlight.default <- function(x, grouping, temp.dir=NULL, pathsvm=NULL, del=TRUE, 
    svm.options=NULL, prior=NULL, out=FALSE, ...)
{
    y <- grouping
    if(is.null(prior)) prior <- table(y) / length(y)

### Dummymatrix
    ys <- as.factor(y)
    tys <- table(ys)
    ymat <- matrix(-1, nrow = nrow(x), ncol = length(tys))
    lev <- levels(ys)
    ymat[cbind(seq(along = ys), sapply(ys, function(x) which(x == lev)))] <- 1
    counts <- as.vector(tys)
    svm.model <- list()
    cmd <- file.path(pathsvm, "svm_learn")
    J <- 1:ncol(ymat)
    train.filename <- paste(temp.dir, "_train_", J, ".dat", sep = "")
    model.filename <- paste(temp.dir, "_model_", J, ".txt", sep = "")
    PWin <- .Platform$OS.type == "windows"
    for (j in J){
        train <- svmlight.file(cbind(ymat[,j], x), train = TRUE)    
        write.table(train, file = train.filename[j], row.names = FALSE, 
            col.names = FALSE, quote = FALSE)
        if (PWin) 
            system(paste(cmd, svm.options, train.filename[j], model.filename[j]), 
                show.output.on.console = out)
        else 
            system(paste(cmd,svm.options, train.filename[j], model.filename[j]))
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

svmlight.data.frame<-function (x, ...) 
{
   res <- svmlight.matrix(structure(data.matrix(x), class = "matrix"), 
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("svmlight")
    res$call <- cl
    res
}


predict.svmlight <- function(object, newdata, ...)
{
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


    x <- svmlight.file(x, train = FALSE)
    test.filename <- paste(object$temp.dir, "_test_.dat", sep = "")
    write.table(x, file = test.filename, row.names = FALSE, 
        col.names = FALSE, quote = FALSE)
    werte <- NULL
    cmd <- file.path(object$pathsvm, "svm_classify")
    J <- seq(along = object$counts)
    model.filename <- paste(object$temp.dir, "_model_", J, ".txt", sep = "")
    pred.filename <- paste(object$temp.dir, "_pred_", J, ".txt", sep = "")
    for (j in J){
        writeLines(object$svm.model[[j]], model.filename[j])
        system(paste(cmd, test.filename, model.filename[j], pred.filename[j]))
        prognose <- read.table(pred.filename[j], header = FALSE)[ , 1]
        werte <- cbind(werte, prognose)
    }
    if (object$del) 
        file.remove(c(test.filename, pred.filename, model.filename))
    classes <- factor(max.col(werte), levels = seq(along = object$lev), 
        labels = object$lev)
    colnames(werte)<-object$lev
    return(list(posterior = werte, class = classes))
}

svmlight.file <- function(x, train = FALSE,...)
{
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

#svm.light.file<-function(x,train=FALSE)
#{
#if (is.vector(x)) x<-t(x)
#erg<-x
#for (i in 1:nrow(x))
#   {
#   if(!train) erg[i,1]<-paste("1:",x[i,1],sep="")
#   if (ncol(x)>1)     
#   {j <- 2:ncol(x)
#   erg[i,j]<-paste((j-train), x[i,j], sep=":")}
#   }
#return(erg)
#} 
