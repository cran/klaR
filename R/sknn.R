sknn<-function (x, ...) 
    UseMethod("sknn")

sknn.default <- function(x,grouping,k=3,...)
{
cl <- match.call()
cl[[1]] <- as.name("sknn")
structure(list(learn = x, grouping = grouping, lev = levels(grouping), k=k, call = cl), class = "sknn")
}


### sknn bei verschiedenen Eingabeformaten:
sknn.formula<-function (formula, data = NULL, ..., subset, na.action = na.fail) 
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
    res <- sknn.default(x, grouping, ...)
    res$terms <- Terms
    cl <- match.call()
    cl[[1]] <- as.name("sknn")
    res$call <- cl
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
        res$na.action <- attr(m, "na.action")
    res
}

sknn.matrix<-function (x, grouping, ..., subset, na.action = na.fail) 
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
    res <- NextMethod("sknn")
    cl <- match.call()
    cl[[1]] <- as.name("sknn")
    res$call <- cl
    res
}

sknn.data.frame<-function (x, ...) 
{
   res <- sknn.matrix(structure(data.matrix(x), class = "matrix"), 
        ...)
    cl <- match.call()
    cl[[1]] <- as.name("sknn")
    res$call <- cl
    res
}


predict.sknn<-function(object, newdata,...)
{
spsknn <- function(x,object)
    {
        abstand<-apply(object$learn, 1, function(y) sum((y-x)^2))
        kdach<-(object$grouping[order(abstand)][1:object$k])
        return(table(kdach))
    } 

if (!inherits(object, "sknn")) 
        stop("object not of class sknn")
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
                newdataa <- eval.parent(parse(text = paste(deparse(object$call$x, 
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

werte<-t(apply(x,1,spsknn,object=object))
werte<-werte/object$k
classes <- factor(max.col(werte), levels = seq(along = object$lev), 
        labels = object$lev)
result<-list(posterior=werte,class=classes)
return(result)
}
