NaiveBayes <- function (x, ...) 
    UseMethod("NaiveBayes")

NaiveBayes.formula <- function (formula, data, ..., subset, na.action = na.pass) 
{
    call <- match.call()
    Yname <- as.character(formula[[2]])
    if (is.data.frame(data)) {
        m <- match.call(expand = FALSE)
        m$... <- NULL
        m$na.action <- na.action
        m[[1]] <- as.name("model.frame")
        m <- eval(m, parent.frame())
        Terms <- attr(m, "terms")
        if (any(attr(Terms, "order") > 1)) 
            stop("NaiveBayes cannot handle interaction terms")
        Y <- model.extract(m, "response")
        X <- m[ , -attr(Terms, "response")]
        return(NaiveBayes(X, Y, ...))
    }
    else stop("NaiveBayes formula interface handles data frames only")
}


NaiveBayes.default <- function (x, grouping, prior = NULL, usekernel = FALSE, ...) 
{
    x <- data.frame(x)
    if (is.null(prior)) apriori <- table(grouping) / length(grouping)
    else apriori <- as.table(prior / sum(prior))
    call <- match.call()
    #Yname <- deparse(substitute(grouping))
    Yname<- "grouping"
    est <- function(var) 
      if(is.numeric(var)) {
        if (usekernel)
            lapply(split(var, grouping), FUN = function(xx) density(xx, ...))
        else
            cbind(tapply(var, grouping, mean), tapply(var, grouping, sd))
      }
    else {
        tab <- table(grouping, var)
        tab / rowSums(tab)
    }
    tables <- lapply(x, est)
    #for (i in 1:length(tables)) 
    #    {names(dimnames(tables[[i]])) <- c(Yname, colnames(x)[i])}
    names(dimnames(apriori)) <- Yname
    structure(list(apriori = apriori, tables = tables, levels = levels(grouping), 
        call = call, x = x, usekernel = usekernel, varnames = colnames(x)), 
        class = "NaiveBayes")
}
