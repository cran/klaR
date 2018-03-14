#################################################################
###
### WOE
###
### Functions to perform WOE transforms factor ~> numeric 
### (for categorical variables and a binary target)
###
###############################################################


### (sub-)functions:
# compute woe - coefficients for a single factor variable x
computewoes <- function(x, y, weights = NULL, adj){
  if(! is.factor(x)) stop("Input variable x must be a factor!")
  if(! is.factor(y)) stop("Target variable y must be a factor!")
  if(sum(table(x,y)==0) > 0) cat("At least one empty cell (class x level) does exists. Zero adjustment applied!\n")
  
  # class wise event rates 
  if(is.null(weights)) xtab <- table(y, x)
  if(!is.null(weights)) xtab <- wtd.table(y, x, weights)
  # correction for empty class levels
  xtab[which(xtab == 0)] <- xtab[which(xtab == 0)] + adj
  if(is.null(weights)) ncl  <- table(x)
  if(!is.null(weights)) ncl <- wtd.table(x, weights = weights)
  
  # compute woes for alle classes
  fxgegy <- xtab
  fy  <- rowSums(fxgegy)
  for(i in 1:2) fxgegy[i,] <- fxgegy[i,] / fy[i] 
  woes <- log(fxgegy[1,]/fxgegy[2,])
  if(any(fxgegy[2,] == 0)) warning("Empty cells result in infinite woes, zeroadj should be specified > 0!")
  
  # add IV 
  # difference of class wise relative frequencies
  bdiff <- fxgegy[1,] - fxgegy[2,]
  # calculate information value
  IV <- sum(bdiff * woes)
  woes <- c(woes, IV)
}


# apply woe transformation to one single variable using prespecified woe coefficients 
applywoes <- function(woe.obj, x.vec){
  # woe.obj: woe coefficients as returned from compute woes => a single element of train.woes
  if ( sum(sapply(levels(x.vec), function(x) return(sum(x == names(woe.obj)) == 0))) > 0 ) stop("Factor Levels do not match!")

  # check whether same level order in woe.obj and x.vec  
  if(any(names(woe.obj) != levels(x.vec))) stop("Some levels of woe object and new data do not match or are not in the same order!")
  xwoe <- woe.obj[as.integer(x.vec)]
  return(as.numeric(xwoe))
}


woe.default <- function(x, grouping, weights = NULL, zeroadj = 0, ids = NULL, appont = TRUE, ...){
  if (is.factor(x)){ 
    warning("Only one single input variable. Variable name in resulting object$woe is only conserved in formula call.")
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x)) stop("x should be of type data frame.")
  if(is.numeric(ids)){
    if(max(ids) > ncol(x)) stop("Uncorrect coloumn ids specified.")    
    fact.ids <- ids
  } 
  
  if(is.character(ids)){
    fact.ids <- which(colnames(x) %in% ids)
    if (length(fact.ids) < 1) stop("Uncorrect variable names specified!")
  }
  
  if(is.null(ids)) fact.ids <- which(sapply(x,is.factor))
  
  # check for unique factors -> no woes 
  if(is.null(weights)) keep.fids <- sapply(x[,fact.ids], function(z) return(sum(table(z)!=0)>1))
  if(!is.null(weights)){
    if(length(weights) != nrow(x)) stop("Lengths of weights and x differ!")
    #require(questionr) # for function wtd.table()
    keep.fids <- sapply(x[,fact.ids], function(z) return(sum(wtd.table(z, weights=weights)!=0)>1))
  }
  if(sum(keep.fids) == 0) stop("All factors with unique levels. No woes calculated!")
  if(any(!keep.fids)){
    cat("Factor(s) with unique (or zero weight) levels (no woes calculated):\n")
    print(names(x)[fact.ids[!keep.fids]])
    fact.ids <- fact.ids[keep.fids]
  }
  if (length(fact.ids) < 1) stop("No factor variables to be transformed!")
  if(!all(sapply(x[,ids], is.factor))) stop("ids should specify only factor variables")
  if (length(table(grouping)) != 2) stop("WOE transformation is only for binary targets!")
  
  
  
  x.fact <- x[, fact.ids]
  # case of more than one factor variable to be transformed
  if(length(fact.ids)>1){
    x.woes <- lapply(x.fact, computewoes, y = grouping, weights = weights, adj = zeroadj)
    # separate woes and IVs
    IVs <-  sapply(x.woes, function(x) return(x[length(x)]))
    x.woes <- lapply(x.woes, function(x) return(x[-length(x)]))
  }
  # case of only one factor variable to be transformed
  if(length(fact.ids)==1){
    x.woes <- computewoes(x=x.fact, y=grouping, weights = weights, adj = zeroadj)
    # separate woes and IVs
    IVs <-  x.woes[length(x.woes)]
    x.woes <- x.woes[-length(x.woes)]
    x.woes <- list(x.woes)
    names(x.woes) <- names(x)[fact.ids]
  }
  
  res <- list("woe" = x.woes, "IV" = IVs)
  class(res) <- "woe"
  
  if(appont){
    xnew <- predict(res, x, ...)
    res$xnew <- xnew
  }
  
  return(res)
}


# definition of method woe
woe <- function (x, ...) 
    UseMethod("woe")


# formula interface
woe.formula <- function(formula, data = NULL, weights = NULL, ...)
{
  res <- NULL
  if(!is.null(weights)){
    if(is.character(weights)){
      # identify variable id of coarse object
      findvar <- sapply(names(data), function(z) return(weights == substr(z,1,nchar(weights))))
      if(sum(findvar) == 0) stop("No variable name matches weights!")
      if(sum(findvar) > 1) stop("Multiple matches for weights!")
      vid <- which(findvar)
      weights <- data[,vid]
      data <- data[,-vid]
      res <- woe(formula = formula, data = data, weights = weights)
    }
  }
  if(is.null(res)){
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval.parent(m$data))) 
      m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.response(m)
    x <- model.frame(Terms, m)
    #cat(as.character(attr(Terms, "variables"))[1:2])    
    x <- x[,-which(names(x) == as.character(attr(Terms, "variables"))[2])] 
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
    #return(list(x, grouping, xvars))
    res <- woe(x, grouping, weights = weights, ...)
    res$terms <- Terms
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    # if x is a single factor variable its name has to be stored within model
    if(is.factor(x)) names(res$woe)[1] <- names(xlev)
    attr(res, "na.message") <- attr(m, "na.message")
    if (!is.null(attr(m, "na.action"))) 
      res$na.action <- attr(m, "na.action") 
  }
  res
}


# predict woes to a data set
predict.woe <- function(object, newdata, replace = TRUE, ...){
  
  if (!is.data.frame(newdata)) stop("newdata should be of type data frame.")
  # object:
  object <- object$woe
  if (sum(sapply(newdata, is.factor)) == 0) stop("No factor variables to be transformed!")
  #fact.ids <- unlist(sapply(names(newdata), function(x) return(which(names(object) == x))))
  fact.ids <- unlist(sapply(names(object), function(x) return(which(names(newdata) == x))))
  
  if (length(fact.ids) < 1) stop("No factor variables to be transformed!")  
  
  fwoe <- names(object)
  fnew <- names(newdata)[which(sapply(newdata, is.factor))]
  fwoeonly <- setdiff(fwoe, intersect(fwoe, fnew))
  fnewonly <- setdiff(fnew, intersect(fwoe, fnew))
  if(length(fnewonly) > 0) cat("No woe model for variable(s):", fnewonly,"\n")
  if(length(fwoeonly) > 0) cat("Variable(s):", fwoeonly," not in newdata.\n")
  
  x.fact <- newdata[, fact.ids]
  
  if(length(fact.ids) > 1){
    x.woes <- as.data.frame(sapply(seq_along(fact.ids), 
                                   function(i) return(applywoes(object[[i]], x.fact[,which(names(x.fact) == names(object)[i])])))) 
    
    names(x.woes) <- paste("woe", names(fact.ids),sep="_")
    x <- data.frame(newdata, x.woes)
  }
  # special case of only one variable (vector) to be transformed
  if(length(fact.ids) == 1){
    x.woes <- applywoes(object[[1]], x.fact) 
    x <- data.frame(newdata, x.woes)
    names(x)[ncol(x)] <- paste("woe", names(fact.ids),sep="_")
  }
  
  if (replace) x <- x[,-fact.ids]  
  return(x)
}



plot.woe <- function(x, type = c("IV", "woes"), ...){
  type <- match.arg(type)
  if(type=="IV"){
    x <- x$IV
    barplot(height = x, ylab = "Information value", ...)
    xmax <- length(x) * 1.2
    lines(c(0, xmax), rep(0.02, 2), col = "red") 
    lines(c(0, xmax), rep(0.10, 2), col = "yellow") 
    lines(c(0, xmax), rep(0.30, 2), col = "green") 
  }
  if(type=="woes"){
    if(!any(names(x) == "xnew")) stop("Plot of type 'woe' only possible for appont == TRUE.")
    vids <- which(substr(names(x$xnew),1,4) == "woe_")
    if(length(vids) == 0) stop("Data contains no woe variables to be displayed.")
    
    par(ask="TRUE")
    n <- nrow(x$xnew)
    for(i in vids){
      tab <- table(x$xnew[,i])/n
      plot(as.numeric(names(tab)), tab, type = "h", xlab = names(x$xnew)[i], 
           ylim=c(0,1), yaxt = "n", ylab ="Relative Frequency", ...)
      axis(2, at = seq(0, 1, 0.2))  
      }
    par(ask="FALSE")
    invisible()
  }
}


print.woe <- function(x, ...){
    cat("Information values of transformed variables: \n\n")
    IV   <- sort(x$IV, decreasing = TRUE)
    names <- names(IV)
    print(cbind(IV))
}
