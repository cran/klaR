# version vom 4.5.'03                       
# jetzt wird nicht mehr durch Null geteilt. 

stepclass <- function(x, ...)
{
  UseMethod("stepclass")
}


stepclass.formula <- function(formula, data, method, ...)
{
  variables <- dimnames(attributes(terms(formula))$factors)[[1]]
  response <- variables[1]
  discriminators <- variables[-1]
  if(any(discriminators == ".")) {
    exclude <- c(response, discriminators[discriminators != "."])
    discriminators <- colnames(data)[!is.element(colnames(data), exclude)]
  }
  result <- stepclass(x=data[, discriminators], grouping=data[, response], 
                      method=method, ...)
  result$call <- match.call()
  return(result)
}
 

stepclass.default <- function(x, grouping, method, prior=1, 
                              improvement=0.95, maxvar=Inf, start.vars=NULL, direction=c("both","forward","backward"), 
                              fold=10, cv.groups=NULL,
                              output=TRUE, ...)
#  Forward/backward variable selection for classification using any specified classification     
#  function and selecting by estimated (cross-validated) error rate.                             
#  The classification function (e.g. "lda") must have its own `predict' function                 
#  (like "predict.lda" for "lda") that either returns a vector of classifications or a list with 
#  an element "$class" containing that vector instead.                                           
#                
# ARGUMENTS:     
#   data         : data matrix (rows=cases, columns=variables)                                     
#   grouping     : class indicator vector (a factor)                                               
#   method       : classification function (e.g. lda)                                              
#   prior        : class priors for error estimation.                                              
#                  prior==1 indicates equally likely classes (default).                            
#   improvement  : least improvement factor of error rate desired to include any new variable      
#   maxvar       : maximum number of variables in model                                            
#   start.vars   : set variables to start with (indices or names)                                  
#   direction    : "forward", "backward" or "both" (default)                                       
#   fold         : parameter for cross-validation (default is 10-fold);                             
#                  omitted if "cv.groups" is specified.                                            
#   cv.groups    : vector of group indicators for cross-validation. by default assigned automatically.
#   output       : indicator (T/F) for textoutput during computation                               
#   ...          : further parameters passed to classification function ("method"), e.g. priors etc.  
#                
# VALUE:         
#   a list with elements                                                                                   
#    $method          : classification function used (e.g."lda")                                           
#    $start.variables : vector of starting variables.                                                      
#    $process         : data frame showing selection process (included/excluded variables and error rates) 
#    $model           : the final model: data frame with 2 columns; indices & names of variables.          
#    $error.rate      : estimated error rate for final model.                                              
{
  textoutput <- function(rate, variables, into.model=NA, out.of.model=NA, variablenames=varnames)
  # function for output during computation (if requested)
  {
    if (is.na(into.model)) {
      if (!is.na(out.of.model)) in.out <- paste('  out: "',variablenames[out.of.model],'"; ',sep="")
      else in.out <- "  starting"
    }
    else in.out <- paste('  in: "',variablenames[into.model], '"; ', sep="")
    cat(paste("error rate: ", round(rate,5), ";",in.out," ",sep="")) # first line
    if (length(variables)==1) cat(paste("variables (1):", variablenames[variables], "\n"))
    else cat(paste("variables (",length(variables),"):", sep=""), 
             paste(variablenames[variables[-length(variables)]], coll=",", sep=""), 
             variablenames[variables[length(variables)]],"\n")
    invisible()
  }

  cv.rate <- function(vars, data=data, grouping=grouping, method=method, 
                      prior=prior, fold=fold, cv.groups=cv.groups, ...)
  # the actual goal function to be minimized 
  # (cross-validated mean error rate)        
  {
    predicted <- numeric(length(grouping)) # vector of (class number) predictions. 
    for (f in 1:fold){
      train <- (cv.groups != f)
      test  <- which(!train)
      train <- which(train)
      traindat <- data[train, vars]
      if (is.vector(traindat))
        if(is.data.frame(data)) { # turn into data frame 
          traindat <- as.data.frame(matrix(traindat, ncol=length(vars)))
          names(traindat) <- names(data)[vars]
        }
        else { # turn into matrix 
          traindat <- matrix(traindat, ncol=length(vars))
        }
      # MODEL FITTING: 
      
      object <- try(do.call(method, list(traindat, grouping[train], ...)), silent=TRUE)
      if (class(object)!="try-error") { # (no error) 
        testdat <- data[test,vars]
        if (is.vector(testdat))
          if (is.data.frame(data)) { # turn into data frame 
            testdat <- as.data.frame(matrix(testdat, ncol=length(vars)))
            names(testdat) <- names(data)[vars]
          }
          else { # turn into matrix
            testdat <- matrix(testdat, ncol=length(vars))
          }
        # PREDICTION: 
        classi <- try(predict(object, testdat), silent=TRUE)
        if (class(classi)!="try-error") { # (no error) 
          if (is.list(classi)) classi <- classi$class
          predicted[test] <- as.numeric(classi)
        }
      }
    }
    # "predicted" now contains class numbers where classification worked 
    # and zeroes otherwise.                                              
    # Zeroes count as wrong classifications.                             
    group.rates <- NULL
    for (lev in 1:length(levels(grouping))){
      group.rates <- c(group.rates, mean(predicted[as.numeric(grouping)==lev] != lev))
    }
    if (any(predicted==0)) 
        warning(" Error(s) in modeling/prediction step.\n")
    return(as.numeric(t(prior) %*% group.rates))
  }

  min.sec <- function(seconds) 
  {
    seconds <- if (length(seconds)>=3) seconds[3]
               else seconds[1]
    result <- c(seconds %/% 3600)
    result <- c("hr" =result, 
                "min"=(seconds - result * 3600) %/% 60, 
                "sec"=(seconds - result * 3600) %% 60)
    return(result) 
  }  
  
  data <- x
  rm("x")
  
  # Let's guess whether we should load a package:
  switch(method,
     "lda" = require("MASS"),
     "qda" = require("MASS"),  
     "rpart" = require("rpart"),
     "naiveBayes" = require("e1071"),
     "svm" = require("e1071"))
    
  
  stopifnot(dim(as.data.frame(data))[1] == length(grouping))
  runtime <- proc.time()[3]
  direction <- match.arg(direction)
  grouping <- factor(grouping)
  g <- length(levels(grouping)) # number of groups
  if (is.finite(maxvar)) improvement <- 100 # arbitrary number greater than one. 
  # set some variables... 
  fwd <- (direction=="forward") || (direction=="both")
  bwd <- (direction=="backward") || (direction=="both")
  if ((length(prior)==1) && (prior==1)) prior <- rep(1/g,g)
  if (!is.null(dimnames(data))) varnames <- dimnames(data)[[2]]        # names for 
  else varnames <- paste("var", as.character(1:dim(data)[2]), sep=".") # variables 
  if ((direction=="backward") && is.null(start.vars)) 
    start.vars <- 1:dim(data)[2] #include all
  # `model' = variables IN model  (vector of indices (numbers))
  if (is.character(start.vars)) {
    model <- seq(along = varnames)[is.element(varnames, start.vars)]
    start.vars <- model
  }
  else model <- start.vars
  # `out' = variables NOT IN model  (vector of indices)
  out <- 1:length(varnames); out <- out[!is.element(out,model)]
  finished <- FALSE
  # assign groups for cross-validation: 
  if (is.null(cv.groups)){
    if (fold > length(grouping)) fold <- length(grouping)
    cv.groups <- rep(0,dim(data)[1])
    groupsizes <- c(0,cumsum(summary(grouping)))
    numbers <- c(rep(1:fold, length(grouping) %/% fold), sample(fold, length(grouping) %% fold))
    for (lev in 1:g) {
      index <- which(grouping == factor(levels(grouping)[lev], levels=levels(grouping))) #indices of class "lev"
      cv.groups[index] <- sample(numbers[(groupsizes[lev]+1):groupsizes[lev+1]])
    }
  }
  else {
    cv.groups <- as.numeric(factor(cv.groups[1:length(grouping)]))
    fold <- max(cv.groups)}
  # some text output... 
  if (output) {
    cat(paste(" `stepwise classification', using ", fold, 
        "-fold cross-validated error rate of method `", 
        method, "'.\n", sep=""))
    cat(dim(data)[1], "observations of", dim(data)[2], 
        "variables in", g, "classes; direction:", direction, "\n")
    if (!is.finite(maxvar)) 
        cat("stop criterion: error improvement less than",
            round((1 - improvement) * 100, 2), "%.\n")
    else cat("stop criterion: assemble", maxvar, "best variables.\n")
    if(.Platform$OS.type == "windows") flush.console()
  }
  if (is.null(start.vars)) old.rate <- 1 - max(prior) # error rate achieved by pure guessing
  else { # compute error rate for starting constellation 
    old.rate <- cv.rate(vars=start.vars, data=data, grouping=grouping, method=method, 
                        prior=prior, fold=fold, cv.groups=cv.groups, ...)
    if(output) {
      textoutput(old.rate, model)
      if(.Platform$OS.type == "windows") flush.console()
    }
  }
  result.e <- old.rate
  result.v <- matrix(c("start", "0"),ncol=2)   # errors & variables 
  last.changed <- NA  # last changed variable doesn't need to be tried in following step.
  
  while (!finished) { # start iterations:
    error.rates <- NULL # compute error rates for...
    if(fwd && (length(out[!is.element(out,last.changed)])>=1)) 
        for (tryvar in out[!is.element(out,last.changed)]) { # ...variables OUT OF model 
               newrate <- cv.rate(vars=c(model,tryvar), data=data, grouping=grouping, 
                                  method=method, prior=prior, fold=fold, cv.groups=cv.groups, ...)
               error.rates <- rbind(error.rates, c(var=tryvar, rate=newrate))
             }
    else if(fwd && (length(out[!is.element(out,last.changed)]) == 0))
        error.rates <- rbind(error.rates, c(var=1, rate=1))
    if(bwd && (length(model[!is.element(model,last.changed)])>=2)) # ...variables IN model 
             for (outvar in model[!is.element(model,last.changed)]) {
               trymodel <- model[!is.element(model, outvar)]
               newrate <- cv.rate(trymodel, data=data, grouping=grouping, method=method, 
                                  prior=prior, fold=fold, cv.groups=cv.groups, ...)
               error.rates <- rbind(error.rates, c(var=outvar, rate=newrate))
             }
    else if(bwd && (length(model[!is.element(model,last.changed)]) == 1)) 
        error.rates <- rbind(error.rates, c(var=1, rate=1))
    best <- order(error.rates[,2])[1]
    last.changed <- error.rates[best,1]
    if (error.rates[best,2] >= improvement*old.rate) {
      if (error.rates[best,2] >= old.rate) 
        finished <- TRUE  # no improvement at all (worse than before)
      else {
        if (is.element(error.rates[best,1],model)){
          # `best' variable was in model, so kicking it out at least doesn't make model worse. 
          # (note that it is easier to get kicked out of model than to get included.) 
          index <- is.element(model, error.rates[best,1])
          out <- c(out, model[index])
          model <- model[!index]
          result.v <- rbind(result.v, c("out",out[length(out)]))
          result.e <- c(result.e, error.rates[best,2])
          if (output) {
            textoutput(error.rates[best,2], model, out.of.model=error.rates[best,1])
            if(.Platform$OS.type=="windows") flush.console()
          }
          old.rate <- error.rates[best,2]
        }
        else finished <- TRUE  # improvement not great enough. 
      }
    }
    else { # improvement could be achieved, so variable is either in- or excluded: 
      if (is.element(error.rates[best,1], model)) {  # kick variable out of model  
        index <- is.element(model, error.rates[best,1])
        out <- c(out, model[index])
        model <- model[!index]
        result.v <- rbind(result.v, c("out",out[length(out)]))
        result.e <- c(result.e, error.rates[best,2])
        if (output) {
          textoutput(error.rates[best,2], model, out.of.model=error.rates[best,1])
          if(.Platform$OS.type == "windows") flush.console()
        }
      }
      else {  # include variable in model
        model <- c(model, error.rates[best,1])
        result.v <- rbind(result.v, c("in",model[length(model)]))
        result.e <- c(result.e, error.rates[best,2])
        out <- out[!is.element(out,error.rates[best,1])]
        if (output) {
          textoutput(error.rates[best,2], model, into.model=error.rates[best,1])
          if(.Platform$OS.type=="windows") flush.console()
        }
      }
      old.rate <- error.rates[best,2]
    }
    # if "maxvar" is specified, model size is checked: 
    if (is.finite(maxvar)) finished <- (length(model)>=maxvar)
  }
  runtime <- min.sec(proc.time()[3] - runtime)
  if (output) {cat("\n"); print(runtime); cat("\n")}
  # selection step finished. 
  # compute apparent error rate (aper)... 
  
  object <- try(do.call(method, list(data[,model], grouping, ...)), silent=TRUE)
  if (class(object)!="try-error") { # (no error) 
    classi <- try(predict(object, data[,model]), silent=TRUE)
    if (class(classi)!="try-error") { # (no error) 
      if (is.list(classi)) classi <- classi$class
      group.rates <- NULL
      for (lev in levels(grouping)) 
        group.rates <- c(group.rates, mean(classi[grouping==lev] != lev))
      aper <- as.numeric(t(prior) %*% group.rates)
    }
    else aper <- NA
  }
  else aper <- NA

  model <- sort(model)
  result <- list("call"=match.call(),
                 "method"=method, "start.variables"=start.vars, 
                 "process"=cbind.data.frame("step"=result.v[,1], 
                                            "var"=as.numeric(result.v[,2]), 
                                            "varname"=c("--",varnames[as.numeric(result.v[-1,2])]),
                                            "error.rate"=result.e),
                 "model"=cbind.data.frame("nr"=model, "name"=I(varnames[model])), 
                 "error.rate"=c("crossval"=result.e[length(result.e)], "apparent"=aper),
                 "runtime"=runtime)#, "cv.groups"=cv.groups)
  rownames(result$process) <- as.character(0:(length(result.v[,1])-1))
  if (length(result$start.variables)>0) names(result$start.variables) <- varnames[result$start.variables]
  class(result) <- "stepclass"
  return(result)
}


print.stepclass <- function(x,...)
{
  kommalist <- function(charvec)
  {
    if (length(charvec)==1) cat(charvec)
    else cat(paste(charvec[-length(charvec)], ",", sep=""), charvec[length(charvec)])
  }
  cat(" method      :", x$method, "\n")
  cat(" final model : ")
  kommalist(x$model$name)
  cat("\n")
  cat(" error rate  :", as.character(signif(x$error.rate[1],4)), "\n")
  invisible(x)
}


plot.stepclass <- function(x, ...)
{
  signum <- rep("-", length(x$process$var)-1)
  signum[as.character(x$process$step[-1])=="in"] <- "+"
  change <- c("START", paste(signum, x$process$varname[-1]))
  par(mar=c(10, 4, 4, 2) + 0.1)
  plot(seq(along = x$process[,1]), x$process$error.rate, type="b",
       xlab="", ylab="estimated error rate", xaxt="n", ...)
  axis(1, at = seq(along = x$process$error.rate), labels=change, las=3, ...)
  invisible(x$progress)
}
