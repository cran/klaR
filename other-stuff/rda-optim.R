rda <- function(x, ...) 
{
  UseMethod("rda")
}


rda.formula <- function(formula, data, ...)
{
  variables <- dimnames(attributes(terms(formula, data=data))$factors)[[1]]
  response <- variables[1]
  discriminators <- variables[-1]
  if (any(discriminators==".")) {
    exclude <- c(response, discriminators[discriminators!="."])
    discriminators <- colnames(data)[!is.element(colnames(data), exclude)]
  }
  result <- rda(x=data[, discriminators], grouping=data[, response], ...)
  result$call <- match.call()
  result$varnames <- discriminators # necessary if  length(disriminators)==1 
  return(result)
}


rda.default <- function(x, grouping=NULL, prior=NULL, 
                        gamma=NA, lambda=NA, regularization=c("gamma"=gamma, "lambda"=lambda),
                        crossval=TRUE, fold=10, train.fraction=0.5, 
                        estimate.error=TRUE, output=FALSE,
                        startsimplex=NULL, max.iter=100, trafo=TRUE,
                        simAnn=FALSE, schedule=2, T.start=0.1, halflife=50, zero.temp=0.01, alpha=2, K=100, ...)
{
  classify <- function(dataset, mu, sigma, pooled, gamma, lambda, g, p, n, prior)
  # mu            : matrix with group means as columns 
  # sigma         : (p x p x g)-Array of group covariances 
  # pooled        : pooled covariance matrix 
  # gamma, lambda : regularization parameters 
  # g, p, n       : numbers of classes, variables & observations
  {
    # compute likelihood of each datum for each group: 
    likelihood <- matrix(NA, nrow=n, ncol=g)
    dataset <- matrix(dataset,ncol=p)
    i <- 1
    singu <- FALSE
    while ((i <= g) & (!singu)){ # compute likelihoods group-wise
      reg.cov <- (1-lambda)*sigma[,,i] + lambda*pooled                   # shift towards pooled cov. 
      reg.cov <- if (is.matrix(reg.cov)) # cov is a matrix
                   (1-gamma)*reg.cov + gamma*mean(diag(reg.cov))*diag(p) # shift towards identity 
                 else # one variable  =>  cov is 1x1-matrix:
                   (1-gamma)*reg.cov + gamma*reg.cov                     # shift towards identity 
      # try to invert covariance matrix: 
      inv.trial <- try(solve(reg.cov), silent=TRUE)
      singu <- ((! is.matrix(inv.trial)) || (!is.finite(ldc<-log(1/det(inv.trial)))))
      if (!singu) likelihood[,i] <- prior[i]*dmvnorm(dataset, mu=mu[,i], inv.sigma=inv.trial, logDetCov=ldc)
      i <- i+1
    }
    if (!singu) result <- apply(likelihood,1,function(x){order(x)[g]}) # return classifications
    else result <- rep(0,nrow(dataset)) # return definitely false classifications 
    return(result)
  }

  goalfunc <- function(paramvec=c(0.5,0.5)) 
  # depends on a (2-dim.) parameter vector, so it can be handled by minimization function. 
  # first element: gamma,  second element: lambda. 
  # returns mean misclassification rate for the bootstrap samples. 
  {
    #if (any(paramvec>1)) error.rates<-1+max(abs(paramvec[1]-1) , abs(paramvec[2]-1))
    #else if (any(paramvec<0)) error.rates<-1+abs(min(paramvec[1],paramvec[2]))
    #else
    paramvec<-linkfunc(paramvec)
    {error.rates <- numeric(fold)
    for (i in 1:fold) {
      siggi <- array(covariances[,,,i],c(p,p,g))
      mumu <- array(means[,,i], c(p,g))
      prediction <- classify(data[-train[,i],], 
                             mu=mumu, sigma=siggi, pooled=covpooled[,,i],
                             gamma=paramvec[1], lambda=paramvec[2], g=g, p=p, n=n-sum(train[,i]!=0), prior=prior)
      errors <- prediction != grouping[-train[,i]]
      group.rates <- tabulate(grouping[-train[,i]][errors],g) / test.freq[,i] #error rate by group
      group.rates[test.freq[,i]==0] <- (g-1)/g # conservative estimate for "empty" groups 
      error.rates[i] <- t(prior)%*%group.rates
    }}
    return(mean(error.rates))
  }
  
  crossval.sample <- function(grouping, fold=10)
  # returns more or less equally sized cross-validation-samples. 
  {
    grouping <- factor(grouping)
    g <- length(levels(grouping)) #number of groups
    if (fold > length(grouping)) fold <- length(grouping)
    cv.groups <- rep(0,length(grouping))
    groupsizes <- c(0,cumsum(summary(grouping)))
    numbers <- c(rep(1:fold, length(grouping) %/% fold), 
                 sample(fold, length(grouping) %% fold)) # group numbers to be assigned 
    for (lev in 1:g) {
      index <- which(grouping==factor(levels(grouping)[lev], levels=levels(grouping))) # indices of class "lev"
      cv.groups[index] <- sample(numbers[(groupsizes[lev]+1):groupsizes[lev+1]])
    }
    return(cv.groups)
  }

#                                                      #
#  Beginning of  _M_A_I_N_ _P_R_O_C_E_D_U_R_E_  (RDA)  #
#                                                      #
  data <- x
  rm(x)
  # if `grouping' vector not given, first data column is taken. 
  if (is.null(grouping)) {
    grouping <- data[,1]
    data <- data[,-1]
  }
  data <- as.matrix(data)
  stopifnot(dim(data)[1]==length(grouping))
  grouping <- factor(grouping)
  classes <- levels(grouping)
  if (!is.null(dimnames(data))) varnames <- dimnames(data)[[2]]
  else varnames <- NULL
  grouping <- as.integer(grouping)
  if (is.null(prior)) { 
    prior <- tabulate(grouping)
    prior <- prior / sum(prior) 
  } # frequencies as prior
  else if (all(prior == 1)) 
    prior <- rep(1 / length(classes), length(classes))   # uniform prior
  names(prior) <- classes
  dimnames(data) <- NULL
  g <- max(grouping)             # number of groups 
  p <- ncol(data)                # number of variables 
  n <- length(grouping)          # number of observations 
  if (all(is.finite(regularization)))  # no optimization 
    opti <- list(minimum=regularization, conv=FALSE, iter=0) 
  else { # optimization 
    bothpar <- (!any(is.finite(regularization)))
    n.i <- rep(round(train.fraction*n), fold) # number of observations in bootstrap (training-)samples 
    if (output) {
      cat(" - RDA -\n")
      cat(n, "observations of", p, "variables in", g, "classes,\n")
      if (crossval) cat(paste(fold, "-fold cross-validation.", sep=""),"\n")
      else cat(fold, "bootstrap samples of", n.i[1], "observations each.\n")
      cat("Class names: ", paste(classes[1:(length(classes)-1)], col=",", sep=""),
          classes[length(classes)], "\n")
      if(.Platform$OS.type == "windows") flush.console()
    }  
    # draw bootstrap/crossval samples (row indices): 
    train <- NULL
    test.freq <- NULL
      if (crossval) { # cross-validation 
        indi <- crossval.sample(grouping, fold)
        tabu <- length(indi)-tabulate(indi)
        train <- matrix(0, nrow=max(tabu), ncol=fold)
        for (i in 1:fold) {
          train[1:tabu[i],i] <- which(indi != i)
          test.freq <- cbind(test.freq, tabulate(grouping[-train[,i]], g))
        }
        n.i <- apply(train, 2, function(x) sum(x > 0))
      }
    else { # no cross-validation, but bootstrapping 
      for (i in 1:fold) {
        new <- NULL
        for (j in 1:g) new <- c(new, sample(which(grouping == j), 2))
        new <- c(new, sample((1:n)[-new], n.i - 2 * g))
        train <- cbind(train, new)
        test.freq <- cbind(test.freq, tabulate(grouping[-train[,i]], g))
        # each sample now contains at least 2 elements from each group. 
        # (samples = columns) 
      }
      remove("new")
    }
    # compute parameter estimates (mu & Sigma) for each group in each training sample: 
    means <- covariances <- covpooled <- NULL
    for (i in 1:fold) {
      means <- array(c(means, 
                       as.vector(t(as.matrix(aggregate(data[train[,i],], 
                            by = list(grouping[train[,i]]), mean)[,-1])))),
                     c(p, g, i))
      # (p x g x i)-Array ... each "slice" contains g group means as column vectors. 
      new.covar <- array(unlist(by(data[train[,i],], grouping[train[,i]], var)), c(p,p,g))
      # (p x p x g)-Array, each slice is covariance matrix of one group. 
      covariances <- array(c(covariances,new.covar), c(p,p,g,i))
      # (p x p x g x i)-Array, each "hyperslice" contains g covariance matrices, as above. 
      new.cp <- array(new.covar, c(p*p,g))
      weights <- (tabulate(grouping[train[,i]])-1) / (n.i[i]-g) 
      # weights proportional to fraction of group in sample 
      covpooled <- array(c(covpooled, matrix(new.cp %*% weights, p, p)), c(p, p, i))
      # (p x p x i)-Array, each slice contains pooled covariance for i-th bootstrap sample. 
    }
    remove(list=c("new.covar","new.cp","weights"))
    if (bothpar) { # optimization over both parameters 
      if (is.null(startsimplex)) {
        #startsimplex <- matrix(rbeta(6,1,1),ncol=2) # Beta-RVs for Startsimplex
        startsimplex <- cbind(c(runif(1,1/11,4/11),runif(1,4/11,7/11),runif(1,7/11,10/11)),
                              c(runif(1,1/11,4/11),runif(1,4/11,7/11),runif(1,7/11,10/11)))
        perm <- cbind(c(1,3,2), c(2,1,3), c(3,1,2), c(2,3,1))
        startsimplex[,2] <- startsimplex[perm[,sample(4,1)],2]
      }
      if (trafo) {  # use transformation in Nelder-Mead. 
        linkfunc <- function(x)
            return(1/(1+exp(-x))) # sigmoidal transformation function (for both parameters) 
        linkinverse <- function(x)
            return(-log(1/x-1))   # inverse of link function 
        #startsimplex <- matrix(rnorm(6,0,1),ncol=2) # Normal-RVs for Startsimplex 
        startsimplex <- linkinverse(startsimplex) # transform startsimplex
        mini <- c(-Inf, -Inf)
        maxi <- c(Inf, Inf)
      }
      else {        # do not use transformation. 
        linkfunc <- function(x)
            return(x) # identity function 
        #linkinverse <- linkfunc
        mini=c(0,0)
        maxi=c(1,1)
      }
      dimnames(startsimplex) <- list(NULL, c("gamma","lambda"))
      # NELDER-MEAD #
      print(startsimplex[1,])
      minimize <- optim(startsimplex[1,], goalfunc, method = "Nelder-Mead")
      opti <- list(minimum=linkfunc(minimize$par), 
                     value=minimize$value, conv=minimize$convergence, iter=-1)
     
      
#      nelder.mead(goalfunc, startsimplex, mini=mini, maxi=maxi, 
#        link=linkfunc, max.iter=max.iter, best.possible=0, out=output, 
#        simAnn=simAnn, schedule=schedule, T.start=T.start, 
#        halflife=halflife, zero.temp=zero.temp, alpha=alpha, K=K)
        
        
        
        
    }
    else { # optimization over single parameter 
      logit <- function(x)  return(1/(1+exp(-x)))
      tryval <- logit(seq(-4,4,le=12)) # values to try first 
      if (is.na(regularization[1])) {
        if (output) {
          cat("Optimizing gamma...\n")
          if(.Platform$OS.type == "windows") flush.console()
        }
        goalfu2 <- function(x)
            return(goalfunc(c(x, regularization[2])))
        err <- apply(matrix(tryval, ncol=1), 1, goalfu2)
        fromto <- c(0, tryval, 1)[which.min(err) + c(0,2)]
        minimize <- optimize(goalfu2, fromto)
        opti <- list(minimum=c(minimize$minimum, regularization[2]), 
                     value=minimize$objective, conv=TRUE, iter=-1)
      }
      else {
        if (output) {
          cat("Optimizing lambda...\n")
          if(.Platform$OS.type == "windows") flush.console()
        }
        goalfu2 <- function(x)
        {return(goalfunc(c(regularization[1], x)))}
        err <- apply(matrix(tryval,ncol=1),1,goalfu2)
        fromto <- c(0,tryval,1)[which.min(err)+c(0,2)]
        minimize <- optimize(goalfu2,fromto)
        opti <- list(minimum=c(regularization[1], minimize$minimum), 
                     value=minimize$objective, conv=TRUE, iter=-2)     
      }
    }
  }
  opt.par <- opti$minimum; names(opt.par) <- c("gamma","lambda")
  if (output) {
    cat("Regularization parameters:\n gamma:", round(opt.par[1],5), 
        "  lambda:", round(opt.par[2],5), "\n")
    if(.Platform$OS.type == "windows") flush.console()
  }
  # compute parameters for complete data: 
  means <- t(as.matrix(aggregate(data,by=list(grouping),mean)[,-1]))
  dimnames(means) <- list(varnames,classes)
 # covariances <- array(unlist(by(data[train[,i],],grouping[train[,i]],var)),
 #                      c(p,p,g), dimnames=list(varnames,varnames,classes))
  covariances <- array(unlist(by(data,grouping,var)),
                       c(p,p,g), dimnames=list(varnames,varnames,classes))
 # weights <- (tabulate(grouping)-1) / (n.i-g)
  weights <- (tabulate(grouping)-1) / (n-g)
  covpooled <- matrix(array(covariances, c(p*p,g)) %*% weights, p, p,
                      dimnames = list(varnames, varnames))
  # predict training data, compute apparent error rate: 
  if (estimate.error) {
    errors <- classify(data, means, covariances, covpooled, opt.par[1], opt.par[2], g=g, p=p, n=n, prior=prior) != grouping
    group.rates <- tabulate(grouping[errors],g) / tabulate(grouping)
    APER <- as.vector(t(prior) %*% group.rates)
    if (output) 
        cat("Apparent error rate (APER) for training data:", 
            round(APER * 100, 3), "%\n")  
  }
  else APER <- NA
  if (crossval) err <- c("APER"=APER, "crossval"=opti$value)
  else err <- c("APER"=APER, "bootstrap"=opti$value)
  result <- list(call=match.call(), 
                 regularization=opt.par, classes=classes, prior=prior, error.rate=err,
                 varnames=varnames,
                 means=means, covariances=covariances, covpooled=covpooled, 
                 converged=opti$conv, iter=opti$iter)
  class(result) <- "rda"
  return(result)
}

predict.rda <- function(object, newdata, posterior=FALSE, aslist=FALSE)
{
  classify <- function(dataset, mu, sigma, pooled, gamma, lambda, g, p, n, prior)
  # difference to `classify'-function above is that LIKELIHOODS are returned instead of classifications
  # mu            : matrix with group means as columns 
  # sigma         : (p x p x g)-Array of group covariances 
  # pooled        : pooled covariance matrix 
  # gamma, lambda : regularization parameters
  # g, p, n       : numbers of classes, variables & observations
  {
    # compute likelihood of each datum for each group: 
    likelihood <- matrix(nrow=n, ncol=g)
    for (i in 1:g){ 
      reg.cov <- (1-lambda) * sigma[,,i] + lambda * pooled               # shift towards pooled cov. 
      reg.cov <- if (is.matrix(reg.cov)) # cov is a matrix
                   (1-gamma)*reg.cov + gamma*mean(diag(reg.cov))*diag(p) # shift towards identity 
                 else # one variable  =>  cov is 1x1-matrix:
                   (1-gamma)*reg.cov + gamma*reg.cov                     # shift towards identity 
      likelihood[,i] <- prior[i] * dmvnorm(dataset, mu = mu[,i], inv.sigma = solve(reg.cov))
    }
    return(likelihood)
  }

  p <- dim(object$means)[1]
  g <- dim(object$means)[2]
  if (!any(is.null(colnames(newdata)),is.null(object$varnames))) { # (both colnames & varnames are given) 
    if(all(is.element(object$varnames,colnames(newdata)))){        # (varnames is a subset of colnames)   
      newdata <- as.matrix(newdata[,object$varnames])
    }
  }
  if(is.vector(newdata)) newdata <- matrix(newdata, ncol=1)
  n <- dim(newdata)[1]
  likeli <- classify(newdata, object$means, object$covariances, object$covpooled,
                     object$regul[1], object$regul[2], g=g, p=p, n=dim(newdata)[1], 
                     prior=object$prior)
  colnames(likeli) <- object$classes
  classi <- apply(likeli, 1, function(x) order(x)[g])
  classi <- factor(object$classes[classi], levels = object$classes)
  postmat <- if (posterior) likeli / rowSums(likeli)
             else NULL
  if (aslist) result <- list("class" = classi, "posterior" = postmat)
  else { 
    result <- classi
    attr(result, "posterior") <- postmat
  }
  return(result)
}


print.rda <- function(x,...)
{
  #cat(" - RDA - \n")
  cat("Call:\n")
  print(x$call)
  cat("\nRegularization parameters:\n")
  #cat("gamma:", round(x$regu[1],5), " lambda:", round(x$regu[2],5), "\n")
  print(x$regu)
  #cat("\nClass prior:\n")
  cat("\nPrior probabilities of groups:\n")
  print(x$prior)
  cat("\nMisclassification rate:\n")
  cat("       apparent:",
    ifelse(is.na(x$error.rate[1]), "--", as.character(round(x$error.rate[1] * 100, 3))), "%\n")
  if (length(x$error.rate) > 1)
    cat(ifelse(names(x$error.rate)[2] == "crossval", 
        "cross-validated:", "   bootstrapped:"), 
        as.character(round(x$error.rate[2] * 100, 3)), "%\n")
  invisible(x)
}


plot.rda <- function(x, textplot=FALSE, ...)
{
  parpty <- par("pty")
  par(pty="s")
  if(textplot) {
    plot(c(0,1),c(0,1), type="n", axes=FALSE, xlab="groups", ylab="covariances",...)
    textcol <- "darkgrey"
    textsize <- 1.5
    textshift <- 0.02
    text(0+textshift, 0+textshift, "QDA", adj=c(0,0), cex=textsize, col=textcol)
    text(1-textshift, 0+textshift, "LDA", adj=c(1,0), cex=textsize, col=textcol)
    #text(0+textshift, 1-textshift, "cond. indep.", adj=c(0,1), cex=textsize, col=textcol)
    #text(1-textshift, 1-textshift, "nearest mean", adj=c(1,1), cex=textsize, col=textcol)
    text(0.5, 1-textshift, "i.i.d. variables", adj=c(0.5,1), cex=textsize, col=textcol)
    axis(1, at=c(0,1), labels=c("unequal", "equal"))
    axis(2, at=c(0,1), labels=c("correlated", "diagonal"))
  }
  else{
    plot(c(0,1),c(0,1), type="n", axes=FALSE, xlab=expression(lambda), ylab=expression(gamma),...)
    axis(1)
    axis(2)
  }
  lines(c(0,1,1,0,0), c(0,0,1,1,0), col="grey")
  lines(c(0,1), rep(x$regu[1],2), col="red1", lty="dotted")
  lines(rep(x$regu[2],2), c(0,1), col="red1", lty="dotted")
  points(x$regu[2], x$regu[1], pch=18, col="red2")
  par(pty=parpty)
  invisible(x$regu)
}


dmvnorm <- function(x, mu=NA, inv.sigma=NA, logDetCov=NA)    
# Density of a Multivariate Normal Distribution 
# works for matrices as well as for vectors     
# (matrices are evaluated row-wise)             
#   !!   supply inverse of covariance   !!      
# `logDetCov' = log(det(Cov)) = log(1/det(inv.Cov)) = -1*log(det(inv.Cov))
{
  if (is.vector(x)) x <- t(x) # x is treated as 1 obsevation of (length(x)) variables
  if (is.na(logDetCov)) logDetCov <- -log(det(inv.sigma))
  singledens <- function(x, M=mu, IS=inv.sigma, ldc=logDetCov)  # density function for a single vector 
  {
    xm <- x - M
    return(as.numeric(exp(- 0.5 * 
           (length(M) * log(2*pi) 
            + ldc
            + (t(xm) %*% IS %*% xm)))))
  } 
  return(apply(x, 1, singledens))#, M=mu, IS=inv.sigma, ldc=logDetCov))
}  
