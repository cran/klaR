predict.NaiveBayes<-function (object, newdata, 
    threshold = 0.001, ...) 
{
    if (missing(newdata)) newdata<-object$x
    if (!any(is.null(colnames(newdata)),is.null(object$varnames))) { # (both colnames & varnames are given) 
        if(all(is.element(object$varnames,colnames(newdata)))){        # (varnames is a subset of colnames)   
            newdata <- data.frame(newdata[,object$varnames])
        }
    }
    
    nattribs <- ncol(newdata)
    isnumeric <- sapply(newdata, is.numeric)
    L <- sapply(1:nrow(newdata), function(i) {
        ndata <- as.numeric(newdata[i, ])
        L <- object$apriori * apply(sapply(1:nattribs, function(v) {
            nd <- ndata[v]
            if (is.na(nd)) 
                rep(1, length(object$apriori))
            else {
                prob <- if (isnumeric[v]) {
                  msd <- object$tables[[v]]
                  if (object$usekernel) sapply(msd,FUN=function(y){dkernel(x=nd,kernel=y,...)})
                  else dnorm(nd, msd[, 1], msd[, 2])
                }
                else object$tables[[v]][, nd]
                prob[prob == 0] <- threshold
                prob
            }
        }), 1, prod)
        L/sum(L)
    })
    posterior<-t(L)
    classdach<-factor(object$levels[apply(L, 2, which.max)], levels = object$levels)
    colnames(posterior)<-object$levels
    rownames(posterior)<-names(classdach)<-rownames(newdata)
    return(list(class=classdach,posterior=posterior))
}
