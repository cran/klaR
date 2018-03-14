corclust <- function(x, cl = NULL, method = "complete"){ 
    
    # x         : Data frame or matrix of variables to be grouped
    # cl        : optional vector containing class labels
    # mincor    : minimum required correlation between the variables of a cluster (a line is drawn)
    # prnt      : should distance matrix of absolute correlations be printed
    # method    : hierarchical clustering algorithm to be used (single, average, ward, default is complete linkage)

    #if(is.matrix(x)) if(!is.numeric(x)) stop("Attributes of x must be numerical!")
    #if(is.data.frame(x)) if(!all(sapply(x,is.numeric))) stop("Attributes of x must be numerical!")
    if(!is.null(cl)) if(!is.factor(cl))  stop("If classes are defined cl must be a factor!") 
    if(length(dim(x)) < 2)  stop("x must be either a matrix or a data frame of numeric attributes!")
    if(min(dim(x))==1)  stop("x must contain more than one variable and more than one observation!")
    if(!is.null(cl)){ 
        if(min(table(cl)) < 2) stop("At least one class consist of less than 2 observations!")
        if(min(table(cl)) < dim(x)[2]) warning("At least one of the classes consist of less observations than variables!")
        }
    
    nx <- colnames(x)
    x <- as.data.frame(x)
    names(x) <- nx
    
    numids <- which(sapply(x, is.numeric))
    facids <- which(sapply(x, is.factor))

    if(length(numids)==0){cat("No numeric variables!\n\n"); numres <- NULL; kor <- NULL}
    if(length(numids)==1){cat("Only one numeric variable. No cluster model for numerics returned!\n\n"); numres <- NULL; kor <- NULL}
    if(length(numids)>1){
      if(is.null(cl)){
        kor <- 1-abs(cor(x[,numids],use="pairwise.complete.obs")) # not yet the matrix of distances!
      }
      else{
        kor <- matrix(1,ncol(x[,numids]),ncol(x[,numids]))
        for(k in seq(along=levels(cl))){
          cls <- levels(cl)[k]
          oldkor <- kor
          kor <- 1-abs(cor(x[which(cl==cls),numids],use="pairwise.complete.obs")) # not yet the matrix of distances!
          newkor <- (oldkor - kor) > 0
          kor[newkor] <- oldkor[newkor]  
        }    
      }
      
      distances <- NULL
      for (i in 1:(ncol(kor)-1)) distances <- c(distances,kor[(i+1):nrow(kor),i]) # only lower triangular matrix, without diagonal elements
      attr(distances,"Labels") <- colnames(x[,numids])
      attr(distances,"Size") <- ncol(x[,numids])
      attr(distances,"Metric") <- "absolute correlations"
      class(distances) <- "dissimilarity" 
      
      #if(prnt) print(kor)
      #variables <- distances
      numres <- hclust(distances, method=method)
        
    }  
    
    if(length(facids)==0) {cat("No factor variables!\n\n"); facres <- NULL; crv <- NULL}
    if(length(facids)==1){cat("Only one factor variable. No cluster model for factors returned!\n\n"); facres <- NULL; crv <- NULL}
    if(length(facids)>1){
      N <- nrow(x)
      
      crv <- function(X1,X2,n=N){ # fucntion to calculate cramer's V statistic
        l <- nlevels(X1)
        m <- nlevels(X2)
        chi2 <- chisq.test(X1,X2,correct=F)
        V <- sqrt(chi2$statistic/(n*min(l-1,m-1)))
        return(V)
      }
      
      if(is.null(cl)){
        distances <- matrix(0, nrow=ncol(x[,facids]), ncol=ncol(x[,facids]))
        rownames(distances) <- colnames(x[,facids])
        colnames(distances) <- colnames(x[,facids])
        for (i in 1:(nrow(distances)-1)){
          for (j in (i+1):ncol(distances)){
            distances[i,j] <- distances[j,i] <- 1-crv(x[,facids[i]],x[,facids[j]])
          }
        }
        # for (i in 1:(nrow(distances))){
        #   distances[i,] <- sapply(1:nrow(distances), function(j) return(crv(x[,facids[i]],x[,facids[j]])))
        # } # ...no speedup!
      }
      else{
        distances <- matrix(0,nrow=ncol(x[,facids]),ncol=ncol(x[,facids]))
        rownames(sim) <- colnames(x[,facids])
        colnames(sim) <- colnames(x[,facids])
        for(k in seq(along=levels(cl))){
          cls <- levels(cl)[k]
          olddist <- distances
          for (i in 1:(nrow(distances)-1)){
            for (j in (i+1):ncol(distances)){
              distances[i,j] <- distances[j,i] <- 1-crv(x[which(cl==cls),facids[i]],x[which(cl==cls),facids[j]])
            }
          }
          newdist <- (olddist - distances) > 0
          distances[newdist] <- olddist[newdist]  
        }
      }
      
      crv <- 1-distances
      distances <- as.dist(distances)
      facres <- hclust(distances, method=method)
    }  
    
    
    result <- list(cor=1-kor, crv=crv, cluster.numerics=numres, cluster.factors=facres, id.numerics=numids, id.factors=facids)
    class(result) <- "corclust"
    return(result)
    }


plot.corclust <- function(x, selection = "both",  mincor = NULL, ...){

  if(is.null(x$cluster.numerics)) selection <- "factor"
  if(is.null(x$cluster.factors)) selection <- "numeric"
  if(selection == "both"){par(ask=TRUE)}
  
  if(selection == "both" | selection == "numeric"){
  ylb <- "1 - absolute correlation within cluster"
  if(x$cluster.numerics$method == "complete") ylb <- "1 - minimum absolute correlation within cluster"
  if(x$cluster.numerics$method == "average") ylb <- "1 - average absolute correlation within cluster"
  if(x$cluster.numerics$method == "single") ylb <- "1 - maximum absolute correlation within cluster"
  plot(x$cluster.numerics, ylab = ylb) # , main=paste("1 - absolute correlations between variables")
  if(!is.null(mincor)) abline(h=1-mincor, col=2)
  }

  if(selection == "both" | selection == "factor"){
    ylb <- "1 - Cramer's V within cluster"
    if(x$cluster.factors$method == "complete") ylb <- "1 - minimum Cramer's V within within cluster"
    if(x$cluster.factors$method == "average") ylb <- "1 - average Cramer's V within within cluster"
    if(x$cluster.factors$method == "single") ylb <- "1 - maximum Cramer's V within within cluster"
    plot(x$cluster.factors, ylab = ylb) # , main=paste("1 - absolute correlations between variables")
    if(!is.null(mincor)) abline(h=1-mincor, col=2)
  }
  par(ask=FALSE)
}

### function that extract cluster IDs for variables from object fo class corclust 
cvtree <- function(object, k = 2, mincor = NULL, ...){
  
  if(class(object) != "corclust") stop("Object must be of class corclust!")
  numfac <- !c(is.null(object$cluster.numerics),is.null(object$cluster.factors))
  if(!is.null(mincor)){
    k <- NULL
    h <- 1-mincor
    if(length(h) == 1) h <- rep(h,2)
    }
  if(is.null(k)){if(is.null(h)) stop("If (is.null(k)) mincor must not be NULL!")}
  if(length(k) == 1) k <- rep(k,2)
  if(!is.null(k)){
    if(is.null(object$cluster.numerics)) k[1] <- 0
    if(is.null(object$cluster.factors)) k[2] <- 0
  }
    
  # extract clusters of numeric variables
  knum <- 0
  if(!is.null(object$cluster.numerics)){
    if(!is.null(k)){
      knum <- k[1] 
      cluster.num <- cutree(object$cluster.numerics, knum)
    }
    else{
      hnum <- h[1]
      cluster.num <- cutree(object$cluster.numerics, h = hnum)
      knum <- length(table(cluster.num))
    }
    # compute summary of average correlation to cluster / closest other cluster forall vars.
    if(knum > 1){
      cors <- object$cor 
      diag(cors) <- NA 
      # compute average correlation of the variables to all clusters
      avcor2cl <- function(clid){apply(cors[,which(cluster.num==clid),drop=FALSE], 1, mean, na.rm=TRUE)}
      cors2cls <- sapply(1:knum, avcor2cl)
      within <- sapply(1:nrow(cors2cls), function(i) cors2cls[i,cluster.num[i]])# ...to own clusters
      # check for single variable clusters
      tab <- table(cluster.num)
      if(any(tab==1)){singles <- as.integer(names(tab)[tab == 1]); within[cluster.num %in% singles] <- 1} 
      # ...to other clusters 
      rest <- cors2cls
      rest <- t(sapply(1:nrow(cors2cls), function(i){a<-cors2cls[i,]; a[cluster.num[i]]<-0; return(a)})) 
      rownames(rest) <- rownames(cors)
      closest <- apply(rest, 1, which.max)
      cor2closest <- apply(rest, 1, max)
      comparison <- cbind(av.within.cor=within, av.cor2closest=cor2closest, closest)
      
      idspercl <- by(comparison[,1], cluster.num, function(x) return(sort(x, decreasing=TRUE, index.return= TRUE)$ix))
      comparison <- cbind(cluster = cluster.num, comparison)
      comparison <- comparison[sort(cluster.num, index.return=TRUE)$ix,]
      for(i in 1:knum){
        a <- comparison[which(comparison[,1]==i),,drop = F]
        a <- a[idspercl[[i]],]
        comparison[which(comparison[,1]==i),] <- a
      }
      comparison.num  <- comparison
    } else comparison.num <- NULL
    
  } else {cluster.num <- NULL; comparison.num <- NULL}
  
  # extract clusters of factor variables  
  if(!is.null(object$cluster.factors)){
    if(!is.null(k)){
      kfac <- k[2] 
      cluster.fac <- knum + cutree(object$cluster.factors, kfac)
    }
    else{
      hfac <- h[2]
      cluster.fac <- knum + cutree(object$cluster.factors, h = hfac)
      kfac <- length(table(cluster.fac))
    }
    # compute summary of average correlation to cluster / closest other cluster forall vars.
    if(kfac>1){
      crv <- object$crv 
      diag(crv) <- NA 
      avcor2cl <- function(clid){apply(crv[,which(cluster.fac==(clid+knum)), drop=FALSE], 1, mean, na.rm=TRUE)}
      # compute average CrVs to all clusters
      cors2cls <- sapply(1:kfac, avcor2cl)
      within <- sapply(1:nrow(cors2cls), function(i) cors2cls[i,(cluster.fac[i]-knum)]) # ...to own clusters
      # check for single variable clusters
      tab <- table(cluster.fac)
      if(any(tab==1)){singles <- as.integer(names(tab)[tab == 1]); within[cluster.fac %in% singles] <- 1} 
      # ...to other clusters
      rest <- cors2cls
      rest <- t(sapply(1:nrow(cors2cls), function(i){a <- cors2cls[i,]; a[(cluster.fac[i]-knum)] <- 0; return(a)})) 
      rownames(rest) <- rownames(crv)
      closest <- apply(rest, 1, which.max) + knum
      cor2closest <- apply(rest, 1, max)
      comparison <- cbind(av.within.cor=within, av.cor2closest=cor2closest, closest)
      
      idspercl <- by(comparison[,1], cluster.fac, function(x) return(sort(x, decreasing=TRUE, index.return= TRUE)$ix))
      comparison <- cbind(cluster = cluster.fac, comparison)
      comparison <- comparison[sort(cluster.fac, index.return=TRUE)$ix,]
      for(i in 1:kfac){
        a <- comparison[which(comparison[,1]==(knum+i)),,drop = F]
        a <- a[idspercl[[i]],]
        comparison[which(comparison[,1]==(knum+i)),] <- a   
      }
      comparison.fac <- comparison
    } else comparison.fac <- NULL
  } else {cluster.fac <- NULL; comparison.fac <- NULL}
  
  clusters <- c(cluster.num, cluster.fac)
  reorder  <- sort(c(object$id.numerics, object$id.factors), index.return=TRUE)$ix
  clusters <- clusters[reorder]
  
  comparison <- rbind(comparison.num, comparison.fac)
  res <- list("cluster" = clusters, "correlations" = comparison)
  class(res) <- "cvtree"
  return(res)
}



# convenience function to apply variable selction based on an obect of class cvtree (or the summary of a CLV object)
xtractvars <- function(object, data, thres = 0.5){
  if(class(object) != "cvtree") {if(!any("groups" %in% names(object))) stop("Object must be of class cvtree!")}
                                 # not only for application on cvtree output but also on summary of CLV objects          
  
  if(class(object) == "cvtree"){
    clids <- unique(object$correlations[,1])
    vsel <- NULL
    for(i in clids){
      clus <- object$correlations[which(object$correlations[,1]==i), , drop=F]
      vsel <- c(vsel, rownames(clus)[unique(c(1, which(abs(clus[,2]) < thres)))])
    }
  }  
  else{
    clids <- 1:length(object$groups)
    vsel <- NULL
    for(i in clids){
      vsel <- c(vsel, rownames(object$groups[[i]])[unique(c(1,which(abs(object$groups[[i]][,1]) < thres)))])
    }
  }
  
  return(data[,which(colnames(data) %in% vsel)])
}




# sonderfall: nur 1 numeric / factor
# in clv: k=1 oder mincor, so dass k=1! ~> abfangen
# suggest ISLR
