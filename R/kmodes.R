kmodes <- function(data, modes, iter.max = 10, weighted = FALSE, fast = TRUE) {
    
    if(!is.data.frame(data)) data <- as.data.frame(data)  
    # check data for coloumn types
    isnumeric <- sapply(data, is.numeric)
    isfactor <- sapply(data, is.factor)
    if(any(isfactor)){
      levs <- vector("list", ncol(data))
      for(j in which(isfactor)){ 
        levsj <- levels(data[,j]) 
        data[,j] <- levsj[data[,j]] # replace factors by characters (for computation of distances if weighted == TRUE)
        levs[[j]] <- levsj # store levels for back transform of modes to factor (otherwise some levels might get lost)
      }
    }
      
    # check for numerics where the values potentially do not represent categories 
    if(any(isnumeric)) {
      lengths <- sapply(data[,isnumeric], function(z) return(length(unique(z))))
      if(any(lengths > 30)) warning("data has numeric coloumns with more than 30 different levels!")
    }
    
    ### updates the mode of cluster "num",
    ### is called everytime an object switches clusters
    update_mode <- function(num, num_var, data, cluster) {
        ## gather all objects of cluster "num"
        clust <- data[which(cluster == num),] 
        apply(clust, 2, function(cat) {
            ## compute the most frequent category for each variable
            cat <- table(cat)
            names(cat)[which.max(cat)]
        })
    } # end update mode


    ### computes the weighted distance between an object "obj" and the mode of a cluster
    ### the frequencies of categories in "data" are used for the weighting
    distance <- function(mode, obj, weights){
        if(is.null(weights))
            return(sum(mode != obj))
            
        #obj <- as.integer(obj) # may lead to wrong calc. of ln.51/52 (no match)
        obj <- as.character(obj)
        mode <- as.character(mode)
        different <- which(mode != obj) 
        n_mode <- n_obj <- numeric(length(different))
        for (i in seq(along = different)) { 
            ## frequencies are computed only if neccessary (not if distance is 0 anyway)
            weight <- weights[[different[i]]]
            names <- names(weight)
            n_mode[i] <- weight[which(names==mode[different[i]])]
            n_obj[i] <- weight[which(names==obj[different[i]])]
        }
        dist <- sum((n_mode + n_obj) / (n_mode * n_obj))
        return(dist)
    } # end distance


    n <- nrow(data)
    num_var <- ncol(data)
    data <- as.data.frame(data)
    
    cluster <- numeric(n) # the actual partition
    names(cluster) <- 1:n
    
    if (missing(modes)) 
        stop("'modes' must be a number or a data frame.")
    if (iter.max < 1) 
        stop("'iter.max' must be positive.")

    if (length(modes) == 1) { 
        ## if the desired number "k" of clusters is given, choose k object randomly from "data"
        k <- modes
        modes <- unique(data)[sample(nrow(unique(data)))[1:k],] 
        for (i in 1:k)
            cluster[which(rownames(data) == rownames(modes)[i])] <- i
    } else { 
        ## use given modes if appropriate
        if (any(duplicated(modes))) 
            stop("Initial modes are not distinct.")
        if (ncol(data) != ncol(modes)) 
            stop("'data' and 'modes' must have same number of columns")
        #modes <- as.matrix(modes) # replaced
        modes <- as.data.frame(modes)
        if(any(isfactor)){ # here further checks on could be added...
          if(!all(sapply(modes[,isfactor], is.factor))) stop("Types of modes do not match data!")
          for(j in which(isfactor)) modes[,j] <- levels(modes[,j])[modes[,j]]
        }
        for(j in 1:num_var) if(weighted) if(!all(modes[,j] %in% unique(data[,j]))) stop("For weighted call values of modes must exist in data!")
        k <- nrow(modes)
    }
    if (k > nrow(unique(data))) 
        stop("More cluster modes than distinct data points.")

    if(weighted){
        ## compute the frequencies of each category for each variable
        weights <- vector("list", num_var)
        for (i in 1:num_var) 
            weights[[i]] <- table(data[,i])
    } else {
        weights <- NULL
    }
    
    # original implementation
    if(!fast){
      for (j in which(cluster==0)) { 
          ## first put all not yet assigned objects into the cluster, which has the nearest mode
          dist <- apply(modes, 1, function(x) distance(x, data[j,], weights))
          cluster[j] <- which.min(dist)
          modes[cluster[j],] <- update_mode(cluster[j], num_var, data, cluster) 
          ## update the mode of the cluster the object has been assigned to
      }
      for (i in 1:iter.max) {
          continue <- FALSE
          for (j in 1:n) { 
              ## run through all objects and assign them to the cluster with the nearest mode
              ## (or leave everything as is, if the object's current cluster's mode is still the nearest)
              dist <- apply(modes, 1, function(x) distance(x, data[j,], weights))
              clust_new <- which.min(dist)
              clust_old <- cluster[j]
              if (clust_new != clust_old) { 
                  ## update the modes of old and new cluster, if object has switched from "clust_old" to "clust_new"
                  cluster[j] <- clust_new
                  modes[clust_new,] <- update_mode(clust_new, num_var, data, cluster)
                  modes[clust_old,] <- update_mode(clust_old, num_var, data, cluster)
                  continue <- TRUE
              }
          }
          ## break, if during a complete iteration no object has changed clusters
          if (!continue) break
      } # end iterations
    }

    ### fast version of the algorithm 
    ### ...running through the entire data set at once before updating clusters which allows for apply
    if(fast){
      ## first put all not yet assigned objects into the cluster, which has the nearest mode
      # compute distances: either weighted or unweighted
      dists <- matrix(NA, nrow=n, ncol = k)
      if(!weighted){
        for(i in 1:k){
          di <- sapply(1:ncol(data), function(j) return(data[,j] != rep(modes[i,j], n)) )
          di <- rowSums(di)
          dists[,i] <- di
        }
      }
      if(weighted){
        # compute frequencies (for weighted computation only once necessary)
        n_obj <- matrix(NA, nrow=n, ncol = ncol(data))
        for(j in 1:ncol(data)) n_obj[,j]  <- weights[[j]][sapply(as.character(data[,j]), function(z) return(which(names(weights[[j]])==z)))]    

        # necessary for updated modes after each iteration
        n_mode <- matrix(NA, nrow=nrow(modes), ncol = ncol(data))
        for(j in 1:ncol(data)) n_mode[,j]  <- weights[[j]][sapply(as.character(modes[,j]), function(z) return(which(names(weights[[j]])==z)))]          
        for(i in 1:k){
          di <- sapply(1:ncol(data), function(j) return(data[,j] != rep(modes[i,j], n)) )
          wts <- (n_mode[rep(i,n),] + n_obj) / (n_mode[rep(i,n),] * n_obj)
          di <- rowSums(di*wts)
          dists[,i] <- di
        }
      }
      # assign clusters 
      #old.cluster  <- cluster
      cluster      <- apply(dists, 1, function(z) {a <- which.min(z); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
      #min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])      
      # update modes
      for(j in 1:nrow(modes)) modes[j,] <- update_mode(j, num_var, data, cluster)        
      
      # start iterations
      for (i in 1:iter.max) {
        continue <- FALSE
        # compute distances: either weighted or unweighted
        dists <- matrix(NA, nrow=n, ncol = k)
        if(!weighted){
          for(i in 1:k){
            di <- sapply(1:ncol(data), function(j) return(data[,j] != rep(modes[i,j], n)) )
            di <- rowSums(di)
            dists[,i] <- di
          }
        }
        if(weighted){
          # necessary for updated modes after each iteration
          n_mode <- matrix(NA, nrow=nrow(modes), ncol = ncol(data))
          for(j in 1:ncol(data)) n_mode[,j]  <- weights[[j]][sapply(as.character(modes[,j]), function(z) return(which(names(weights[[j]])==z)))]          
          for(i in 1:k){
            di <- sapply(1:ncol(data), function(j) return(data[,j] != rep(modes[i,j], n)) )
            wts <- (n_mode[rep(i,n),] + n_obj) / (n_mode[rep(i,n),] * n_obj)
            di <- rowSums(di*wts)
            dists[,i] <- di
          }
        }
        # assign clusters 
        old.cluster  <- cluster
        cluster      <- apply(dists, 1, function(z) {a <- which.min(z); if (length(a)>1) a <- sample(a,1); return(a)}) # sample in case of multiple minima
        #min.dists     <- apply(cbind(clusters, dists), 1, function(z) z[z[1]+1])      
        # update modes
        for(j in 1:nrow(modes)) modes[j,] <- update_mode(j, num_var, data, cluster)        
        
        ## break, if during a complete iteration no object has changed clusters
        if(any(old.cluster != cluster)) continue <- TRUE
        if (!continue) break
      } # end iterations
    }

    # compute results
    cluster.size <- table(cluster) # sizes of all clusters
    if (length(cluster.size) < k) 
        warning("One or more clusters are empty.")
      
    # compute final matrix of distances: either weighted or unweighted
    dists <- matrix(NA, nrow=n, ncol = k)
    if(weighted){
      n_mode <- matrix(NA, nrow=nrow(modes), ncol = ncol(data))
      for(j in 1:ncol(data)) n_mode[,j]  <- weights[[j]][sapply(as.character(modes[,j]), function(z) return(which(names(weights[[j]])==z)))]          
    }
    for(i in 1:k){
      di <- sapply(1:ncol(data), function(j) return(data[,j] != rep(modes[i,j], n)) )
      if(!weighted) di <- rowSums(di)        
      if(weighted){
        if(!fast){# ...for fast == FALSE n_obj ist computed within calls of distance()
          n_obj <- matrix(NA, nrow=n, ncol = ncol(data))
          for(j in 1:ncol(data)) n_obj[,j]  <- weights[[j]][sapply(as.character(data[,j]), function(z) return(which(names(weights[[j]])==z)))]    
        }
        wts <- (n_mode[rep(i,n),] + n_obj) / (n_mode[rep(i,n),] * n_obj)
        di <- rowSums(di*wts)        
      }
      dists[,i] <- di
    }
    
    ## compute for each cluster the sum of the distances of all objects of a cluster to the cluster's mode:
    diffs <- numeric(k)
    for (i in seq_along(cluster.size)) diffs[i] <- sum(dists[cluster == i,i])
      #diffs[i] <- sum(apply(data[cluster == i,], 1, function(x) sum(x != modes[i,]) ))
    
    rownames(modes) <- 1:k
    colnames(modes) <- colnames(data)
    
    # convert coloumn modes to original coloumn types of data
    if(any(isfactor))  for(j in which(isfactor)) modes[,j] <- factor(modes[,j], levels = levs[[j]])
    if(any(isnumeric)) for(j in which(isnumeric)) modes[,j] <- as.numeric(modes[,j])
        
    result <- list(cluster = cluster, size = cluster.size, modes = modes, 
        withindiff = diffs, iterations = i, weighted = weighted)
    class(result) <- "kmodes"
    return(result)
}


print.kmodes <- function(x, ...)
{
    cat("K-modes clustering with ", length(x$size), " clusters of sizes ",
        paste(x$size, collapse=", "), "\n", sep="")
    cat("\nCluster modes:\n")
    print(x$modes, ...)
    cat("\nClustering vector:\n")
    print(x$cluster, ...)
    cat("\nWithin cluster simple-matching distance by cluster:\n")
    print(x$withindiff, ...)
    cat("\nAvailable components:\n")
    print(names(x))
    invisible(x)
}
