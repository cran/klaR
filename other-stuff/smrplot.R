smrplot <- function(posterior, trueclass = NULL, center = FALSE, 
                col = rainbow(ncol(posterior)), pchest = rep(21, ncol(posterior)),
                pchwrong = rep(20, ncol(posterior)), pchclass = rep(19, ncol(posterior)), ...)
{
    ncp <- ncol(posterior)
    if (ncp == 4) 
        require(scatterplot3d)
    else if (ncp != 3)
        stop(paste(ncp, "-classes not supported\n"))
    
    estclass <- factor(max.col(posterior), levels = seq(along = colnames(posterior)), 
        labels = colnames(posterior))

    if (!is.null(trueclass)) 
        wrongest <- trueclass != estclass

  simpprob <- function(x)
  {
    if(is.vector(x)) x <- t(x)
    if (ncol(x) == 3)
        return(x %*% matrix(c(1, 0.5, 0, 0, sqrt(3)/2, 0), nrow = 3))
    else if (ncol(x)==4)
        return(x %*% 
            matrix(c(1, 0.5, 0.5, 0, 
                0, sqrt(3)/2, sqrt(3)/6, 0, 
                0, 0, sqrt(6)/3, 0), nrow = 4))
  }

  baurand <- function(dimi = 4, center = FALSE)
  {
    punkte <- matrix(0, nrow = 2*sum(1:(dimi-1)), ncol = dimi)
    schwerpunkt <- rep(1/dimi, dimi)
    k <- 1
    for (i in 1:(dimi-1))
        for (j in (i+1):dimi)
            {
                punkte[k, i] <- 1
                punkte[(k+1), j] <- 1
                k <- k + 2
            }
    if (center)
    {
        for (i in 1:(dimi-1))
            dummy1 <- numeric(dimi)
            dummy1[i] <- 0.5
            for (j in i:dimi)
            {
                dummy <- dummy1
                dummy[j] <- 0.5
                punkte <- rbind(punkte, dummy, schwerpunkt)
            }
    }    
    return(punkte)
  }

  if (ncol(posterior) == 4)
    {
        s3d <- scatterplot3d(simpprob(baurand(4, center = center)), 
            type = "l", axis = FALSE, grid = FALSE, angle = 75, scale.y = 0.5, ...)
        s3d$points3d(simpprob(diag(4)), col = col, pch = pchclass, cex = 1.5)
        s3d$points3d(simpprob(posterior), col = col[estclass], 
            pch = pchest[estclass], cex = 1)
        if(!is.null(trueclass)) 
            s3d$points3d(simpprob(posterior[wrongest, ]), col = col[trueclass[wrongest]], 
                pch = pchwrong[trueclass[wrongest]], cex = 1)
        legend(s3d$xyz.convert(0, 0, 1), legend = colnames(posterior), 
            col = col, pch = pchclass, cex = 1.2)
    }
    else if(ncol(posterior) == 3)
    {
        s2d <- plot(simpprob(baurand(3, center = center)),  
            type = "l", axes = FALSE, ann = FALSE, ...)
        points(simpprob(diag(3)), col = col, pch = pchclass, cex = 1.5)
        points(simpprob(posterior), col = col[estclass], pch = pchest[estclass], cex = 1)
        if (!is.null(trueclass)) 
            points(simpprob(posterior[wrongest, ]), 
                col = col[trueclass[wrongest]], pch = pchwrong[trueclass[wrongest]], cex = 1)
        legend(0.7, 0.8, legend = colnames(posterior), col = col, 
            pch = pchclass, cex = 1.2)
    }
}
