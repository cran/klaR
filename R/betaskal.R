b.scal <- function(member, grouping, dis = FALSE, eps = 0.0001)
{
betaregion <- function(betaobj, lev, dis = FALSE, eps = 0.0001)
    {
    member <- data.matrix(betaobj[,1:length(lev)])
    colnames(member) <- lev
    kdach <- lev[max.col(member)]
    CR <- mean(kdach == betaobj$grouping)
    k <- which(kdach[1] == lev)
    memvec <- member[ , k]
    NT <- length(memvec)
    pm <- mean(memvec)
    S <- var(memvec)  
    if(S < eps)
        {
           NM <- min((pm * (1-pm) / S) - 1, 1e4)
           N <- floor(min(NT, NM))
        }
    else 
        {
        NM <- (pm * (1-pm) / S) - 1
        NU <- floor(min(NT, NM))
        NO <- floor(max(NT, NM))
        if (!dis) N <- NU 
        else # choose optimal N
            {
            N <- NU:NO
            DC <- numeric(length(N))
            for (i in seq(along = N))
                {
                memvecneu <- qbeta(pbeta(memvec, pm * NM, (1-pm) * NM), CR * N[i] , (1-CR) * N[i]) 
                member[,-k] <- ifelse(rep(memvec, length(lev) - 1) > (1 - 1e-6), 
                                        (1 - memvecneu) / (length(lev) - 1), 
                                        member[ , -k] * (1 - memvecneu) / (1 - memvec))
                member[,k] <- memvecneu
                ergi <- ucpm(member, betaobj$grouping)
                DC[i] <- ergi$CR * NT + ergi$AC ### optimize CR+AC!!!!
                }
            N <- N[which.max(DC)]
            }
        }
    return(list(CR=CR, N=N, pm=pm, NM=NM, S=S))
    }  
lev <- levels(grouping)
memberdata <- data.frame(member, grouping = grouping)
memberlist <- split(memberdata, max.col(member))
res <- list()
res$model <- lapply(memberlist, betaregion, lev=lev, dis=dis, eps=eps)
res$eps <- eps
res$member <- betascale(res, member)
return(res)
}



betascale <- function(betaobj, member)
{
if (missing(member)) member <- betaobj$member
else 
    {
    eps <- betaobj$eps
    betaobj <- betaobj$model
    lev <- names(betaobj)
    memberalt <- member
    for (i in seq(along = lev))
        {
        ti <- which(max.col(memberalt) == lev[i])
        memvec <- member[ti, i]
        bets <- betaobj[[i]]
        if (bets$S < eps) {
            memvecneu <- memvec - bets$pm + bets$CR
            memvecneu[memvecneu > 1] <- 1
        }
        else{
            memvecneu <- with(bets, qbeta(pbeta(memvec, pm * NM, (1 - pm) * NM), 
                            CR * N, (1 - CR) * N))
        }
        member[ti,-i] <- ifelse(rep(memvec, length(lev)-1) > (1-1e-6), 
                                (1-memvecneu)/(length(lev)-1), 
                                member[ti,-i]*(1-memvecneu)/(1-memvec))
        member[ti,i] <- memvecneu
        }
    }
return(member)
}
