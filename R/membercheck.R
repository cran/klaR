membercheck <- function(member){
    member <- as.matrix(member)
    if(any(member) < 0 || any(member) > 1) 
        stop("Membership values (posterior probabilities) must be in [0,1].")
    if(!identical(all.equal(rowSums(member), rep(1, nrow(member))), TRUE))
        stop("Membership values (posterior probabilities) must sum up to 1.")
}
