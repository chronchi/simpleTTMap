f1prime2 <-
function(xprime, mprime, eprime = 1) {
    xprime[abs(xprime - mean(mprime)) > 
    (2 * sd(mprime) / sqrt(ncol(mprime)))] <- NA
    return(xprime);
}
