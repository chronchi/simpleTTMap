f1prime3 <-
function(xprime, mprime, eprime=1){
    xprime[abs(xprime-apply(
    cbind(xprime, mprime), 1, function(x){
    mean(x[x > 0])})) > 2 * eprime] <- NA
    return(xprime);
}
