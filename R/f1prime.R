# f1prime <-
# function(xprime, mprime, eprime = 1) {
    # xprime[abs(xprime - apply(cbind(xprime, mprime), 1,
    # function(x){median(x[x > 0])})) > eprime] <- NA
    # return(xprime);
# }
f1prime <-
function(xprime, mprime, eprime = 1) {
    xprime[abs(xprime - apply(mprime, 1,
    function(x){median(x[x > 0])})) > eprime] <- NA
    return(xprime);
}
