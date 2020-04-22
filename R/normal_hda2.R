normal_hda2 <-
function(Dmat, Nmodel, new.cnames = "Norm") {
    alpha <- function(Da, N) {
        Da2 <- ifelse(N[,1] <= Da, Da, N[,1])
        Da2 <-  ifelse(Da2 <= N[,2], Da2, N[,2])
        if(dim(as.matrix(Nmodel))[2] == 1){
            Da2 <- as.matrix(Da2)
            rownames(Da2) <- rownames(N)
        }
        return(Da2)}
    if(dim(as.matrix(Nmodel))[2] == 1){ 
        N <- cbind(as.data.frame(Nmodel), as.data.frame(Nmodel))
        rownames(N) <- names(Nmodel)
    }
    else{
        options(warn=-1)
        Nmodel[Nmodel == 0] <- NA
        min <- apply(Nmodel, 1, min, na.rm = TRUE)
        max <- apply(Nmodel, 1, max, na.rm = TRUE)
        min[min == Inf] <- 0
        max[max == -Inf] <- 0
        options(warn=0)
        N <- cbind(min,max)}
    mat <- alpha(as.matrix(Dmat), N);
    colnames(mat) <- new.cnames;
    return(mat)
}
