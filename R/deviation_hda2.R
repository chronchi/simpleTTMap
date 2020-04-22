deviation_hda2 <-
function(Dmat, Nc.Dmat, new.cnames="Dis"){
    alpha <- function(Da, Nc) {
    Da2 <- Da-Nc
    return(Da2)}
    mat <- alpha(as.matrix(Dmat[
    rownames(Nc.Dmat), ]), as.matrix(Nc.Dmat));
    mat <- ifelse(as.matrix(Dmat[rownames(Nc.Dmat), ]) == 0, 0, mat)
    colnames(mat) <- new.cnames
    return(mat)
}
