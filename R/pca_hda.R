pca_hda <- function (mat, j) 
{
    td <- function(vec) {
        ifelse(length(vec) > 1.5, y <- diag(vec),
            y <- as.matrix(vec))
        return(y)}
    pn <- function(vec) {
        vn <- vec
        vn[-(seq_len(j))] <- rep(0, length(vec) - j)
        return(vn)}
    mat.svd <- svd(mat)
    U.mat <- mat.svd$u
    V.mat <- mat.svd$v
    D.mat <- td(pn(vec = mat.svd$d))
    matt <- U.mat %*% D.mat %*% t(V.mat)
    attributes(matt) <- attributes(mat)
    return(matt)}
