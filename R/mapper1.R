mapper1 <- function(mat, e = 0.1){
    # print(paste("e_map =", e))
    m <- mat_1_0(mat, e=e)
    i <- 0
    while(sum(m) > dim(m)[1] & i < dim(m)[1]){
        i <- i + 1
        m <- mat_d(m)
        m <- mat_d_col(m)
    }
    return(m)
}
