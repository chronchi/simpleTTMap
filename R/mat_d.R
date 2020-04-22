mat_d <-
function(mat){
    for(i in seq_len(dim(mat)[1])){
        s <- grep(1,mat[i, ])
        for(j in s){
            if(j == i) mat <- mat
            else{
                w <- grep(1, mat[j,])
                mat[i, w] <- rep(1, length(w))
                mat[j, ] <- rep(0, dim(mat)[1])
            }
        }
    }
    return(mat)}
