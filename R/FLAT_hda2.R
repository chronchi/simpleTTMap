FLAT_hda2 <-
function(mat, el = 1, metho = 0, p = 1.1,
method_mean_or_median = 1) {
    if(dim(mat)[2] == 1){
        d <- list(m = 0, n = 0, mat = mat, na_numbers = 0)
    }
    else{
        bigindx <- rbind(seq_len(ncol(mat)));  
        colnames(bigindx) <- colnames(mat);
        matt <- apply(bigindx, 2, function(j){
            if(method_mean_or_median == 1){
                f1prime(mat[, j], mat[, -(j)], eprime = el)}
            else{f1prime3(mat[, j], mat[, -(j)], eprime = el)}
        })
        rownames(matt) <- rownames(mat)
        m <- colSums(is.na(matt))
        na_numbers <- as.matrix(
        (rowSums(is.na(matt))) / dim(matt)[2])
        n <- dim(matt)[1] - dim(matt[
        rowSums(is.na(matt)) < p * dim(matt)[2], ])[1]
        if(metho == 0){
            mattt <- as.matrix(impute_median(matt[
            rowSums(is.na(matt)) < p * dim(matt)[2], ],
            apply(mat, 1, median)))
        }else{
            mattt <- as.matrix(impute_median(matt[
            rowSums(is.na(matt)) < p * dim(matt)[2],],
            metho))
        }
        d <- list(m = m, n = n, mat = mattt, na_numbers = na_numbers)
    }
    return(d)
}
