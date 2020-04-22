impute_median <-
function(ma, med){
    ma <- t(apply(ma, 1, function(ma) ifelse(is.na(ma), 
    median(ma, na.rm = TRUE), ma)))
    ma[is.na(rowSums(ma)) == TRUE,] <-
    med[is.na(rowSums(ma)) == TRUE]
    return(ma)
}
