mat2pcl <-
function(mat, tag) {
    mat.df <- as.data.frame(rbind(1, mat));
    rownames(mat.df)[[1]] <- "EWEIGHT";
    mat.pcl <- cbind(tag, mat.df);
    return(mat.pcl)}
