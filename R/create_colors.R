create_colors <-
function(q, muc, size){
    min <- min(muc)
    max <- max(muc)
    k <- (max - min) / 50
    f <- unlist(lapply(seq_len(length(size)), function(i){
        (sum(muc[q[[i]]])) / (length(q[[i]]))}))
    w <- round((f - min) / k)
    names(f) <- names(size)
    w[w == 0] <- 1
    g<- lapply(seq_len(length(size)), function(i){
        (col2rgb(matlab.like(50)[w[i]])) / 255
    })
    junk <- list(col = g, average = f)
    return(junk)
}
