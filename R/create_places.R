create_places <-
function(size){
    fa <- rep(0, length(size))
    if(length(size) > 1){
        for(i in 2:length(size)){
            fa[i] <- fa[i - 1] + size[[i - 1]] + size[[i]] + 6}}
    return(fa)
}
