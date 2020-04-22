generate_p_val_correlation <-
function(ttmap_part1_hda, select){
    dd <- as.dist((1 - cor(ttmap_part1_hda$Dc.Dmat[select, ])) / 2)
    dd <- as.matrix(dd)
    for(i in seq_len(dim(dd)[1])){
        for(j in seq_len(dim(dd)[2])){
            dd[i, j] <- cor.test(ttmap_part1_hda$Dc.Dmat[select, i],
            ttmap_part1_hda$Dc.Dmat[select, j],
            alternative = "greater")$p.value
        }
    }
    return(dd)
}
