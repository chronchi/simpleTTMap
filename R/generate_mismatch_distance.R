generate_mismatch_distance <-
function(ttmap_part1_hda, select, alpha = 1){
    D2 <- ttmap_part1_hda$Dc.Dmat[select, ]
    Da <- sign(D2)
    dd <- matrix(rep(0, (dim(Da)[2])^2), nrow = dim(Da)[2])
    ind <- combn(ncol(ttmap_part1_hda$Dc.Dmat[select, ]), 2)
    dd2 <- lapply(as.list(seq_len(ncol(ind))),function(x){
        W <- Da[, ind[1, x]]-Da[, ind[2, x]]
        names(W) <- rownames(Da)
        select1 <- names(W[W[]%in%c(-1, 1)])
        W2 <- abs(D2[, ind[1, x]] - D2[, ind[2, x]])
        names(W2) <- rownames(D2)
        select2 <- names(W2[W2 > alpha])
        s <- length(intersect(select1, select2))
        select1 <- names(W[W[]%in%c(-2, 2)])
        W2 <- abs(D2[,ind[1, x]])
        names(W2) <- rownames(D2)
        select2 <- names(W2[W2 > alpha])
        W2 <- abs(D2[, ind[2, x]])
        names(W2) <- rownames(D2)
        select3 <- names(W2[W2 > alpha])
        s <- s + length(intersect(select1, union(select2, select3)))
        s / (dim(ttmap_part1_hda$Dc.Dmat[select, ])[1])
    })
    dd[lower.tri(dd)] <-  unlist(dd2)
    dd[upper.tri(dd)] <- t(dd)[upper.tri(dd)]
    rownames(dd) <- colnames(D2)
    colnames(dd) <- colnames(D2)
    return(dd)
}
