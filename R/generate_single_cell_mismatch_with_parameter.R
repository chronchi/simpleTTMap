generate_single_cell_mismatch_with_parameter <-
function(ttmap_part1_hda, select, alpha = 1){
    D2 <- ttmap_part1_hda$Dc.Dmat[select, ]
    Da <- sign(D2)
    dd <- matrix(rep(0,
    (dim(Da)[2])^2), nrow = dim(Da)[2])
    ind <- combn(ncol(ttmap_part1_hda$Dc.Dmat[select, ]), 2)
    dd2 <- lapply(as.list(seq_len(ncol(ind))), function(x){
        W2 <- abs(D2[, ind[1, x]])
        W3 <-  abs(D2[, ind[2, x]])
        names(W2) <- rownames(D2)
        names(W3) <- rownames(D2)
        select2 <- names(W2[W2 > alpha])
        select3 <- names(W3[W3 > alpha])
        select2 <- union(select2, select3)
        W <- Da[, ind[1, x]] * Da[, ind[2, x]]
        names(W) <- rownames(Da)
        select1 <- names(W[W[] == (-1)])
        length(intersect(select1, unique(select2))) / length(W)
    })
    dd[lower.tri(dd)] <- unlist(dd2)
    dd[upper.tri(dd)] <- t(dd)[upper.tri(dd)]
    rownames(dd) <- colnames(D2)
    colnames(dd) <- colnames(D2)
    return(dd)
}
