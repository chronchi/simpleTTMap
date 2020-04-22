ttmap_sgn_genes_inter2 <-
function(q, ttmap_part1_hda, alpha = 0){
    A <- ttmap_part1_hda$Dc.Dmat[,
    colnames(ttmap_part1_hda$Dc.Dmat) %in% q]
    if(length(q) > 1){
        A <- as.matrix(A)
        w <- apply(abs(A), 1, max)
        w2 <- apply(A, 1, max)
        w3 <- apply(A, 1, min)
        w2 <- sign(w2*w3)
    }
    else{
        A <- as.matrix(A)
        rownames(A) <- rownames(ttmap_part1_hda$Dc.Dmat)
        w <- abs(A[A != 0, ])
        w2 <- sign(w)
    }
    select2 <- names(w2[w2[] > 0])
    #we just don't want those which are different
    select <- names(w[w[] > alpha])
    #we will only consider those which are different than alpha
    s <- intersect(select, select2)
    if(length(s) == 1){
        A <- as.matrix(t(A[s, ]))
        rownames(A) <- s
    }
    else{A <- A[s, ]}
    return(A)
}
