calcul_e_single<-
function(dd5, pvalcutoff=0.95, tt1, 
alpha = 1, S=colnames(tt1$Normal.mat)){
    er <- tt1$Normal.mat[,
    colnames(tt1$Normal.mat)%in%S]
    er[er == 0] <- NA
    w <- apply(er, 1, mean, na.rm = TRUE)
    w2 <- apply(er, 1, sd, na.rm = TRUE)
    w <- ifelse(is.na(w) == TRUE, 0, w)
    w2 <- ifelse(is.na(w2) == TRUE, 0, w2)
    m_f <- apply(tt1$flat.Nmat$mat[
    rownames(tt1$Normal.mat),
    colnames(tt1$flat.Nmat$mat)%in%S], 1, min)
    m2_f <- apply(tt1$flat.Nmat$mat[
    rownames(tt1$Normal.mat),
    colnames(tt1$flat.Nmat$mat)%in%S], 1, max)
    p2 <- lapply(seq_len(length(m_f)), function(i){
        pnorm(m_f[i], w[i], w2[i], lower.tail = TRUE)
    })
    p2 <- unlist(p2)
    p2p <- lapply(seq_len(length(m_f)), function(i){
        pnorm(m_f[i] - alpha, w[i], w2[i], lower.tail = TRUE)})
    p2p <- unlist(p2p)
    p3 <- lapply(seq_len(length(m_f)), function(i){ 
        pnorm(m2_f[i], w[i], w2[i], lower.tail = TRUE)})
    p4 <- lapply(seq_len(length(m_f)), function(i){ 
        pnorm(m_f[i], w[i], w2[i], lower.tail = TRUE)})
    p3 <- unlist(p3)
    p4 <- unlist(p4)
    p3 <- p3 - p4
    p1 <- lapply(seq_len(length(m_f)), function(i){ 
        pnorm(m2_f[i], w[i], w2[i], lower.tail = FALSE)})
    p1p <- lapply(seq_len(length(m_f)), function(i){
        pnorm(m2_f[i] + alpha, w[i], w2[i], lower.tail = FALSE)})
    p1 <- unlist(p1)
    p1p <- unlist(p1p)
    a3 <- 2 * ((p1 * p2p) + (p2 * p1p) - (p1p * p2p))
    sl <- lapply(seq_len(length(as.vector(
    dd5[upper.tri(dd5)]))), function(i){ 
        ele <- as.vector(dd5[upper.tri(dd5)])[i]
        #s <- pnorm(ele*length(m_f), mean = sum(a3,na.rm=TRUE),
        # sd = sqrt(sum(a3*(1-a3),na.rm=TRUE)), lower.tail = TRUE)
        s <- ppois(ele * length(m_f), sum(a3), 
        lower.tail = TRUE, log.p = FALSE)
        ifelse(s < pvalcutoff, ele, 0)
    })
    s <- max(unlist(sl))
    if(s <= 0){
        sl <- lapply(seq_len(length(
        as.vector(dd5[upper.tri(dd5)]))),function(i){ 
            ele <- as.vector(dd5[upper.tri(dd5)])[i]
            s <- ppois(ele * length(m_f), sum(a3),
            lower.tail = TRUE, log.p = FALSE)
            ifelse(s < 1-.Machine$double.neg.eps, ele, 0)
        })
        s <- max(unlist(sl))}
    return(s)
}
