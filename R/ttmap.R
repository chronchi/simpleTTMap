ttmap <- function(ttmap_part1_hda, 
                  m1,
                  select = row.names(ttmap_part1_hda$Dc.Dmat), 
                  ddd, 
                  e,
                  filename = "TEST", n = 3, ad = 0, bd = 0, piq = 1,
                  dd = generate_mismatch_distance(ttmap_part1_hda = ttmap_part1_hda,
                  select = select),
                  mean_value_m1 = "N", 
                  ni = 2,
                  output_directory = getwd(),
                  plot_graph = 0,
                  verbose = FALSE) {    

    write.table(dd, 
                file = file.path(output_directory, "distance_matrix.txt"),
                sep = "\t", quote = FALSE, 
                row.names = TRUE, col.names = NA)

    minespace <- 50
    beg <- ls()
    annot <- function(q, n = n) {
        q <- lapply(seq_len(length(q)), function(i){
        ddd[q[[i]], n]
        })
        return(q)
    }
    
    qmin <- min(m1)
    q25 <- quantile(m1, 0.25)
    q5 <- quantile(m1, 0.5)
    q75 <- quantile(m1, 0.75)
    qmax <- max(m1)
    
    low <- m1[m1[] < q25]
    mid1 <- m1[m1[] < q5 & m1[] >= q25]
    mid2 <- m1[m1[] >= q5 & m1[] < q75]
    high <- m1[m1[] >= q75]
    f <- e
    
    if(length(low) == 0){
        low_map <- as.matrix(0)
    }
    else {
        low_map <- mapper1(as.matrix(
        dd[names(low),names(low)]), e = f)
    }
    
    if(length(mid1) == 0){
        mid1_map <- as.matrix(0)
    }
    else{
        mid1_map <- mapper1(as.matrix(dd[names(mid1),
        names(mid1)]), e = f)
    }
    
    if(length(mid2) == 0){
        mid2_map <- as.matrix(0)
    }
    else{
        mid2_map <- mapper1(as.matrix(dd[names(mid2),
        names(mid2)]), e = f)
    }
    
    if(length(high) == 0){
        high_map <- as.matrix(0)
    }
    else{
        high_map <- mapper1(as.matrix(dd[names(high),
        names(high)]), e = f)
    }
    
    if(length(all) == 0){
        all <- as.matrix(0)
    }
    else{
        all <- mapper1(as.matrix(dd), e = f)
    }
    
    if(bd != 0){
        low_map <- cutoff_low(low_map, 
        para = piq, text = "low")
        mid1_map <- cutoff_low(mid1_map, 
        para = piq, text = "mid1")
        mid2_map <- cutoff_low(mid2_map, 
        para = piq, text = "mid2")
        high_map <- cutoff_low(high_map, 
        para = piq, text = "high")
        all <- cutoff_low(all, para = piq, text = "all")
    }   

    if(dim(all)[1] == 1){
        if(all != 0){
            q_all <- names(all)
            q_all <- as.list(q_all)
            names(q_all) <- names(all)
        }    
        else{
            q_all <- c()
        }
    }
    else{
        q_all <- apply(all, 1, grep, pattern = 1)
        q_all <- as.list(q_all)
        q_all <- lapply(seq_len(length(q_all)), function(i){
        colnames(as.matrix(dd))[q_all[[i]]]
        })
    }

    q1_all <- q_all[lapply(q_all, length) > 0]
    size_all <- lapply(q1_all, length)   

    if(dim(mid1_map)[1] == 1){
        if(mid1_map != 0){
            q_mid1 <- names(mid1)
            q_mid1 <- as.list(q_mid1)
            names(q_mid1) <- names(mid1)
        } 
        else{
            q_mid1 <- c()}
        }
    else{
        q_mid1 <- apply(mid1_map, 1, grep, pattern = 1)
        q_mid1 <- as.list(q_mid1)
        q_mid1 <- lapply(seq_len(length(q_mid1)),  function(i){
        colnames(as.matrix(dd)[names(mid1), names(mid1)])[
        as.vector(q_mid1[[i]])]
        })
    }
    
    q1_mid1 <- q_mid1[lapply(q_mid1, length) > 0]
    size_mid1 <- lapply(q1_mid1,length)
    if(dim(mid2_map)[1] == 1){
        if(mid2_map != 0){
            q_mid2 <- names(mid2)
            q_mid2 <- as.list(q_mid2)
            names(q_mid2) <- names(mid2)
        }
        else{
            q_mid2 <- c()
        }
    }
    else{
        q_mid2 <- apply(mid2_map, 1, grep, pattern = 1)
        q_mid2 <- as.list(q_mid2)
        q_mid2 <- lapply(seq_len(length(q_mid2)), function(i){
        colnames(as.matrix(dd)[names(mid2), names(mid2)])[
        as.vector(q_mid2[[i]])]
        })
    }
    
    q1_mid2 <- q_mid2[lapply(q_mid2, length) > 0]
    size_mid2 <- lapply(q1_mid2, length)
    if(dim(low_map)[1] == 1){
        if(low_map != 0){
            q_low <- names(low)
            q_low <- as.list(q_low)
            names(q_low) <- names(low)
        }
        else{
            q_low <- c()
        }
    }
    else{
        q_low <- apply(low_map, 1, grep, pattern = 1)
        q_low <- as.list(q_low)
        q_low <- lapply(seq_len(length(q_low)),  function(i){
        colnames(as.matrix(dd)[names(low), names(low)])[
        as.vector(q_low[[i]])]
        })
    }
    
    q1_low <- q_low[lapply(q_low, length) > 0]
    size_low <- lapply(q1_low, length)
    if(dim(high_map)[1] == 1){
        if(high_map != 0){
            q_high <- names(high)
            q_high <- as.list(q_high)
            names(q_high) <- names(high)
        }
        else{
            q_high <- c()
        }
    }
    else{
        q_high <- apply(high_map, 1, grep, pattern = 1)
        q_high <- as.list(q_high)
        q_high <- lapply(seq_len(length(q_high)), function(i){
        colnames(as.matrix(dd)[names(high), names(high)])[
        as.vector(q_high[[i]])]
        })
    }

    q1_high <- q_high[lapply(q_high, length) > 0]
    size_high <- lapply(q1_high, length)
    
   squize <- function(q1_all, m1, size_all, n = n){
        r_all <- create_colors(q1_all, m1, size_all)
        sort_r <- sort(r_all$average, index.return = TRUE)
        size_all <- size_all[sort_r$ix]
        q1_all <- q1_all[sort_r$ix]
        q1_all_a <- annot(q1_all, n = n)
        f_all <- create_places(size_all)
        out1 <- list(f = f_all, q1_all = q1_all)
        return(out1)
    }
    p <- squize(q1_all, m1, size_all, n = n)
    #lans <- lapply(seq_len(length(size_all)), function(i){
    #    spheres3d(p$f[i], 0, 0, radius = p$s[i],
    #    color=rgb((p$r)$col[[i]][1], (p$r)$col[[i]][2],
    #    (p$r)$col[[i]][3], alpha = 1))
    #    text3d(p$f[i], max(as.matrix(unlist(p$s))) + 5, 0, i)
    #    if(ad == 0){text3d(p$f[i], 
    #        -max(as.matrix(unlist(p$s))) - 5, 0, 
    #        paste(unique(p$q1_all_a[[i]]), collapse=""))
    #    }
    #})
    m <- 0
    if(dim(low_map)[1] == 1 && low_map == 0){
        du <- max(as.matrix(unlist(size_all))) + minespace
        d_low <- du
        du <- du + minespace
        size_low <- c()
        p_low <- list()
        p_low$f <- 0
    }
    else{
        p_low <- squize(q1_low, m1, size_low, n=n)
        du <- max(as.matrix(unlist(size_all))) + 
        max(as.matrix(unlist(size_low))) + minespace
        d_low <- du
        #lans <- lapply(seq_len(length(p_low$s)), function(i){
        #spheres3d(p_low$f[i], du, 
        #0, radius = p_low$s[i], 
        #color = rgb((p_low$r)$col[[i]][1],
        #(p_low$r)$col[[i]][2], (p_low$r)$col[[i]][3],alpha = 1))
        #text3d(p_low$f[i], 
        #du + max(as.matrix(unlist(size_low))) + 5, 0,
        #length(size_all) + i)
        #if(ad == 0){text3d(p_low$f[i], 
        #    (du - max(as.matrix(unlist(size_low))) - 5), 0,
        #    paste(unique(p_low$q1_all_a[[i]]), collapse = ""))
        #}
        #})
        #create_links(q = p$q1_all, q1 = p_low$q1_all, f = p$f,
        #f1 = p_low$f, m, du, 0)
        du <- du + max(as.matrix(unlist(size_low))) + minespace
    }
    if(dim(mid1_map)[1] == 1 && mid1_map == 0){
        d_mid1 <- du
        du <- du + minespace
        size_mid1 <- c()
        p_mid1 <- list()
        p_mid1$f <- 0}
    else{
        m <- max(p_low$f)
        p_mid1 <- squize(q1_mid1, m1, size_mid1, n = n)
        d_mid1 <- du + max(as.matrix(unlist(size_mid1)))
        du <- d_mid1
        #l <- 10
        #lans <- lapply(seq_len(length(size_mid1)), function(i){
        #    spheres3d(m + p_mid1$f[i], du, l,
        #    radius = p_mid1$s[i], color = rgb((p_mid1$r)$col[[i]][1],
        #    (p_mid1$r)$col[[i]][2], (p_mid1$r)$col[[i]][3],alpha = 1))
        #    text3d(m + p_mid1$f[i], du + max(as.matrix(
        #    unlist(size_mid1))) + 5, l,
        #    length(size_all) + length(size_low) + i)
        #    if(ad == 0){text3d(m + p_mid1$f[i],
        #        (du - max(as.matrix(unlist(size_mid1))) -5), l,
        #        paste(unique(p_mid1$q1_all_a[[i]]), collapse = ""))}
        #})
        #create_links(p$q1_all, p_mid1$q1_all, p$f, p_mid1$f, m, du, l)
        du <- du + max(as.matrix(unlist(size_mid1)))+minespace
    }
    m <- max(p_mid1$f) + m
    if(dim(mid2_map)[1] == 1 && mid2_map == 0){
        d_mid2 <- du
        du <- du + minespace
        p_mid2 <- list()
        p_mid2$f <- 0
        size_mid2 <- c()}
    else{
        du <- du + max(as.matrix(unlist(size_mid2)))
        d_mid2 <- du
        p_mid2 <- squize(q1_mid2, m1, size_mid2, n = n)
        #l <- 20
        #lans <- lapply(seq_len(length(size_mid2)),function(i){
        #    spheres3d(m + p_mid2$f[i],
        #    du, l, radius = p_mid2$s[i], 
        #    color = rgb((p_mid2$r)$col[[i]][1],
        #    (p_mid2$r)$col[[i]][2], (p_mid2$r)$col[[i]][3],alpha = 1))
        #    text3d(m + p_mid2$f[i],
        #    du + max(as.matrix(unlist(size_mid2))) + 5,
        #    l,length(size_all) + length(size_low) +
        #    length(size_mid1) + i)
        #    if(ad == 0){text3d(m + p_mid2$f[i],
        #        (du - max(as.matrix(unlist(size_mid2))) -5), l,
        #        paste(unique(p_mid2$q1_all_a[[i]]), 
        #        collapse = ""))}
        #})
        #create_links(p$q1_all, p_mid2$q1_all, p$f, p_mid2$f, m, du, l)
        du <- du + minespace +  max(as.matrix(unlist(size_mid2)))
    }
    m <- max(p_mid2$f) + m   
    if(dim(high_map)[1] == 1 && high_map == 0){
        d_high <- du
        p_high <- list()
        p_high$f <- 0
        size_high <- c()}
    else{
        p_high <- squize(q1_high, m1, size_high, n=n)
        du <- du + max(as.matrix(unlist(size_high)))
        d_high <- du
        #l <- 30
        #lans <- lapply(seq_len(length(size_high)), function(i){
        #    spheres3d(m + p_high$f[i], du, l,
        #    radius = p_high$s[i], color = rgb((p_high$r)$col[[i]][1],
        #    (p_high$r)$col[[i]][2], (p_high$r)$col[[i]][3], alpha=1))
        #    text3d(m + p_high$f[i],
        #    du + max(as.matrix(unlist(size_high))) + 5, l,
        #    length(size_all) + 
        #    length(size_low) + length(size_mid1) +
        #    length(size_mid2) + i)
        #    if(ad == 0){text3d(m + p_high$f[i],
        #        (du - max(as.matrix(unlist(size_high))) -5), l,
        #        paste(unique(p_high$q1_all_a[[i]]), 
        #        collapse = ""))}
        #})
        #create_links(p$q1_all, p_high$q1_all, p$f, p_high$f, m, du, l)
    }
    #annot_right(p_high, p, l, d_high, d_mid2, d_mid1, d_low, m)
    #u <- min(- minespace,-max(as.matrix(unlist(size_all))) - 10)
    #lans <- lapply(seq_len(100), function(i){
    #    segments3d(c((((i - 1) / 100)*(max(m + max(p_high$f),
    #    max(p$f)))), ((i / 100)*(max(m + max(p_high$f),
    #    max(p$f))))), c(u, u), c(l, l),
    #    col= matlab.like2(100)[i], lwd = 20)
    #})
    #text3d(-1-100, u, l, "Mean Deviation")
    #if(mean_value_m1 == "N"){
    #    text3d(-40, u, l, round(qmin / length(select),2))
    #    text3d(max(m + max(p_high$f),
    #    max(p$f)) + 40, u, l, 
    #    round(qmax / length(select),2))
    #}
    #else{
    #    text3d(-40, u, l, round(qmin, 2))
    #    text3d(max(m + max(p_high$f), max(p$f)) + 40, 
    #    u, l, round(qmax, 2))
    #}
    
    # descriptions files
    e <- list()
    if(length(p$q1_all) != 0){for(i in seq_len(length(size_all))){
            e_1 <- paste(p$q1_all[[i]][1],
            paste(ddd[p$q1_all[[i]][1], ni], ")"), sep=" (")
            if(length(p$q1_all[[i]]) > 1){
                for(j in 2:length(p$q1_all[[i]])){
                    e_1 <- paste(paste(e_1, 
                    p$q1_all[[i]][j], sep=","),
                    paste(ddd[p$q1_all[[i]][j], ni], ")"), sep=" (")
                }
            }
            e[[i]] <- paste(i, e_1, sep = ":" )
    }}
    
    d <- length(size_all)
    e[[d + 1]]<- "The first quartile"
    if(length(p_low$q1_all) != 0){for(i in seq_len(length(size_low))){
            e_1 <- paste(p_low$q1_all[[i]][1],
            paste(ddd[p_low$q1_all[[i]][1], ni], ")"),sep="(")
            if(length(p_low$q1_all[[i]]) > 1){
                for(j in 2:length(p_low$q1_all[[i]])){
                    e_1 <- paste(paste(e_1,
                    p_low$q1_all[[i]][j], sep=","),
                    paste(ddd[p_low$q1_all[[i]][j], ni], ")"), sep="(")
                }
            }
            e[[i + d + 1]] <- paste(i + d - 1, e_1, sep = ":" )
        }
    }   
    
    d <- length(size_all) + length(size_low) + 1
    e[[d+1]]<- "The second quartile"
    if(length(p_mid1$q1_all) != 0){
        for(i in seq_len(length(size_mid1))){
            e_1 <- paste(p_mid1$q1_all[[i]][1],
            paste(ddd[p_mid1$q1_all[[i]][1], ni], ")"),sep="(")
            if(length(p_mid1$q1_all[[i]]) > 1){
                for(j in 2:length(p_mid1$q1_all[[i]])){
                e_1 <-paste(paste(e_1,p_mid1$q1_all[[i]][j], 
                sep = ","),
                paste(ddd[p_mid1$q1_all[[i]][j], ni],")"),
                sep="(")
                }
            }
            e[[i + d + 1]]<- paste(i + d - 2, e_1, sep = ":" )
        }
    }  
    
    d <- length(size_all) + length(size_low) + length(size_mid1) + 2
    e[[d + 1]]<- "The third quartile"
    if(length(p_mid2$q1_all) != 0){
        for(i in seq_len(length(size_mid2))){
            e_1 <- paste(p_mid2$q1_all[[i]][1],
            paste(ddd[p_mid2$q1_all[[i]][1], ni], ")"),
            sep = "(")
            if(length(p_mid2$q1_all[[i]]) > 1){
                for(j in 2:length(p_mid2$q1_all[[i]])){
                e_1 <- paste(paste(e_1, p_mid2$q1_all[[i]][j], 
                sep = ","),
                paste(ddd[p_mid2$q1_all[[i]][j], ni], ")"), 
                sep = "(")
                }
            }
            e[[i + d + 1]]<-paste(i + d - 3, e_1, sep = ":" )
        }
    }
    d <- length(size_all) + length(size_low) +
        length(size_mid1) + length(size_mid2) + 3
    
    e[[d + 1]]<- "The fourth quartile"
    if(length(p_high$q1_all) != 0){
        for(i in seq_len(length(size_high))){
            e_1 <- paste(p_high$q1_all[[i]][1],
            paste(ddd[p_high$q1_all[[i]][1], ni], ")"), 
            sep = "(")
            if(length(p_high$q1_all[[i]]) > 1){
                for(j in 2:length(p_high$q1_all[[i]])){
                e_1 <- paste(paste(e_1, p_high$q1_all[[i]][j], 
                sep = ","),
                paste(ddd[p_high$q1_all[[i]][j], ni], ")"), 
                sep = "(")
                }
            }
            e[[i + d + 1]]<-paste(i + d - 4, e_1, sep = ":" )
        }
    }
    
    write.table(e, file=file.path(output_directory, paste(filename,
        "description.txt", sep = "_")),
        quote = FALSE, sep = "\n", row.names = FALSE,
        col.names = FALSE)
    
    end_out <- list(low_map = p_low$q1_all,
                    mid1_map=p_mid1$q1_all,
                    mid2_map=p_mid2$q1_all, 
                    high_map = p_high$q1_all,
                    all = p$q1_all)
    
    vari <- setdiff(ls(), list(beg)[[1]])
    vari <- setdiff(list(vari)[[1]], "end_out")
    rm(list = vari)
    
    return(end_out)
}
