hyperrectangle_deviation_assessment <- function(x, 
                                                dataname, 
                                                normalname, 
                                                output_directory = getwd(),
                                                k = dim(x$Normal.mat)[2])
{
    u <- meshRows_hda(x$tag.pcl, x$flat.Nmat$mat)
    tag.pcl <- u[[1]]
    tag.pcl <- as.data.frame(rbind(1, tag.pcl)); 
    rownames(tag.pcl)[[1]] <- "EWEIGHT"
    s <- meshRows_hda(x$Normal.mat, x$flat.Nmat$mat)
    t <- meshRows_hda(s[[1]], x$Disease.mat)
    Normal.mat <- s[[1]];
    Disease.mat <- t[[2]];
    flat.Nmat <- s[[2]][rownames(Disease.mat),];
    rm(s, t, u)
    B <- x$B
    U1 <- x$U1
    U2 <- x$U2
    rm(x);
    
    if(k == dim(Normal.mat)[2]){
        Normal.model <- flat.Nmat
    }
    else{
        Normal.model <- pca_hda(mat = flat.Nmat, j = k); 
        print("K!=dim..")
    }

    Nc.Dmat <- lapply(seq_len(length(unique(B))), function(i){
        normal_hda2(Dmat = Disease.mat[, U2[[i]]],
        Nmodel = Normal.model[, U1[[i]]],
        new.cnames = paste(U2[[i]], "Norm", sep = ".")); 
    })
    
    Nc.Dmat <- do.call(cbind, Nc.Dmat)
    rownames(Nc.Dmat) <- rownames(Normal.mat)

    # calculate the deviation components given the normalized control group
    Dc.Dmat <- deviation_hda2(Dmat = Disease.mat[,unlist(U2)], 
                              Nc.Dmat, new.cnames = paste(unlist(U2), 
                              "Dis", sep = ".")); 
    rownames(Dc.Dmat) <- rownames(Normal.mat)
    
    # save deviation components
    Dc.Dpcl <- mat2pcl(mat = Dc.Dmat, tag = tag.pcl);
    write_pcl(Dc.Dpcl, 
              paste(dataname, "Tdis", sep = "."), 
              output_directory = output_directory);
   
    # save normal components of test samples
    Nc.Dpcl <- mat2pcl(mat = Nc.Dmat, tag = tag.pcl);
    write_pcl(Nc.Dpcl, 
              paste(dataname, "Tnorm", sep = "."),
              output_directory = output_directory);
    
    # save the control group 
    Normal.model <- mat2pcl(mat = Normal.model, tag = tag.pcl);
    write_pcl(Normal.model,
              paste(normalname, "NormalModel", sep = "."),
              output_directory = output_directory); 
    
    # total absolute deviation 
    sumabs <- function(vec){ 
        y <- abs(vec) 
        z <- sum(y)
        return(z) 
    }
    
    # apply the total deviation sum to every sample 
    m <- apply(Dc.Dpcl[, -(seq_len(3))], 2, sumabs)

    end_out <- list(Dc.Dmat = Dc.Dmat, m = m);
    
    return(end_out)
};
