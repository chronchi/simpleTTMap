setGeneric("make_matrices",function(mat, col_ctrl, col_test,
NAME, CLID, GWEIGHT = rep(1,
dim(mat)[1]), EWEIGHT = 0)standardGeneric("make_matrices"))
setMethod("make_matrices", signature(mat = "matrix"),
    function(mat, col_ctrl, col_test,
    NAME, CLID, GWEIGHT = rep(1, dim(mat)[1]), EWEIGHT = 0){
        experiment <- list()
        experiment$CTRL <- as.data.frame(mat[, col_ctrl])
        experiment$CTRL <- cbind(CLID, NAME, GWEIGHT, experiment$CTRL)
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        if(length(EWEIGHT) == 1){
            EWEIGHT<- rep(1, dim(mat[, col_ctrl])[2] + 3)
        }
        experiment$CTRL <- rbind(EWEIGHT, experiment$CTRL)
        rownames(experiment$CTRL)[1] <- "EWEIGHT"
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        experiment$CTRL[, 1] <- as.character(experiment$CTRL[,1])
        experiment$CTRL[1, 1] <- "EWEIGHT"
        colnames(experiment$CTRL)[c(1, 2)] <- c("CLID", "NAME")
        experiment$TEST <- as.data.frame(mat[, col_test])
        experiment$TEST <- cbind(CLID, NAME, GWEIGHT, experiment$TEST)
        experiment$TEST <- as.data.frame(experiment$TEST)
        if(length(EWEIGHT) == 1){
            EWEIGHT <- rep(1, dim(mat[, col_test])[2] + 3)
        }
        experiment$TEST <- rbind(EWEIGHT, experiment$TEST)
        rownames(experiment$TEST)[1] <- "EWEIGHT"
        experiment$TEST <- as.data.frame(experiment$TEST)
        experiment$TEST[, 1] <- as.character(experiment$TEST[,1])
        experiment$TEST[1, 1] <- "EWEIGHT"
        colnames(experiment$TEST)[c(1, 2)] <- c("CLID", "NAME")
        return(experiment)
    }
)
setMethod("make_matrices", signature(mat = "data.frame"), 
    function(mat, col_ctrl, col_test,
    NAME, CLID, GWEIGHT = rep(1,dim(mat)[1]), EWEIGHT = 0){
        experiment <- list()
        experiment$CTRL <- mat[, col_ctrl]
        experiment$CTRL <- cbind(CLID, NAME, GWEIGHT, experiment$CTRL)
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        if(length(EWEIGHT) == 1){
            EWEIGHT<- rep(1, dim(mat[,col_ctrl])[2] + 3)
        }    
        experiment$CTRL <- rbind(EWEIGHT, experiment$CTRL)
        rownames(experiment$CTRL)[1] <- "EWEIGHT"
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        experiment$CTRL[, 1] <- as.character(experiment$CTRL[, 1])
        experiment$CTRL[1, 1] <- "EWEIGHT"
        colnames(experiment$CTRL)[c(1, 2)] <- c("CLID", "NAME")
        experiment$TEST <- as.data.frame(mat[, col_test])
        experiment$TEST <- cbind(CLID, NAME, GWEIGHT, experiment$TEST)
        experiment$TEST <- as.data.frame(experiment$TEST)
        if(length(EWEIGHT) == 1){
            EWEIGHT <- rep(1, dim(mat[, col_test])[2] + 3)
        }
        experiment$TEST <- rbind(EWEIGHT, experiment$TEST)
        rownames(experiment$TEST)[1] <- "EWEIGHT"
        experiment$TEST <- as.data.frame(experiment$TEST)
        experiment$TEST[, 1] <- as.character(experiment$TEST[, 1])
        experiment$TEST[1, 1] <- "EWEIGHT"
        colnames(experiment$TEST)[c(1, 2)] <- c("CLID", "NAME")
        return(experiment)
    }
)
setMethod("make_matrices", signature(mat = "ExpressionSet"),
    function(mat, col_ctrl, col_test,
    NAME, CLID, GWEIGHT = rep(1, dim(mat)[1]), EWEIGHT = 0){mat <- 
            as(mat, "RangedSummarizedExperiment")
        f <- (SummarizedExperiment::assay(mat))
        f <- lapply(seq_len(dim(f)[1]), function(x){ w <- unlist(f[x, ]); 
            w[is.na(w) == TRUE] <- median(w, na.rm = TRUE);return(w)})
        f <- matrix(unlist(f), ncol = length(colnames(mat)))
        rownames(f) <- rownames(mat)
        colnames(f) <- colnames(mat)
        mat <- f
        experiment <- list()
        experiment$CTRL <- as.data.frame(mat[, col_ctrl])
        experiment$CTRL <- cbind(CLID, NAME, GWEIGHT, experiment$CTRL)
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        if(length(EWEIGHT) == 1){
            EWEIGHT<- rep(1, dim(mat[,col_ctrl])[2] + 3)
        }
        experiment$CTRL <- rbind(EWEIGHT, experiment$CTRL)
        rownames(experiment$CTRL)[1] <- "EWEIGHT"
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        experiment$CTRL[, 1] <- as.character(experiment$CTRL[,1])
        experiment$CTRL[1, 1] <- "EWEIGHT"
        colnames(experiment$CTRL)[c(1, 2)] <- c("CLID", "NAME")
        experiment$TEST <- as.data.frame(mat[, col_test])
        experiment$TEST <- cbind(CLID, NAME, GWEIGHT, experiment$TEST)
        experiment$TEST <- as.data.frame(experiment$TEST)
        if(length(EWEIGHT) == 1){
            EWEIGHT <- rep(1, dim(mat[, col_test])[2] + 3)
        }
        experiment$TEST <- rbind(EWEIGHT, experiment$TEST)
        rownames(experiment$TEST)[1] <- "EWEIGHT"
        experiment$TEST <- as.data.frame(experiment$TEST)
        experiment$TEST[, 1] <- as.character(experiment$TEST[,1])
        experiment$TEST[1, 1] <- "EWEIGHT"
        colnames(experiment$TEST)[c(1, 2)] <- c("CLID", "NAME")
        return(experiment)
    }
)
setMethod("make_matrices", signature(mat = "SummarizedExperiment"),
    function(mat, col_ctrl, col_test,
    NAME, CLID, GWEIGHT = rep(1,dim(mat)[1]), EWEIGHT = 0){
            mat <- as(mat, "RangedSummarizedExperiment")
        f <- (SummarizedExperiment::assay(mat))
        f <- lapply(seq_len(dim(f)[1]), function(x){ w <- unlist(f[x, ]); 
            w[is.na(w) == TRUE] <- median(w, na.rm = TRUE);return(w)})
        f <- matrix(unlist(f), ncol = length(colnames(mat)))
        rownames(f) <- rownames(mat)
        colnames(f) <- colnames(mat)
        mat <- f
        experiment <- list()
        experiment$CTRL <- as.data.frame(mat[, col_ctrl])
        experiment$CTRL <- cbind(CLID,NAME,GWEIGHT,experiment$CTRL)
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        if(length(EWEIGHT) == 1){
            EWEIGHT<- rep(1, dim(mat[,col_ctrl])[2] + 3)
        }
        experiment$CTRL <- rbind(EWEIGHT, experiment$CTRL)
        rownames(experiment$CTRL)[1] <- "EWEIGHT"
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        experiment$CTRL[, 1] <- as.character(experiment$CTRL[,1])
        experiment$CTRL[1, 1] <- "EWEIGHT"
        colnames(experiment$CTRL)[c(1, 2)] <- c("CLID","NAME")
        experiment$TEST <- as.data.frame(mat[, col_test])
        experiment$TEST <- cbind(CLID, NAME, GWEIGHT, experiment$TEST)
        experiment$TEST <- as.data.frame(experiment$TEST)
        if(length(EWEIGHT) == 1){
            EWEIGHT <- rep(1, dim(mat[, col_test])[2] + 3)
        }
        experiment$TEST <- rbind(EWEIGHT, experiment$TEST)
        rownames(experiment$TEST)[1] <- "EWEIGHT"
        experiment$TEST <- as.data.frame(experiment$TEST)
        experiment$TEST[, 1] <- as.character(experiment$TEST[,1])
        experiment$TEST[1, 1] <- "EWEIGHT"
        colnames(experiment$TEST)[c(1, 2)] <- c("CLID", "NAME")
        return(experiment)
    }
)
setMethod("make_matrices", signature(mat = "RangedSummarizedExperiment"),
    function(mat, col_ctrl, col_test,
    NAME, CLID, GWEIGHT = rep(1, dim(mat)[1]),  EWEIGHT = 0){
        f <- (SummarizedExperiment::assay(mat))
        f <- lapply(seq_len(dim(f)[1]), function(x){ w <- unlist(f[x, ]); 
            w[is.na(w) == TRUE] <- median(w, na.rm = TRUE);return(w)})
        f <- matrix(unlist(f), ncol = length(colnames(mat)))
        rownames(f) <- rownames(mat)
        colnames(f) <- colnames(mat)
        mat <- f
        experiment <- list()
        experiment$CTRL <- as.data.frame(mat[, col_ctrl])
        experiment$CTRL <- cbind(CLID, NAME, GWEIGHT, experiment$CTRL)
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        if(length(EWEIGHT) == 1){
            EWEIGHT<- rep(1, dim(mat[,col_ctrl])[2] + 3)
        }
        experiment$CTRL <- rbind(EWEIGHT, experiment$CTRL)
        rownames(experiment$CTRL)[1] <- "EWEIGHT"
        experiment$CTRL <- as.data.frame(experiment$CTRL)
        experiment$CTRL[, 1] <- as.character(experiment$CTRL[, 1])
        experiment$CTRL[1, 1] <- "EWEIGHT"
        colnames(experiment$CTRL)[c(1, 2)] <- c("CLID","NAME")
        experiment$TEST <- as.data.frame(mat[, col_test])
        experiment$TEST <- cbind(CLID, NAME, GWEIGHT, experiment$TEST)
        experiment$TEST <- as.data.frame(experiment$TEST)
        if(length(EWEIGHT) == 1){
            EWEIGHT <- rep(1, dim(mat[, col_test])[2] + 3)
        }
        experiment$TEST <- rbind(EWEIGHT, experiment$TEST)
        rownames(experiment$TEST)[1] <- "EWEIGHT"
        experiment$TEST <- as.data.frame(experiment$TEST)
        experiment$TEST[, 1] <- as.character(experiment$TEST[, 1])
        experiment$TEST[1, 1] <- "EWEIGHT"
        colnames(experiment$TEST)[c(1, 2)] <- c("CLID","NAME")
        return(experiment)
    }
)