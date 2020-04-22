\name{make_matrices}
\alias{make_matrices}
\title{
    Prepares the matrices for 
    \code{\link[TTMap]{control_adjustment}}
}
\description{
    \code{\link[TTMap]{make_matrices}} generates the control 
    and the test matrice in the 
    right format
}
\usage{
    make_matrices(mat, col_ctrl, col_test, NAME, CLID,
    GWEIGHT = rep(1, dim(mat)[1]), EWEIGHT = 0)
}
\arguments{
    \item{mat}{the gene expressions can be \bold{matrix}, \bold{data.frame},
    \code{"\link[=RangedSummarizedExperiment-class]{
    RangedSummarizedExperiment}"},
    \code{"\link[=ExpressionSet-class]{
    ExpressionSet}"} format}
    \item{col_ctrl}{the columns in the matrix "mat" of the control 
    samples}
    \item{col_test}{the columns in the matrix "mat" of the test
    samples}
    \item{NAME}{Name of genes,or annotation, e.g. WNT4}
    \item{CLID}{Identities of genes,e.g. ENSMUSG00000000001}
    \item{GWEIGHT}{the weight for each gene}
    \item{EWEIGHT}{the weight for each experiment}
}
\details{
    \code{\link[TTMap]{make_matrices}} generates the test matrix 
    and the control matrix in
    the format accepted by 
    \code{\link[TTMap]{control_adjustment}} 
    from a matrix object
}
\value{
    \item{junk}{A list containing $CTRL and $TEST the matrices 
    to impute 
    in \code{\link[TTMap]{control_adjustment}}}
}
\author{
    Rachel Jeitziner}
\seealso{
    \code{\link[TTMap]{control_adjustment}},
    \code{\link[TTMap]{hyperrectangle_deviation_assessment}},
    \code{\link[TTMap]{ttmap_sgn_genes}},
    \code{"\link[=RangedSummarizedExperiment-class]{
    RangedSummarizedExperiment}"}
}
\examples{
    ##--
    ##--
    Aa = 6
    B1 = 3
    B2 = 3
    C0 = 100
    D0 = 10000
    a0 = 4
    b0 = 0.1
    a1 = 6
    b1 = 0.1
    a2 = 2
    b2 = 0.5
    ALPHA = 1
    E = 1
    Pw = 1.1
    Bw = 0
    RA <- matrix(rep(0, Aa * D0), nrow = D0)
    RB1 <- matrix(rep(0, B1 * D0), nrow = D0)
    RB2 <- matrix(rep(0, B2 * D0), nrow = D0)
    RA <- lapply(seq_len(D0 - C0), function(i) rnorm(Aa, 
    mean = a0, sd = sqrt(b0)))
    RA<-do.call(rbind, RA)
    RB1<- lapply(seq_len(D0 - C0), function(i) rnorm(B1, 
    mean = a0, sd = sqrt(b0)))
    RB1 <- do.call(rbind, RB1)
    RB2 <- lapply(seq_len(D0 - C0), function(i) rnorm(B2, 
    mean = a0, sd = sqrt(b0)))
    RB2 <- do.call(rbind, RB2)
    RA_c <- lapply(seq_len(C0), function(i) rnorm(Aa, 
    mean = a0, sd = sqrt(b0)))
    RA_c <- do.call(rbind, RA_c)
    RB1_c <- lapply(seq_len(C0), function(i) rnorm(B1, 
    mean = a1, sd = sqrt(b1)))
    RB1_c <- do.call(rbind, RB1_c)
    RB2_c <- lapply(seq_len(C0), function(i) rnorm(B2, 
    mean = a2, sd = sqrt(b2)))
    RB2_c <- do.call(rbind, RB2_c)
    norm1 <- rbind(RA, RA_c)
    dis <- cbind(rbind(RB1, RB1_c), rbind(RB2, RB2_c))
    colnames(norm1) <- paste("N", seq_len(Aa), sep = "")
    rownames(norm1) <- c(paste("norm", seq_len(D0 - C0), sep = ""),
    paste("diff", seq_len(C0), sep = ""))
    colnames(dis) <- c(paste("B1", seq_len(B1), sep=""),
    paste("B2", seq_len(B2), sep  =""))
    rownames(dis)<-c(paste("norm",
    seq_len(D0 - C0), sep = ""), 
    paste("diff", seq_len(C0), sep = ""))
    the_experiment <- TTMap::make_matrices(cbind(norm1, dis),
    col_ctrl = colnames(norm1),
    col_test = colnames(dis), NAME = rownames(norm1),
    CLID = rownames(norm1))
    ###other example using SummarizedExperiment
    library(airway)
    data(airway)
    airway <- airway[rowSums(assay(airway))>80,]
    assay(airway) <- log(assay(airway)+1,2)
    the_experiment <- TTMap::make_matrices(airway, 
    seq_len(4), seq_len(4) + 4,
    rownames(airway), rownames(airway))
}
