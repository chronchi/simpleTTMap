\name{calcul_e}
\alias{calcul_e}
\alias{calcul_e_single}
\title{
    Calculation of the value of epsilon}
\description{
    Calculation of the value of epsilon}
\usage{
    calcul_e(dd5, pvalcutoff = 0.95, tt1, alpha = 1, S = 
    colnames(tt1$Normal.mat))
    calcul_e_single(dd5, pvalcutoff = 0.95, tt1, alpha = 1, S = 
    colnames(tt1$Normal.mat))
}
\arguments{
    \item{dd5}{
    distance matrix as created by \code{generate_mismatch_distance}
    }
    \item{pvalcutoff}{
    cutoff of 0.05 percent (default) or less}
    \item{tt1}{
    output of \code{\link[TTMap]{control_adjustment}}}
    \item{alpha}{
    a cutoff value for the FC between the group of control and the 
    disease group}
    \item{S}{
    subset of columns to be considered}
}
\value{
    \item{al}{number representing the cutoff to choose for the 
    relatedness with dd5}
}
\author{
    Rachel Jeitziner}
\seealso{
    \code{\link[TTMap]{control_adjustment}},
    \code{\link[TTMap]{hyperrectangle_deviation_assessment}},
    \code{\link[TTMap]{ttmap_sgn_genes}},
    \code{\link{generate_mismatch_distance}}
}
\examples{
    ##--
    library(airway)
    data(airway)
    airway <- airway[rowSums(assay(airway))>80,]
    assay(airway) <- log(assay(airway)+1,2)
    ALPHA <- 1
    the_experiment <- TTMap::make_matrices(airway,
    seq_len(4), seq_len(4) + 4,
    rownames(airway), rownames(airway))
    TTMAP_part1prime <-TTMap::control_adjustment(
    normal.pcl = the_experiment$CTRL,
    tumor.pcl = the_experiment$TEST, 
    normalname = "The_healthy_controls", 
    dataname = "Effect_of_cancer", 
    org.directory = tempdir(), e = 0, P = 1.1, B = 0);
    Kprime <- 4;
    TTMAP_part1_hda <-
    TTMap::hyperrectangle_deviation_assessment(x = 
    TTMAP_part1prime,
    k = Kprime,dataname = "Effect_of_cancer",
    normalname = "The_healthy_controls");
    annot <- c(paste(colnames(
    the_experiment$TEST[,-(seq_len(3))]), "Dis", sep = "."),
    paste(colnames(the_experiment$CTRL[, 
    -seq_len(3)]), "Dis", sep = "."))
    dd5_sgn_only <-TTMap::generate_mismatch_distance(
    TTMAP_part1_hda,
    select=rownames(TTMAP_part1_hda$Dc.Dmat), alpha = ALPHA)
    e <- TTMap::calcul_e(dd5_sgn_only, 0.95, TTMAP_part1prime, 1)
}
