\name{ttmap_sgn_genes}
\alias{ttmap_sgn_genes}
\alias{ttmap_sgn_genes_inter2}
\alias{ttmap_sgn_genes_inter}
\title{
    Gives a list of associated genes per cluster
}
\description{
    \code{\link[TTMap]{ttmap_sgn_genes}} function
}
\usage{   
    ttmap_sgn_genes(ttmap_part2_gtlmap, ttmap_part1_hda,
    ttmap_part1_ctrl_adj, c, n = 2, a = 0, 
    filename = "TEST2", annot = ttmap_part1_ctrl_adj$tag.pcl, 
    col = "NAME", path = getwd(), Relaxed = 1)
    ttmap_sgn_genes_inter2(q, ttmap_part1_hda, alpha = 0)
    ttmap_sgn_genes_inter(q, ttmap_part1_hda, alpha = 0)
}
\arguments{
    \item{ttmap_part2_gtlmap}{
    output of \code{\link[TTMap]{ttmap}}}
    \item{ttmap_part1_hda}{
    output of \code{\link[TTMap]{hyperrectangle_deviation_assessment}}}
    \item{ttmap_part1_ctrl_adj}{
    output of \code{\link[TTMap]{control_adjustment}}}
    \item{c}{annotation file of the samples
    }
    \item{n}{
    column to give the name to the cluster}
    \item{a}{cutoff to be considered different than noise
    }
    \item{filename}{
    Name of the files}
    \item{annot}{
    annotation file}
    \item{col}{
    which column should be considered to annotate the features}
    \item{path}{
    where to put the output files}
    \item{Relaxed}{
    If Relaxed then one allows sample to be as the 
    control and for all 
    the others in one cluster to be 
    going in the same direction (more 
    than alpha) otherwise all 
    the features must be deviating to be 
    considered a significant feature}
    \item{q}{
    The sample in one cluster}
    \item{alpha}{
    cutoff to be considered different than noise inherited by a
    }
}
\details{
    Is giving per cluster the features 
    that vary in the same direction
}
\value{
    generates a file per cluster of significant 
    features with an annotation
}
\author{
    Rachel Jeitziner}
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
    the_experiment$TEST[,-(seq_len(3))]),"Dis", sep = "."),
    paste(colnames(the_experiment$CTRL[, 
    -seq_len(3)]), "Dis", sep = "."))
    annot <- cbind(annot, annot)
    rownames(annot)<-annot[, 1]
    dd5_sgn_only <-TTMap::generate_mismatch_distance(
    TTMAP_part1_hda,
    select=rownames(TTMAP_part1_hda$Dc.Dmat), alpha = ALPHA)
    TTMAP_part2 <-
    TTMap::ttmap(TTMAP_part1_hda, TTMAP_part1_hda$m,
    select = rownames(TTMAP_part1_hda$Dc.Dmat), annot,
    e = TTMap::calcul_e(dd5_sgn_only, 0.95, TTMAP_part1prime, 1), 
    filename = "first_comparison", n =  1, dd = dd5_sgn_only)
    TTMap::ttmap_sgn_genes(TTMAP_part2, TTMAP_part1_hda, 
    TTMAP_part1prime, annot,
    n = 2, a = 1, filename = "first_list_of_genes",
    annot = TTMAP_part1prime$tag.pcl, col = "NAME", 
    path = getwd(), Relaxed = 1)
}
