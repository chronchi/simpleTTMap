\name{hyperrectangle_deviation_assessment}
\alias{hyperrectangle_deviation_assessment}
\title{
    Calculation of deviation components
}
\description{
    \code{\link[TTMap]{hyperrectangle_deviation_assessment}} function 
    calculates the 
    hyperrectangle deviation
    assessment (HDA) that calculates the 
    deviation components using 
    \code{normal_hda2} which calculates 
    the normal component of the test 
    sample and \code{deviation_hda2}
    which calculates the deviation component.
}
\usage{
    hyperrectangle_deviation_assessment(x, 
    k = dim(x$Normal.mat)[2], dataname, 
    normalname,Org.directory = getwd())
}
\arguments{
    \item{x}{
    output object given back 
    by \code{\link[TTMap]{control_adjustment}}, list}
    \item{k}{
    A factor if not all the lines in 
    the control group should be kept}
    \item{dataname}{
    the name of the project}
    \item{normalname}{
    A name for the corrected control files}
    \item{Org.directory}{
    where the outputs should be saved}
}
\details{
    The function performs the hyperrectangle deviation 
    assessment (HDA)
}
\value{
    Outputs
    \item{Tdis.pcl}{The matrix of the deviation components for each 
    test sample}
    \item{Tnorm.pcl}{The matrix of the normal components for each test 
    sample}
    \item{NormalModel.pcl}{The normal model used}
    Values
    \item{Dc.Dmat}{the deviation component matrix composed of the 
    deviation components of all the samples in the test group}
    \item{m}{the values of the filter function per sample in the test 
    group}
}
\author{
    Rachel Jeitziner}
\seealso{
    \code{\link[TTMap]{control_adjustment}},
    \code{\link[TTMap]{hyperrectangle_deviation_assessment}},
    \code{\link[TTMap]{ttmap_sgn_genes}}}
\examples{
    ##a full example can be found in ttmap_sgn_genes
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
    k = Kprime, dataname = "Effect_of_cancer",
    normalname = "The_healthy_controls");
}
