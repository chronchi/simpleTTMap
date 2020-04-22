\name{ttmap}
\alias{ttmap}
\title{
    Visualisation of the clustering
}
\description{
    Enables a quick view on the groups in the dataset (globally) and 
    how locally they differ.
}
\usage{
    ttmap(ttmap_part1_hda, m1, 
    select = row.names(ttmap_part1_hda$Dc.Dmat), 
    ddd, e, filename = "TEST", n = 3, ad = 0, bd = 0, piq = 1, 
    dd = generate_mismatch_distance(ttmap_part1_hda = ttmap_part1_hda, 
    select = select), mean_value_m1 = "N", ni = 2)
}
    %- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ttmap_part1_hda}{
    list output of \code{\link[TTMap]{hyperrectangle_deviation_assessment}}}
    \item{m1}{
    either a user imputed vector whose names are the names of the 
    samples with addition of .Dis. or by default it is the amount of 
    deviation}
    \item{select}{
    Should all the features (default) or only a sublist be considered 
    to calculate the distance}
    \item{ddd}{
    Annotation matrix with rownames the different sample names with 
    addition of .Dis. There can be as many columns as wanted, but only 
    the column n will be selected to annotated the clusters}
    \item{e}{
    integer parameter defining under which value two samples are 
    considered to be close}
    \item{filename}{
    Name for the description file annotating the clusters}
    \item{n}{ The column to be considered to annotate the clusters
    }
    \item{ad}{
    if ad!=0 then the clusters on the output picture will not be 
    annotated}
    \item{bd}{
    if different than 0 (default), the output will be without outliers 
    of the test data set (clusters composed of only "piq" element)}
    \item{piq}{
    parameter used to determine what small clusters are, see bd}
    \item{dd}{
    the distance matrix to be used}
    \item{mean_value_m1}{
    if == "N" the average of the values in m1 divided by the number of 
    the samples are put into the legend (by default represents the 
    average of the samples in a cluster of the mean-deviation of the 
    features) otherwise it will show the average value of the values in 
    m1 (is useful for instance if m1 represents the age of the 
    samples) }
    \item{ni}{
    The column to consider to annotate the samples (is put into 
    parenthesis) for the description file}
}
\value{
    \item{all}{the clusters in the overall group}
    \item{low}{the clusters in the lower quartile group}
    \item{mid1}{the clusters in the first middle quartile group}
    \item{mid2}{the clusters in the second middle quartile group}
    \item{high}{the clusters in the higher quartile group}
}
\details{
    Is the Two-tiers Mapper function. The output is an interactive 
    image of the clusters in the different layers.}
\author{
    Rachel Jeitziner}
\seealso{
    \code{\link[TTMap]{control_adjustment}},
    \code{\link[TTMap]{hyperrectangle_deviation_assessment}},
    \code{\link[TTMap]{ttmap_sgn_genes}}}
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
}
