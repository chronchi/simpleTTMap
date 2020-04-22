\name{control_adjustment}
\alias{control_adjustment}
\title{
    Calculates a corrected control group, discovers outliers in it.
}
\description{
    \code{\link[TTMap]{control_adjustment}} function 
    finds outliers in the control group and 
    removes them
}
\usage{
    control_adjustment(normal.pcl, tumor.pcl, normalname, dataname, 
    org.directory = "", A = 1, e = 0, meth = 0, P = 1.1, B = 0)
}
\arguments{
    \item{normal.pcl}{
    the control matrix with annotation as obtained by $CTRL from 
    \code{\link[TTMap]{make_matrices}}}
    \item{tumor.pcl}{
    the disease/test data matrix with annotation as obtained by $TEST 
    from \code{\link[TTMap]{make_matrices}}} 
    \item{normalname}{
    A name for the corrected control files}
    \item{dataname}{
    the name of the project}
    \item{org.directory}{
    where the outputs should be saved}
    \item{A}{
    integer if A=0 then the difference to the median is calculated 
    otherwise the difference to the mean.
    }
    \item{e}{
    integer giving how far to the median an outlier is at least}
    \item{meth}{
    value or method that defines how to replace outliers, default is 
    set to replace by the median
    }
    \item{P}{
    if more than P percent of features are outliers the feature is 
    removed, by default all are kept}
    \item{B}{
    Batch vector a vector for normal and test samples with a same 
    number corresponding to a same batch}
}
\details{
    \code{\link[TTMap]{control_adjustment}} calculates a 
    corrected control group, discovers 
    outliers in it.
}
\value{
    Several files are created
    \item{paste(org.directory,normalname,".normMesh",sep = "")}{The 
    normal matrix with only common features with the test matrix. This 
    file is only created if the two have different rows}
    \item{paste(org.directory,dataname,".normMesh",sep = "")}{The 
    test matrix with only common features with the normal matrix. This 
    file is only created if the two have different rows.}
    \item{mean_vs_variance.pdf}{A pdf showing a plot of the mean (X 
    axis) against the variances (Y axis) of each feature}
    \item{mean_vs_variance_after_correction.pdf}{A pdf showing a plot 
    of the mean (X axis) against the variances (Y axis) of each feature
    after correction of the control group}
    \item{na_numbers_per_row.txt}{number of outliers per row}
    \item{na_numbers_per_col.txt}{number of outliers per column}
    And values of ttmap_part1_ctrl_adj
    \item{e}{Selected criteria for what is an outlier}
    \item{tag.pcl}{Annotation of features, ID of features and weight}
    \item{Normal.mat}{The control matrix without annotation and only 
    with the common rows with Disease.mat}
    \item{Disease.mat}{The test/disease matrix without annotation and 
    only with the common rows with Disease.mat}
    \item{flat.Nmat}{A list $mat being the corrected control matrix $m 
    a record of the different numbers of removed genes per sample}
    \item{record}{numbers recording the number of columns in 
    Disease.mat and Normal.mat}
    \item{B}{The batch vector B introduced in the begining}
    \item{U1}{The different batches in Normal.mat}
    \item{U2}{The different batches in Disease.mat}
}
\author{
    Rachel Jeitziner}
\seealso{
    \code{\link[TTMap]{hyperrectangle_deviation_assessment}},\code{
    \link[TTMap]{ttmap}}
    \code{
    \link[TTMap]{ttmap_sgn_genes}}}
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
}
