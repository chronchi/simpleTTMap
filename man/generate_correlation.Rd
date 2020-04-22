\name{generate_correlation}
\alias{generate_correlation}
\alias{generate_euclidean}
\alias{generate_mismatch_distance}
\alias{generate_p_val_correlation}
\alias{generate_single_cell_complete_mismatch}
\alias{generate_single_cell_mismatch_with_parameter}
\title{
    Generates different distance matrices
}
\description{
    Single cell complete mismatch distance, single cell complete 
    mismatch distance with a parameter of cutoff, mismatch distance, 
    correlation distance, p-value of correlation test distance and 
    euclidean distance.
}
\usage{
    generate_single_cell_complete_mismatch(ttmap_part1_hda, 
    select, alpha = 1)
    generate_single_cell_mismatch_with_parameter(ttmap_part1_hda,
    select, alpha = 1)
    generate_correlation(ttmap_part1_hda, select)
    generate_euclidean(ttmap_part1_hda, select)
    generate_mismatch_distance(ttmap_part1_hda, select, alpha = 1)
    generate_p_val_correlation(ttmap_part1_hda, select)
}
\arguments{
    \item{ttmap_part1_hda}{
    an object given back 
    by \code{\link[TTMap]{hyperrectangle_deviation_assessment}}
    }
    \item{select}{
    A sublist of rownames of ttmap_part1_hda$Dc.Dmat
    }
    \item{alpha}{
    A real number corresponding to a cutoff
    }
}
\details{
    If one is interested only in clustering samples according to a list 
    of genes belonging to a certain pathway, then this list is provided 
    to the parameter select. 
    Alpha is a cutoff for deviations that should be considered as noise, 
    for gene expression data such as normalised RNA-seq or microarrays 
    for instance a cutoff of 1, corresponding to a two fold change is 
    being chosen.
}
\value{
    Distance matrix
}
\author{
    Rachel Jeitziner}
\examples{
    ttmap_part1_hda <- list()
    ttmap_part1_hda$Dc.Dmat <- matrix(c(-1, 2, 0, -4, 5, 6), nrow = 2)
    rownames(ttmap_part1_hda$Dc.Dmat) <- c("Gene1", "Gene2")
    colnames(ttmap_part1_hda$Dc.Dmat) <- c("A", "B", "C")
    dd <- TTMap::generate_mismatch_distance(ttmap_part1_hda, select = 
    rownames(ttmap_part1_hda$Dc.Dmat))
    dd <- TTMap::generate_euclidean(ttmap_part1_hda, select = 
    rownames(ttmap_part1_hda$Dc.Dmat))
}
