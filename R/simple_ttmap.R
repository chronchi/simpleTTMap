# Wrap all functions from TTMap package into a single function. 
# It should return the clusters, normalized control matrix and deviation
# components. 

simple_ttmap <- function( dataset, 
                          outlier_parameter = 1.0,
                          output_directory  = '.',
                          alpha             = 1.0, 
                          batches           = c(),
                          test_samples      = c(),
                          dataname          = 'ttmap',
                          plot_graph        = FALSE 
                        ){

    # define the experiment matrix. First we separate the control group to the rest
    name_samples = colnames(dataset)
    ctrl_indices = grep("ctrl", tolower(name_samples))
    
    if (length(test_samples) == 0){
        test_indices = setdiff(seq_len(length(name_samples)), ctrl_indices)
    }
    else {
        test_indices = grep(paste(tolower(test_samples), sep = "|"), tolower(name_samples))    
    }

    # create experiment matrix
    experiment <- make_matrices(dataset, 
                                ctrl_indices, 
                                test_indices, 
                                NAME=rownames(dataset), 
                                CLID=rownames(dataset)
                               )

    # adjustment of the control group by batches.
    if(length(batches) == 0){
        batch_numbers = 0
    } 
    else {
        # find the batch correspondents in the ctrl and test groups
        number_of_batches <- length(batches)
        counter = 0
        batch_numbers = c()
        
        batches_correspondence = seq_len(length(batches))-1
        names(batches_correspondence) = batches
        
        # first control
        for (ctrl_indice in ctrl_indices){
            # find the batch for the specific sample
            batch_sample <- match(1, lapply(tolower(batches), grep, tolower(name_samples[ctrl_indice])))
            batch_numbers <- c(batch_numbers,
                               batches_correspondence[[batch_sample]])
        }
        for (test_indice in test_indices){
            # find the batch for the specific sample
            batch_sample <- match(1, lapply(tolower(batches), grep, tolower(name_samples[test_indice])))
            batch_numbers <- c(batch_numbers,
                               batches_correspondence[[batch_sample]])
        }
    }
 
 
    ttmap_ctrl_adjm <-control_adjustment(normal.pcl = experiment$CTRL,
                                         tumor.pcl = experiment$TEST,
                                         normalname = "ctrl", 
                                         dataname = dataname,
                                         output_directory = output_directory,
                                         B = batch_numbers,
                                         e = outlier_parameter,
                                        )
    
    # hyperrectangle deviation assessment 
    ttmap_hda <- hyperrectangle_deviation_assessment(
                                x = ttmap_ctrl_adjm,
                                k = dim(ttmap_ctrl_adjm$Normal.mat)[2],
                                dataname = dataname, 
                                normalname = "ctrl",
                                output_directory = output_directory
                )

    # annotation file for plotting ttmap
    annot <- c(paste(colnames(experiment$TEST[,-seq_len(3)]), "Dis", sep="."),
               paste(colnames(experiment$CTRL[,-seq_len(3)]), "Dis", sep="."))
    annot <- cbind(annot,annot)
    
    # create a third column with names for mapper 
    annot <- cbind(annot, paste("S", seq_len(dim(experiment$CTRL)[2] + dim(experiment$TEST)[2] - 6), sep=""))
    rownames(annot) <- annot[,1]
    colnames(annot)[3] = 'Sample'
    write.csv(annot, file.path(output_directory, 'annot.txt'))

    # calculate mismatch distance among samples. Here we use all genes
    distance_samples <- generate_mismatch_distance(ttmap_hda,
                                                   select=rownames(ttmap_hda$Dc.Dmat), 
                                                   alpha = alpha)

    # calculate clusters and possibly plot
    ttmap_gtlmap <- ttmap(ttmap_hda,
                          ttmap_hda$m,
                          select = rownames(ttmap_hda$Dc.Dmat),
                          annot,
                          e = calcul_e(distance_samples, 0.95, ttmap_ctrl_adjm, alpha),
                          filename = dataname,
                          output_directory = output_directory,
                          n = 3, 
                          dd = distance_samples, 
                          bd = 0,
                          ad = 0,
                          plot_graph = plot_graph
                         )

    # save cluster information 
    output_significant_genes = file.path(output_directory, "clusters")
    ttmap_sgn_genes(ttmap_gtlmap, 
                    ttmap_hda,
                    ttmap_ctrl_adjm, 
                    annot, 
                    n = 3, 
                    filename = "test_samples", 
                    annot = ttmap_ctrl_adjm$tag.pcl, 
                    col = "NAME", 
                    output_directory = output_significant_genes,
                    Relaxed = 1, 
                    a = alpha)

    name_samples_test = experiment$TEST
    
    # create list to store all outputs and return 
    ttmap_all_results <- list(ttmap_ctrl   = ttmap_ctrl_adjm, 
                              ttmap_hda    = ttmap_hda,
                              ttmap_gtlmap = ttmap_gtlmap,
                              samples_name = name_samples_test 
                             ) 

    return(ttmap_all_results) 
}
