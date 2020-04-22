ttmap_sgn_genes <- function(ttmap_part2_gtlmap, 
                            ttmap_part1_hda, 
                            ttmap_part1_ctrl_adj, 
                            c, 
                            n = 2, 
                            a = 0, 
                            filename = "TEST2",
                            annot = ttmap_part1_ctrl_adj$tag.pcl, 
                            col = "NAME",
                            output_directory = ".", 
                            Relaxed = 1) {

    # create subdirectories to save the significant genes
    sub_directories = c("all", "mid1", "mid2", "high", "low")
    for (sub_dir in sub_directories) {
        dir.create( file.path( output_directory, sub_dir ), recursive = TRUE)     
    } 
 
    if(Relaxed != 1) {
        ttmap_sgn_genes_inter <- ttmap_sgn_genes_inter2
    }
    
    junk_low <- list(x = 0)

    # node index 
    node_index = 1
    
    # all tiers
    tiers = list('all', 'low_map', 'mid1_map', 'mid2_map', 'high_map')

    for (tier in tiers) {

        tier_gtlmap = ttmap_part2_gtlmap[[tier]]

        for(i in seq_len(length(tier_gtlmap))){

            samples_in_cluster = tier_gtlmap[[i]]

            A <- ttmap_sgn_genes_inter(samples_in_cluster,
                                       ttmap_part1_hda, alpha = a)
            
            A <- as.data.frame( as.matrix(A) )
            
            if (dim(A)[2] == 1){ 
                colnames(A) <- samples_in_cluster
            }

            # A <- cbind(annot[rownames(A), col], A)
            
            tier_strip = str_remove(tier, '_.*')
            
            name_to_save = file.path(output_directory, 
                                     tier_strip,
                                     paste( toString(node_index),
                                        '.txt', sep=''
                                     )                  
                                    )
            write.table(A, 
                        file = name_to_save, 
                        quote = FALSE, sep = "\t", 
                        row.names = TRUE, col.names = NA)

            node_index = node_index + 1
        }
    }

}
