# extract pvalues from the dataframe. the default key is p.value
extract_pvalues <- function(dataframe){
  # extract keys
  df_names <- names(dataframe)
  # create empty array
  pvalues <- c()
  tvalues <- c()
  # iterate over keys to extract pvalues
  for (n in df_names){
    pvalues <- c(pvalues, dataframe[[n]]$p.value)
    tvalues <- c(tvalues, dataframe[[n]]$statistic)
  }
  stats_ <- list(pvalues, tvalues)
  return(stats_)
}

# create dataframe containing p.values mean deviation component, 
# adjusted p-values and gene names. The input is a table wit t.test results
create_dataframe <- function(ttest_results){
  # first extract p-values
  stats_ <- extract_pvalues(ttest_results)
  # adjust pvalues
  adjus_pvalues <- p.adjust(stats_[[1]], method='BH')
  # create final dataframe with gene name, deviation component, p-value and 
  # adjusted p-value
  genes <- names(ttest_results)
  mean_dev_com <- c()
  for(n in genes){
    mean_dev_com <- c(mean_dev_com, as.numeric(ttest_results[[n]]$estimate))
  }
  final_dataframe = data.frame(Gene=genes, 
                               logFC = mean_dev_com,
                               t = stats_[[2]],
                               P.Value = stats_[[1]],
                               adj.P.Val = adjus_pvalues)
  return(final_dataframe)
}

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
            
            # calculate the p-values associated for each gene. 
            # the null hypothesis is that the genes don't deviate 
            # from the control group, i.e., their deviation components
            # are zero. 
            
            #t_tests <- apply(as.matrix(A), 1, t.test)  
            
            #t_results <- create_dataframe(t_tests)

            # append the p value and adjusted p value columns to the original
            # dataframe
            #A <- cbind('adj_p_value'=t_results$adj.P.Val, A)
            #A <- cbind('p-value'=t_results$P.Value, A)
            
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
