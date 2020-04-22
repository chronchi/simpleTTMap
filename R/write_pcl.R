write_pcl <- function(df, dataname, output_directory = ""){
    
    # path to save the dataframe 
    dir.address <- file.path(output_directory, paste(dataname, ".pcl", sep = ""))
    
    # save dataframe and return it  
    X <- write.table(df, file = dir.address,
                     append = FALSE, quote = FALSE, 
                     sep = "\t", eol = "\n", na = "",
                     row.names = FALSE, col.names = TRUE)

    return(X)
}
