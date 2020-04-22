read_pcl <-
function(filename, na.type = "", Nrows = -1, 
Comment.char = "",...) {
    x.df <- read.table(paste(filename,
    "pcl",sep = "."), header = TRUE, sep = "\t", 
    quote = "\"", na.strings = na.type, skip = 0,
    nrows = Nrows, comment.char = Comment.char,...);
    x.df$CLID <- as.vector(x.df$CLID, mode = "character");
    x.df$NAME <- as.vector(x.df$NAME, mode = "character");
    rownames(x.df) <- x.df$CLID;
    return(x.df)
}
