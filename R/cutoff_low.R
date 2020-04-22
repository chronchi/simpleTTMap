cutoff_low <-
function(low_map, para = 1, text = "low"){
    if(dim(low_map)[1] == 1){if(low_map == 0){
        print(paste(paste("There were not enough", text, sep = " "),
        "classified samples", sep=" "))}}
    else{s <- apply(low_map, 1, sum)
        low_map <- low_map[s > para,]}
    return(low_map);
}
