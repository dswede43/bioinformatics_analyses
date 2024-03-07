#Load read count data
#---
#Function to load read count data and metadata into R

#parameters
#read_count_file: name of the read count data file
#metadata_file: name of metadata for read count data
#factors: names of categorical variables

#load packages
library(dplyr)


#Function
#---
#load data
load_read_count_data = function(read_count_file, metadata_file, factors){
    read_counts = read.csv(read_count_file, header = TRUE, check.names = FALSE, sep = ",", row.names = 1) #gene count matrix

    metadata = read.csv(metadata_file, header = TRUE, sep = ",", row.names = 1) #sample metadata
    metadata = mutate(metadata, across(all_of(factors), as.factor)) #convert variables to factors

    #check for matching sample dimensions and names in both data frames
    if(all(colnames(read_counts) %in% rownames(metadata)) & all(colnames(read_counts) == rownames(metadata)) == TRUE){
        print("Column dimensions and names match.")
    } else{
        print("Column dimensions and names do not match.")
    }
    return(list("read_counts" = read_counts, "metadata" = metadata))
}

