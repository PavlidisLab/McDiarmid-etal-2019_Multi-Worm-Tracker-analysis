library(dplyr)

setwd("/space/grp/mbelmadani/worms/")

morphological <- list.files("/space/grp/mbelmadani/worms/NAll_dat_files/", 
                            recursive = T, 
                            full.names = T)

plates <- sapply(strsplit(morphological, "/"), function(x) x[10])
unique(plates)

lines <- sapply(strsplit(morphological, "/"),function(x) x[9])
unique(lines)

genes <- sub(" .*", "", lines)
unique(genes)

if ((length(plates) == length(morphological) & length(stimuli_induced)) == length(stimuli_induced)) { stop("Error; plates should be unique.") }

add_lines <- function(fs) {
  count_lines <- 1
  all_phenotypes <- data.frame()
  
  for (f in fs) {
    phenotypes <- read.table(file = f)  
    
    phenotypes <- cbind(genes[count_lines], phenotypes) 
    phenotypes <- cbind(lines[count_lines], phenotypes) 
    phenotypes <- cbind(plates[count_lines], phenotypes) 
    
    count_lines <- count_lines + 1
    all_phenotypes <- rbind(all_phenotypes, phenotypes)
  }
  #all_phenotypes <- sapply(X = add_lines(fs), FUN = cbind) # Load plates
  
  return(all_phenotypes)  
}

morphologicalDF <- add_lines(morphological)
names(morphologicalDF) <- morphological_features
head(morphologicalDF)
tail(morphologicalDF)
write.table(morphologicalDF, file = "morphologicalDF-NALL.tsv", sep = "\t") # TODO: Renamed with -NALL to avoid collision.

stimuliDF <- add_lines(stimuli_induced)
names(stimuliDF) <- stimuli_features
head(stimuliDF)
tail(stimuliDF)
write.table(stimuliDF, file = "stimuliDF-NALL.tsv", sep = "\t")  # TODO: Renamed with -NALL to avoid collision.
