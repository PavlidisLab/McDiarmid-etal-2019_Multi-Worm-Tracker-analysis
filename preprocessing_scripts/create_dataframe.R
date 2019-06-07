library(here)
setwd(here::here())


morphological <- list.files(config::get("datasets"), 
                            recursive = T, 
                            full.names = T, 
                            pattern = ".*dat")

morphological_features <- c(c("plate", "line", "gene"), names(read.csv(config::get("morphological-features"), header = T )))

stimuli_induced <- list.files(config::get("datasets"), 
                              recursive = T, 
                              full.names = T,
                              pattern = ".*rev")

stimuli_features <- c(c("plate", "line", "gene"), names(read.csv(config::get("stimuli-features"), header = T )))

## Check that a matching number of .dat and .rev files are loaded.
if (!assertthat::are_equal(length(morphological), length(stimuli_induced))) { stop("[ERROR] Mismatch in the number of `.dat`/`.rev` files.") }

# IMPORTANT: Assuming the file names are structure as LINE/PLATE/FILE.{dat,rev}, as output by Choreography.
plates <- sapply(strsplit(morphological, "/"), function(x) x[length(x) - 1])
unique(plates)
length(unique(plates))

lines <- sapply(strsplit(morphological, "/"),function(x) x[length(x) - 2])
unique(lines)
length(unique(lines))

# Gene is the line split on space.
genes <- sub("\\(.*", "",
          sub(" .*", "", lines)
)
names(table(unique(genes)))
length(unique(genes))

## Check that the plates are not duplicated
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

  return(all_phenotypes)  
}

morphologicalDF <- add_lines(morphological)
names(morphologicalDF) <- morphological_features
head(morphologicalDF)
tail(morphologicalDF)
write.table(morphologicalDF, file = paste0(config::get('data'), 'morphologicalDF.tsv'), sep = "\t")

stimuliDF <- add_lines(stimuli_induced)
names(stimuliDF) <- stimuli_features
head(stimuliDF)
tail(stimuliDF)
write.table(stimuliDF, file = paste0(config::get('data'), 'stimuliDF.tsv'), sep = "\t")
