getDateInformationByPlate  <-
  ## Obtain the date experiment data information
  function(Dataframe) {
    ## Return the day part of the plate string
    days = data.frame(
      plate = rownames(Dataframe),
      day = gsub(
        pattern = "_.*",
        replacement = "",
        gsub(
          pattern = ".*:",
          replacement = "",
          rownames(Dataframe)
        )
      )
    )
    return(days)
  }

cleanAlleleNames <-
  ## Get rid of spaces in rownames, return as vector
  function(X)
    gsub(pattern = " ",
         replacement = "",
         x = rownames(X))

clipPlate <-
  ## Remove plate from string.
  function(x)
    gsub(pattern = ":.*", replacement = "", x)

simpleSquash <-
  ##  Aggregate values `groupVector` by function `f` 
  function(x, groupVector, f=mean) {
    squashed <- data.frame()
    for (group in unique(groupVector) ){
      row <- data.frame(
        apply(X = x[grepl(group, groupVector, fixed=T), ], 
              MARGIN = 2, 
              FUN=function(y) f(y, na.rm = T)
        )
      )
      row <- t(row)
      row.names(row) <- group
      squashed <- rbind(squashed, row)
    }
    return(squashed)
  }


meltClusters <- 
  # Creates something similar to a molten dataframe from clustering
  # Assumes the result of pvclust::pvpick as input X
  # Returns a data frame.
  function(X){
    groups <- data.frame(group=character(0),
                         gene=character(0))
    
    for (i in 1:length(X)){
      for (gene in X[[i]]){
        row <- data.frame(group=as.character(i), gene=gene)
        groups <- rbind(groups, row)
      }
    }
    return(groups)
  }

