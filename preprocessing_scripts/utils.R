library("ggplot2")
library("ggrepel")
library("plyr")
library("dplyr")
library("pheatmap")
library("reshape2")
library("iheatmapr")
library("pvclust")
library("dendextend")
library("parallel")
library("RColorBrewer")
library("Rtsne")
library("plot3D")
library("Cairo")
library("dendsort")

source("preprocessing_scripts/constants.utils.R")
source("preprocessing_scripts/sugar.utils.R")
source("preprocessing_scripts/math.utils.R")
source("preprocessing_scripts/normalization.utils.R")
source("preprocessing_scripts/io.utils.R")
source("preprocessing_scripts/wrangle.utils.R")
source("preprocessing_scripts/plotting.utils.R")

computeProbabilities <-
  ## Compute probabilities for Not Rerversing (NR) and Started Reversing(SR)
  ## Expects a dataframe with
  ## -nWormNR (n of worms not reversing)
  ## -nWormSR (n of worms who started reversing)
  ## -nWormAR (n of worms Already Reversing)
  ## -DurationAvg (average reversal duration)
  ## -DistAvg (average reversal distance)
  function(Dataframe) {
    denominator <- Dataframe$nWormNR + Dataframe$nWormSR
    
    #ratioAR <- Dataframe$nWormAR / denominator # Not sure how one would compute that.
    ratioNR <- Dataframe$nWormNR / denominator
    ratioSR <- Dataframe$nWormSR / denominator
    
    # Merge ratio to dataframe
    return(cbind(
      Dataframe,
      data.frame(
        nWormAR = Dataframe$nWormAR,
        # Let's just return the number of AR.
        probNR = ratioNR,
        probSR = ratioSR,
        revDuration = Dataframe$DurationAvg,
        revDistance = Dataframe$DistAvg
      )
    ))
  }

dropRedundant <-
  ## Drop features that are considered redundant in the analysis.
  function(X) {
    return(X[, setdiff(names(X), REDUNDANT)])
  }

collapsePlatesToMatrix <-
  # Aggregate dataframe by function FUN per variable x plates.
  function(Dataframe, FUN = mean_rm.na) {
    labeled_ <- unique(paste0(Dataframe$line, ":", Dataframe$plate))
    plateMtrx <- aggregateLine(smeltPlate(Dataframe), FUN = FUN)
    
    labeled <- c()
    for (pn in 1:nrow(plateMtrx)) {
      currentlabel <- rownames(plateMtrx)[pn]
      labeled <-
        c(labeled, labeled_[grepl(pattern = currentlabel, x = labeled_)])
      
      if (is.null(currentlabel)) {
        print(paste("removing", pn))
      }
    }
    
    rownames(plateMtrx) <- labeled
    return(plateMtrx)
  }

applyColMeansByPlate <-
  ## Compute column mean per plate on continuous features only
  function(dataframe) {
    DISCRETE_IDX <- match(DISCRETE, names(dataframe))
    DISCRETE_IDX <- DISCRETE_IDX[!is.na(DISCRETE_IDX)]
    dataframe.discrete <- dataframe[, DISCRETE_IDX]
    dataframe.continuous <- dataframe[, -DISCRETE_IDX]
    
    meanDataframe <- data.frame(stringsAsFactors = FALSE)
    for (key in unique(dataframe$plate)) {
      discrete <- dataframe.discrete[(dataframe.discrete$plate == key) ,]
      ss <-
        subset(dataframe.continuous,
               subset = (dataframe.discrete$plate == key))
      
      # Doing things a bit differently here. We can't use rowname since we have duplicates.
      cm <- t(data.frame(colMeans(ss, na.rm = TRUE)))
      geneRowMean <-
        cbind.data.frame(unique(discrete), cm) # In case we'd want to merge the feature back in.
      
      meanDataframe <- rbind.data.frame(meanDataframe, geneRowMean)
    }
    
    return(meanDataframe)
  }

collapseToMatrix <-
  ## Aggregate dataframe by function FUN per variable ~ alleles.
  function(Dataframe, FUN = mean_rm.na) {
    return(aggregateLine(smeltLine(Dataframe), FUN = FUN))
  }

smeltLine <-
  ## Melt dataframe into a (line, variable, value) table.
  function(X) {
    X_ <- cbind(X$line, getContinuous(X, also = COMPUTED))
    X_ <- melt(X_)
    names(X_) <- c("line", "variable", "value")
    return(X_)
  }

smeltPlate <-
  ## Melt dataframe into a (plate, variable, value) table.
  function(X) {
    X_ <- cbind(X$plate, getContinuous(X, also = COMPUTED))
    X_ <- melt(X_)
    names(X_)[1] <-
      c("line") # Fixme. This is just so other functions expecting 'line' as the primary ID won't complain.
    return(X_[, c("line", "variable", "value")])
  }

getContinuous <-
  # Return the Discrete elements of the dataframe
  function(X, also = NULL) {
    if (!is.null(also)) {
      DISCRETE_ <- setdiff(DISCRETE, also)
    } else {
      DISCRETE_ <- DISCRETE
    }
    
    CONTINUOUS_IDX <- match(DISCRETE_, names(X))
    CONTINUOUS_IDX <- CONTINUOUS_IDX[!is.na(CONTINUOUS_IDX)]
    
    if (length(CONTINUOUS_IDX) < 1) {
      return(X)
    }
    return(X[, -CONTINUOUS_IDX])
    
  }

sdOfTwoSlices <-
  # Given two slices of time, compute the SD of their differences
  function(a, b, VAR = NULL) {
    FUN = sd_rm.na
    
    a <- a[, c("line", "plate", VAR)]
    b <- b[, c("line", "plate", VAR)]
    a.melted <- melt(a)
    b.melted <- melt(b)
    
    a.melted.agg <-
      dcast(a.melted,
            value.var = "value",
            fun.agg = mean_rm.na,
            plate ~ variable) # Agregating on plate x variable
    names(a.melted.agg) <- c("plate", "a")
    b.melted.agg <-
      dcast(b.melted,
            value.var = "value",
            fun.agg = mean_rm.na,
            plate ~ variable) # Agregating on plate x variable
    names(b.melted.agg) <- c("plate", "b")
    
    delta <-
      merge(a.melted.agg,
            b.melted.agg,
            by.x = "plate",
            by.y = "plate")
    a.melted$delta <- delta[, "a"] - delta[, "b"]
    
    ds <-
      dcast(a.melted[, c("plate", "variable", "line", "delta")],
            value.var = "delta",
            fun.agg = FUN,
            line ~ variable) # Agregating on line
    
    return(ds[, VAR])
  }

aggregateLine <-
  # Helper method to create an aggregated matrix of value "variable" using funcion "FUN".
  function(X, variable = "value", FUN = mean_rm.na) {
    b <-
      dcast(X, value.var = variable, fun.agg = FUN, line ~ variable) # Agregating on line
    bnames <- b[, 1] # Names
    mtrx <- as.matrix(b[, -1])
    row.names(mtrx) <- bnames
    return(mtrx)
  }
