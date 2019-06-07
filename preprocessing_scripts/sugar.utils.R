'%!in%' <-
  function(x, y) {
    !('%in%'(x, y))
  } # Opposite of the dplyr %in% statement.

getFirstN <-
  ## For each plate; find the first/last N elements based on time
  function(Dataframe,
           N = 3,
           first = T) {
    Ns <- Dataframe[0, ] # Clone structure
    plates <- unique(Dataframe$plate)
    for (plate in plates) {
      plateDF <- (Dataframe[Dataframe$plate == plate, ])
      sortBy <- order(plateDF$time)
      if (!first) {
        sortBy <- rev(sortBy)
      }
      topN <- plateDF[sortBy, ][1:N,]
      Ns <- rbind(Ns, topN)
    }
    
    return(Ns)
  }

getLastN <-
  ## Get last N elements based on time
  function(Dataframe, N = 3) {
    return(getFirstN(
      Dataframe = Dataframe,
      N = N,
      first = F
    ))
  }

fence <-
  # Clip values based on upper/lower bounds.
  function(vec, UB = 3, LB = -3)
    pmax(LB, pmin(vec, UB))

fence.matrix <-
  # Like fence but for matrix
  function(mat,
           UB = 3,
           LB = -3,
           removeWT = F) {
    mat.clipped <-
      rbind.data.frame(apply(
        X = mat,
        MARGIN = 2,
        FUN = function(vec)
          fence(vec, UB = UB, LB = LB)
      ))
    rownames(mat.clipped) <- row.names(mat)
    
    if (removeWT) {
      mat.clipped = mat.clipped[!grepl(rownames(mat.clipped), pattern = "^N2"), ]
    }
    
    return(mat.clipped)
  }

mask.by.pvalues <-
  # Mask a matrix by p-values
  # Assumption: Rows are genes, columns are metrics
  function(src,
           mask,
           pvalue.threshold = 0.05,
           fdr.threshold = 0.1,
           use.fdr = F) {
    if (!(
      assertthat::are_equal(rownames(src), rownames(mask)) &
      assertthat::are_equal(colnames(src), colnames(mask))
    )) {
      stop("[ERROR] Row/Cols do not matchs]")
    }
    
    # Prepare to compute mask
    mask.computed = mask
    src.masked = src
    
    if (use.fdr) {
      for (colname in colnames(src)) {
        mask.computed[, colname] <-
          p.adjust(mask[, colname], method = 'fdr')
      }
      src.masked[mask.computed > fdr.threshold] = 0.0
    } else {
      src.masked[mask.computed > pvalue.threshold] = 0.0
    }
    
    return(src.masked)
  }
