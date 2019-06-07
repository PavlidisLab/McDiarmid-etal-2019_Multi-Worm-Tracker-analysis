scaleCols <-
  # Apply built in scale per column
  function(X,
           center.bool = T,
           scale.bool = T) {
    for (feature in colnames(X)) {
      X[, feature] = scale(X[, feature], center = center.bool, scale = scale.bool)
    }
    return(X)
  }