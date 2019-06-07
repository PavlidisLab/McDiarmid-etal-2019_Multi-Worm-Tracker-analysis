writeObject <- function(Dataframe,
                        filename = NULL,
                        path = NULL,
                        ...) {
  if (is.null(filename)) {
    filename <-  deparse(substitute(Dataframe))
  }
  
  if (is.null(path)) {
    path <- config::get("data")
  }
  
  filepath <- paste0(path, filename)
  write.table(Dataframe, file = filepath, ...)
}

readObjectAsMatrix <- function(filename, path = NULL) {
  if (is.null(path)) {
    path <- config::get("data")
  }
  return(as.matrix(read.table(paste0(path, filename))))
}

readObjectAsTable <- function(filename, path = NULL) {
  if (is.null(path)) {
    path <- config::get("data")
  }
  return(read.table(paste0(path, filename)))
}

readMatrix <- function(path){
  # Doesn't assume the path starts with "data/"
  return(as.matrix(read.table(path)))
}
