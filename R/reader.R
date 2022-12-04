# use data.table because its faster?

readExpression <- function(filePath, header = FALSE){
  x <- read.table(filePath, header = header, row.names = 1)
  return(x)
}

readCoordinates <- function(filePath, header = FALSE){
  x <- read.table(filePath, header - header, row.names = 1)
  field <- x[, 1]
  x <- x[, -1]
  colnames(x) <- c("X", "Y")
  return(list(coords = x, field = field))
}

readGenes <- function(filePath, header = FALSE){
  x <- read.table(filePath, header = header)
  x <- as.character(x[[1]])
  return(x)
}

