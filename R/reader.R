#' @title Reads expression or main matrix from a text file
#' @description The function reads in a expression matrix and
#' assumes that the text file has no header
#' and each column is a sample/spot/cell and each row
#' is a gene/feature. The intersection of row and column
#' is assumed to be a normalized value. The first column
#' in the input file is the gene name or index. Second
#' column onwards are the expression data.
#' @param filePath Specify the file path for the input
#' expression data.
#' @param header Specify if the input file has a header
#' row. Default \code{FALSE}.
#' @return A data.frame from the input file
#' @export
readExpression <- function(filePath, header = FALSE){
  x <- read.table(filePath, header = header, row.names = 1)
  return(x)
}

#' @title Reads spatial coordinates from a text file
#' @description This functions reads in spatial coordinates
#' and fov from a text file. It assumes that text file has no
#' header and that there are three columns and as many
#' rows as there are number of samples/spots/cells. The
#' input text file has four columns. First column are just
#' the indices for the rows which will be ignored while
#' reading in the text file. Second column is the
#' FOV (field of vision) value for a particular sample.
#' The third and fourth columns are the X & Y coordinates for
#' the respective sample. The output is a list of two
#' elements. First element accessed by $coords are
#' the coordinates and the second element is the fov
#' accessed by $field.
#' @param filePath Specify the file path for the input
#' coordinates file.
#' @param header Specify if the input file has a header.
#' Default \code{FALSE}.
#' @return A list of 2 elements where the first one
#' is data.frame of X & Y coordinates and the second
#' element of the list is a numeric vector of FOV.
#' @export
readCoordinates <- function(filePath, header = FALSE){
  x <- read.table(filePath, header - header, row.names = 1)
  field <- x[, 1]
  x <- x[, -1]
  colnames(x) <- c("sdimx", "sdimy")
  return(list(coords = x, field = field))
}

#' @title Reads names of genes/features from a text file
#' @description This function reads in genes or feature names
#' from a text file and outputs a vector of these names.
#' The function assumes that each gene name is one
#' line in the text file.
#' @param filePath Specify the file path of the file
#' containing the gene or feature names.
#' @param header Specify if the input file has a header.
#' Default \code{FALSE}.
#' @return A character vector for gene or feature names.
#' @export
readGenes <- function(filePath, header = FALSE){
  x <- read.table(filePath, header = header)
  x <- as.character(x[[1]])
  return(x)
}
