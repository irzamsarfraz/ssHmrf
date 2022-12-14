% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reader.R
\name{readExpression}
\alias{readExpression}
\title{Reads expression or main matrix from a text file}
\usage{
readExpression(filePath, header = FALSE)
}
\arguments{
\item{filePath}{Specify the file path for the input
expression data.}

\item{header}{Specify if the input file has a header
row. Default \code{FALSE}.}
}
\value{
A data.frame from the input file
}
\description{
The function reads in a expression matrix and
assumes that the text file has no header
and each column is a sample/spot/cell and each row
is a gene/feature. The intersection of row and column
is assumed to be a normalized value. The first column
in the input file is the gene name or index. Second
column onwards are the expression data.
}
