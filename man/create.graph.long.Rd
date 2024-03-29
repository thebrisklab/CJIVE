% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{create.graph.long}
\alias{create.graph.long}
\title{Function for plotting networks with ggplot}
\usage{
create.graph.long(gmatrix, sort_indices = NULL)
}
\arguments{
\item{gmatrix}{square matrix of size p-by-p in which entries represent the strength of (un-directed) edges between the p nodes}

\item{sort_indices}{vector of length p by which nodes are sorted. If NULL, then nodes are not sorted. Default is NULL.}
}
\value{
a data frame of three variables: X1, which represents the row from which the edge comes;  X2, which represents the column from which the edge comes;
       3) value, matrix entry representing the strength of the edge between the nodes represented by X1 and X2
}
\description{
Convert matrix representation of a network for graphical display via ggplot
}
