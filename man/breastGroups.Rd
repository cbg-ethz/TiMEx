% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentDatasets.R
\name{breastGroups}
\alias{breastGroups}
\title{Metagroups of genes in breast cancer}
\format{\code{breastGroups} is a list with 273 elements, where each 
element is a vector of genes with identical alteration patterns as the 
current gene. The numbers indicate the positions of the genes in the input
matrix.}
\source{
Produced with the function \code{\link{doMetagene}}.
}
\description{
Dataset containing the genes with identical alteration patterns in the 
breast cancer dataset \code{\link{breast}} (before preprocessing). It 
is represented as a list of metagenes, with as many elements as input 
genes that had an identical alteration pattern with at least one other 
input gene.
}
