#' TiMEx: A package for finding mutually exclusive groups of alterations in 
#' large cancer datasets
#'
#' @section Overview:
#'The most important function in this package is 
#'\code{\link{TiMEx}}, which identifies all mutually
#'exclusive groups in a binary dataset. \code{\link{TiMEx}} is a procedure 
#'implementing three steps: first, all pairs in the input dataset are tested 
#'for mutual exclusivity. Second, maximal cliques are identified on the basis of a 
#'selected number of pairs. Third, the resulting cliques are tested for mutual
#'exclusivity. The options of \code{\link{TiMEx}} 
#'are thresholds on the significance and intensity of mutual exclusivity 
#'of mutually exclusive pairs (\code{pairMu} and 
#'\code{pairPvalue}) and q-value cutoff on the identification of groups 
#'(\code{groupPvalue}). Unless
#'otherwise specified, \code{\link{TiMEx}} will use the default values of 
#'these options. Alternatively, the users 
#'interested in running separately the three steps of the TiMEx procedure 
#'should run, in this order, the functions \code{\link{analyzePairs}},
#'\code{\link{doMaxCliques}}, and \code{\link{findSignifCliques}}.
#'
#'@section Preprocessing and postprocessing:
#'Moreover, this package provides functions to pre-process the input data 
#'(\code{\link{doMetagene}},\code{\link{removeLowFreqs}}), to post-process the
#'resulting groups (\code{\link{produceTablesSignifGroups}}, 
#'\code{\link{subsampleAnalysis}}, \code{plotGroupByName}, 
#'\code{\link{recoverAllNamesGroups}}), as well as to simulate a dataset 
#'generated from the TiMEx model (\code{\link{simulateGenes}}).
#'
#' @section Datasets:
#' Multiple datasets are available within this package. 
#' \code{\link{gbmDendrix}} is a glioblastoma dataset used by ... in ..., 
#' \code{\link{breast}} and \code{\link{ovarian}} are datasets 
#' downloaded from the cBio Portal and preprocessed as explained in "TiMEx: a
#' ..."
#'
#'@section More:
#'For more in-depth explanations of the TiMEx package, including working 
#'examples, please see the corresponding vignette
#'
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015)

#' @docType package
#' @name TiMEx-package
#' @aliases TiMEx-package
NULL