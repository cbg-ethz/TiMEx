% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentDatasets.R
\name{gbmMuexSubsampling}
\alias{gbmMuexSubsampling}
\title{Stability of mutually exclusive groups in the glioblastoma datased used 
by muex}
\format{\code{gbmMuexSubsampling} is a list with as many elements as 
subsampling frequencies provided. Each element is further a list with as 
many elements as number of sizes of the significantly mutually exclusive 
groups identified. Additionally, \code{bonf} and \code{fdr} are two lists 
corresponding to each of these elements, representing different multiple 
correction methods. Finally, each element is a vector of subsampling 
frequencies of the significant mutually exclusive groups identified. For 
example, \code{gbmMuexSubsampling[[1]][[3]]} represents the relative counts 
of the identified mutually exclusive groups of size 3 for a subsampling 
frequency of 30\%, for both  \code{fdr} and \code{bonf} (bonferroni) 
multiple correction methods.}
\source{
Produced with the function \code{\link{subsampleAnalysis}}, ran with
the inputs 

\code{subsampl<-c(0.3,0.5,0.8)}

\code{noReps<-100}

and the mutually exclusive groups from \code{\link{gbmMuexOutput}}.
}
\description{
Dataset containing the stability of the mutually exclusive groups identified
by TiMEx in the glioblastoma dataset used in Szczurek \emph{et. al}: 
\emph{Modeling mutual exclusivity of 
cancer mutations}. Research in Computational Molecular Biology (2014), 
after subsampling the set of patients at frequencies of \code{30\%}, 
\code{50\%}, and \code{80\%}, 100 times.
}
