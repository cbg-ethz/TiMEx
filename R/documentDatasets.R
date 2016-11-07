############################################################################### 
#' Breast cancer dataset
#'
#' Dataset containing a binary alteration pattern for the breast cancer 
#' dataset downloaded from cBioPortal (TCGA) in July 2014, and preprocessed
#' as described in Constantinescu \emph{et al.}: \emph{TiMEx: A Waiting Time 
#' Model for Mutually Exclusive Cancer Alterations}. Bioinformatics (2015). 
#' Rows represent patients, and columns represent alterations. 
#'
#' @format \code{breast} is a binary matrix with 958 rows and 537 columns.
#' @source \url{http://www.cbioportal.org/study.do?cancer_study_id=brca_tcga}
#' @name breast
#' @aliases breast
NULL



############################################################################### 
#' Metagroups of genes in breast cancer
#'
#' Dataset containing the genes with identical alteration patterns in the 
#' breast cancer dataset \code{\link{breast}} (before preprocessing). It 
#' is represented as a list of metagenes, with as many elements as input 
#' genes that had an identical alteration pattern with at least one other 
#' input gene.
#'
#' @format \code{breastGroups} is a list with 273 elements, where each 
#' element is a vector of genes with identical alteration patterns as the 
#' current gene. The numbers indicate the positions of the genes in the input
#' matrix.
#' @source Produced with the function \code{\link{doMetagene}}.
#' @name breastGroups
#' @aliases breastGroups
NULL



############################################################################### 
#' Mutually exclusive groups in breast cancer
#'
#' Dataset containing the groups identified as significantly 
#' mutually exclusive by TiMEx in breast cancer, together with their 
#' intensities of mutual exclusivity, corrected p-values, and other 
#' information.
#'
#' @format \code{breastOutput} is a list consisting of:
#' \itemize{
#' \item{\code{genesSignif}} {list of significantly mutually exclusive groups, 
#' as gene names, sorted by corrected p-value. The list contains as many 
#' elements as identified lengths of groups. For example, 
#' \code{genesSignif[[2]]} is a list containing the gene names of the 
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent gene names of significantly mutually exclusive groups.}
#' 
#' \item{\code{idxSignif}} {list of significantly mutually exclusive groups, as
#' indices in the input matrix, sorted by corrected p-value. The list 
#' contains as many elements as identified lengths of groups. For example,
#' \code{idxSignif[[2]]} is a list containing the indices of the  
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent indices of significantly mutually exclusive groups.}
#' 
#' \item{\code{pvals}} {list of corrected significant p-values corresponding to
#' the tested cliques, ordered ascendingly. The list contains as many elements 
#' as identified lengths of significant groups. For example, \code{pvals[[2]]}
#' is a list containing the p-values of the significant maximal cliques of 
#' size 2. Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a vector, of length the number of significant 
#' maximal cliques of a given size.} 
#' 
#' \item{\code{posSignif}} {list of positions of the significant groups in the 
#' input list of maximal cliques, ordered ascendingly by corrected p-value. 
#' The list contains as many elements as identified lengths of significant 
#' groups. For example, \code{posSignif[[2]]} is a list containing the 
#' positions of the significant groups of size 2.  Each list of this type 
#' further has two elements, \code{fdr} and \code{bonf}, corresponding to 
#' different multiple correction methods.  Each element is a vector, of length 
#' the number of significant maximal cliques of a given size.}
#' 
#' \item{\code{MusGroup}} {list of inferred mu values corresponding to
#' the tested cliques, ordered ascendingly by the corresponding corrected 
#' p-value. The list contains as many elements as identified lengths of 
#' significant groups. For example, \code{MusGroup[[2]]} is a list containing 
#' the mu values of the significant maximal cliques of size 2. Each list of 
#' this type further has two elements, \code{fdr} and \code{bonf}, 
#' corresponding to different multiple testing correction methods. Each 
#' element is a vector, of length the number of significant maximal cliques 
#' of a given size.} 
#' 
#' \item{\code{mcStruct}} {input structure of maximal cliques to be tested 
#' for mutual exclusivity, as returned by \code{\link{doMaxCliques}}.}
#' 
#' \item{\code{matrix}} {input binary alteration matrix.}
#' 
#' \item{\code{groupPvalue}} {input threshold for the corrected p-value, lower
#' than which cliques are significant.}
#' }
#' 
#' @source Produced with the function \code{\link{TiMEx}}, on the binary matrix
#' in the input dataset \code{\link{breast}}.
#' @name breastOutput
#' @aliases breastOutput
NULL



############################################################################### 
#' Stability of mutually exclusive groups in breast cancer
#'
#' Dataset containing the stability of the mutually exclusive groups identified
#' by TiMEx in the breast cancer dataset \code{\link{breast}}, after 
#' subsampling the set of patients at frequencies of \code{30\%}, 
#' \code{50\%}, and \code{80\%}, for 100 times.
#'
#' @format  \code{breastSubsampling} is a list with as many elements as 
#' subsampling frequencies provided (3 in this case). Each element is further 
#' a list with as many elemenets as number of sizes of the significantly 
#' mutually exclusive groups identified. Aditionally, \code{bonf} and 
#' \code{fdr} are two lists corresponding to each of these elements, 
#' representing different multiple correction methods. Finally, each element 
#' is a vector of relative counts of the significantly mutually exclusive 
#' groups identified. For example, \code{breastSubsampling[[1]][[3]]} 
#' represents the relative counts of the identified mutually exclusive groups 
#' of size 3 for a subsampling frequency of 30\%, for both  \code{fdr} and 
#' \code{bonf} (bonferroni) multiple correction methods.
#' 
#' @source Produced with the function \code{\link{subsampleAnalysis}}, ran with
#' the inputs \code{subsampl<-c(0.3,0.5,0.8)}, \code{noReps<-100}, and the
#' mutually exclusive groups from \code{\link{breastOutput}}.
#' @name breastSubsampling
#' @aliases breastSubsampling
NULL



############################################################################### 
#' Breast cancer subtypes
#'
#' Dataset containing binary alteration patterns for the breast cancer subtypes
#' 'luminalA', 'luminalB', 'Her2', and 'Basal2', downloaded from cBioPortal 
#' (TCGA) in July 2014, and preprocessed as explained in 'TiMEx: A Waiting 
#' Time Model For Mutually Exclusive Cancer Alterations', by Constantinescu 
#' \emph{et al.}  (Bioinformatics, 2016). Rows represent patients, and columns 
#' represent alterations.
#'
#' @format \code{breastSubtypes} is a list with 4 elements, corresponding to 
#' the 4 breast subtypes. Each element is binary matrix with 537 columns, as 
#' follows: 
#' \code{breastSubtypes$luminalA} consists of 222 rows,
#' \code{breastSubtypes$luminalB} consists of 125 rows,
#' \code{Her2} consists of 55 rows, and \code{Basal} consists of 76 rows. 
#' @source \url{http://www.cbioportal.org/study.do?cancer_study_id=brca_tcga}
#' @name breastSubtypes
#' @aliases breastSubtypes
NULL



############################################################################### 
#' Mutually exclusive groups in breast cancer subtypes
#'
#' Dataset containing the groups identified as significantly 
#' mutually exclusive by TiMEx in each of the 4 breast cancer subtypes 
#' 'luminalA', 'luminalB', 'Her2', and 'Basal', together with their 
#' corresponding intensities of mutual exclusivity, corrected p-values, and 
#' other information.
#'
#' @format \code{breastSubtypesOutput} is a list consisting of 4 lists, 
#' corresponding 
#' to the 4 breast subtypes. Each list further consists of:
#' \itemize{
#' \item{\code{genesSignif}} {list of significantly mutually exclusive groups, 
#' as gene names, sorted by corrected p-value. The list contains as many 
#' elements as identified lengths of groups. For example, 
#' \code{genesSignif[[2]]} is a list containing the gene names of the 
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent gene names of significantly mutually exclusive groups.}
#' 
#' \item{\code{idxSignif}} {list of significantly mutually exclusive groups, as
#' indices in the input matrix, sorted by corrected p-value. The list 
#' contains as many elements as identified lengths of groups. For example,
#' \code{idxSignif[[2]]} is a list containing the indices of the  
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent indices of significantly mutually exclusive groups.}
#' 
#' \item{\code{pvals}} {list of corrected significant p-values corresponding to
#' the tested cliques, ordered ascendingly. The list contains as many elements 
#' as identified lengths of significant groups. For example, \code{pvals[[2]]}
#' is a list containing the p-values of the significant maximal cliques of 
#' size 2. Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a vector, of length the number of significant 
#' maximal cliques of a given size.} 
#' 
#' \item{\code{posSignif}} {list of positions of the significant groups in the 
#' input list of maximal cliques, ordered ascendingly by corrected p-value. 
#' The list contains as many elements as identified lengths of significant 
#' groups. For example, \code{posSignif[[2]]} is a list containing the 
#' positions of the significant groups of size 2.  Each list of this type 
#' further has two elements, \code{fdr} and \code{bonf}, corresponding to 
#' different multiple correction methods.  Each element is a vector, of length 
#' the number of significant maximal cliques of a given size.}
#' 
#' \item{\code{MusGroup}} {list of inferred mu values corresponding to
#' the tested cliques, ordered ascendingly by the corresponding corrected 
#' p-value. 
#' The list contains as many elements as identified lengths of significant 
#' groups. 
#' For example, \code{MusGroup[[2]]} is a list containing the mu values of the 
#' significant maximal cliques of size 2. Each list of this type further has 
#' two elements, \code{fdr} and \code{bonf}, corresponding to different 
#' multiple testing correction methods. Each element is a vector, of length 
#' the number of significant maximal cliques of a given size.} 
#' 
#' \item{\code{mcStruct}} {input structure of maximal cliques to be tested 
#' for mutual exclusivity, as returned by \code{\link{doMaxCliques}}.}
#' 
#' \item{\code{matrix}} {input binary alteration matrix.}
#' 
#' \item{\code{groupPvalue}} {input threshold for the corrected p-value, lower
#' than which cliques are significant.}
#' }
#' 
#' @source Produced with the function \code{\link{TiMEx}}, on the four binary
#' matrices in the input dataset \code{\link{breastSubtypes}}.
#' @name breastSubtypesOutput
#' @aliases breastSubtypesOutput
NULL



############################################################################### 
#' Glioblastoma dataset used by Dendrix
#'
#' Dataset containing a binary alteration pattern for the glioblastoma 
#' dataset used by Leiserson \emph{et. al} in 'Simultaneous identification of 
#' multiple driver pathways in cancer' (Plos Computational Biology, 2013). 
#' Rows represent patients, and columns represent alterations. In the names 
#' of alterations, \emph{(D)} represents a copy number deletion, and 
#' \emph{(A)} represents a copy number amplification. 
#'
#' @format \code{gbmDendrix} is a binary matrix with 261 rows and 486 columns.
#'
#' @source \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/
#' journal.pcbi.1003054}
#' @name gbmDendrix
#' @aliases gbmDendrix
NULL



############################################################################### 
#' Mutually exclusive groups in the glioblastoma dataset used by Dendrix
#'
#' Dataset containing the groups identified as significantly 
#' mutually exclusive by TiMEx in the glioblastoma dataset used by Leiserson 
#' \emph{et. al} in 'Simultaneous identification of multiple driver pathways 
#' in cancer' (Plos Computational Biology, 2013), together with their 
#' intensities of mutual exclusivity, corrected p-values, and other 
#' information.
#' 
#' @format \code{gbmDendrixOutput} is a list consisting of:
#' \itemize{
#' \item{\code{genesSignif}} {list of significantly mutually exclusive groups, 
#' as gene names, sorted by corrected p-value. The list contains as many 
#' elements as identified lengths of groups. For example, 
#' \code{genesSignif[[2]]} is a list containing the gene names of the 
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent gene names of significantly mutually exclusive groups.}
#' 
#' \item{\code{idxSignif}} {list of significantly mutually exclusive groups, as
#' indices in the input matrix, sorted by corrected p-value. The list 
#' contains as many elements as identified lengths of groups. For example,
#' \code{idxSignif[[2]]} is a list containing the indices of the  
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent indices of significantly mutually exclusive groups.}
#' 
#' \item{\code{pvals}} {list of corrected significant p-values corresponding to
#' the tested cliques, ordered ascendingly. The list contains as many elements 
#' as identified lengths of significant groups. For example, \code{pvals[[2]]}
#' is a list containing the p-values of the significant maximal cliques of 
#' size 2. Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a vector, of length the number of significant 
#' maximal cliques of a given size.} 
#' 
#' \item{\code{posSignif}} {list of positions of the significant groups in the 
#' input list of maximal cliques, ordered ascendingly by corrected p-value. 
#' The 
#' list contains as many elements as identified lengths of significant groups.
#' For example, \code{posSignif[[2]]} is a list containing the positions of 
#' the 
#' significant groups of size 2.  Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple
#' correction methods.  Each element is a vector, of length the number of
#' significant maximal cliques of a given size.}
#' 
#' \item{\code{MusGroup}} {list of inferred mu values corresponding to
#' the tested cliques, ordered ascendingly by the corresponding corrected 
#' p-value. 
#' The list contains as many elements as identified lengths of significant 
#' groups. 
#' For example, \code{MusGroup[[2]]} is a list containing the mu values of the 
#' significant maximal cliques of size 2. Each list of this type further has 
#' two elements, \code{fdr} and \code{bonf}, corresponding to different 
#' multiple testing correction methods. Each element is a vector, of length 
#' the number of significant maximal cliques of a given size.} 
#' 
#' \item{\code{mcStruct}} {input structure of maximal cliques to be tested 
#' for mutual exclusivity, as returned by \code{\link{doMaxCliques}}.}
#' 
#' \item{\code{matrix}} {input binary alteration matrix.}
#' 
#' \item{\code{groupPvalue}} {input threshold for the corrected p-value, lower
#' than which cliques are significant.}
#' }
#' 
#' @source Produced with the function \code{\link{TiMEx}}, on the binary matrix
#' in the input dataset \code{\link{gbmDendrix}}.
#' @name gbmDendrixOutput
#' @aliases gbmDendrixOutput
NULL



############################################################################### 
#' Stability of mutually exclusive groups in the glioblastoma datased used 
#' by Dendrix
#'
#' 
#' Dataset containing the stability of the mutually exclusive groups identified
#' by TiMEx in the glioblastoma dataset \code{\link{gbmDendrix}}, used by
#' Leiserson \emph{et. al} in 'Simultaneous identification of multiple driver 
#' pathways in cancer' (Plos Computational Biology, 2013), after subsampling 
#' the set of patients at frequencies of \code{30\%}, \code{50\%}, and 
#' \code{80\%}, for 100 times. 
#'
#' @format \code{gbmDendrixSubsampling} is a list with as many elements as 
#' subsampling frequencies provided (3 in this case). Each element is further 
#' a list with as many elemenets as number of sizes of the significantly 
#' mutually exclusive groups identified. Aditionally, \code{bonf} and 
#' \code{fdr} are two lists corresponding to each of these elements, 
#' representing different multiple correction methods. Finally, each element 
#' is a vector of relative counts of the significantly mutually exclusive 
#' groups identified. For example, \code{gbmDendrixSubsampling[[1]][[3]]} 
#' represents the relative counts of the identified mutually exclusive groups 
#' of size 3 for a subsampling frequency of 30\%, for both  \code{fdr} and 
#' \code{bonf} (bonferroni) multiple correction methods.
#'  
#' @source Produced with the function \code{\link{subsampleAnalysis}}, ran with
#' the inputs \code{subsampl<-c(0.3,0.5,0.8)}, \code{noReps<-100}, and the
#' mutually exclusive groups from \code{\link{gbmDendrixOutput}}.
#' 
#' @name gbmDendrixSubsampling
#' @aliases gbmDendrixSubsampling
NULL



############################################################################### 
#' Glioblastoma dataset used by muex
#'
#' Dataset containing a binary alteration pattern for the glioblastoma 
#' dataset used by Szczurek \emph{et. al} in 'Modeling mutual exclusivity of 
#' cancer mutations' (Research in Computational Molecular Biology, 2014). 
#' Rows represent patients, and columns represent alterations. 
#'
#' @format \code{gbmMuex} is a binary matrix with 236 rows and 83 columns.
#'
#' @source \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/
#' journal.pcbi.1003503}
#' @name gbmMuex
#' @aliases gbmMuex
NULL



############################################################################### 
#' Mutually exclusive groups in the glioblastoma dataset used by muex
#'
#' Dataset containing the groups identified as significantly 
#' mutually exclusive by TiMEx in the glioblastoma dataset used by Szczurek
#' \emph{et. al} in 'Modeling mutual exclusivity of cancer mutations' 
#' (Research in Computational Molecular Biology, 2014), together with their 
#' intensities of mutual exclusivity, corrected p-values, and other 
#' information. The dataset which was used as input for producing these 
#' groups can be accessed via \code{data(gbm)} in the R package \code{muex},
#' available at \url{https://www1.ethz.ch/bsse/cbg/software/muex}.
#'
#' @format \code{gbmMuexOutput} is a list consisting of:
#' \itemize{
#' \item{\code{genesSignif}} {list of significantly mutually exclusive groups, 
#' as gene names, sorted by corrected p-value. The list contains as many 
#' elements as identified lengths of groups. For example, 
#' \code{genesSignif[[2]]} is a list containing the gene names of the 
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent gene names of significantly mutually exclusive groups.}
#' 
#' \item{\code{idxSignif}} {list of significantly mutually exclusive groups, as
#' indices in the input matrix, sorted by corrected p-value. The list 
#' contains as many elements as identified lengths of groups. For example,
#' \code{idxSignif[[2]]} is a list containing the indices of the  
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent indices of significantly mutually exclusive groups.}
#' 
#' \item{\code{pvals}} {list of corrected significant p-values corresponding to
#' the tested cliques, ordered ascendingly. The list contains as many elements 
#' as identified lengths of significant groups. For example, \code{pvals[[2]]}
#' is a list containing the p-values of the significant maximal cliques of 
#' size 2. Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a vector, of length the number of significant 
#' maximal cliques of a given size.} 
#' 
#' \item{\code{posSignif}} {list of positions of the significant groups in the 
#' input list of maximal cliques, ordered ascendingly by corrected p-value. 
#' The 
#' list contains as many elements as identified lengths of significant groups.
#' For example, \code{posSignif[[2]]} is a list containing the positions of 
#' the 
#' significant groups of size 2.  Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple
#' correction methods.  Each element is a vector, of length the number of
#' significant maximal cliques of a given size.}
#' 
#' \item{\code{MusGroup}} {list of inferred mu values corresponding to
#' the tested cliques, ordered ascendingly by the corresponding corrected 
#' p-value. 
#' The list contains as many elements as identified lengths of significant 
#' groups. 
#' For example, \code{MusGroup[[2]]} is a list containing the mu values of the 
#' significant maximal cliques of size 2. Each list of this type further has 
#' two elements, \code{fdr} and \code{bonf}, corresponding to different 
#' multiple testing correction methods. Each element is a vector, of length 
#' the number of significant maximal cliques of a given size.} 
#' 
#' \item{\code{mcStruct}} {input structure of maximal cliques to be tested 
#' for mutual exclusivity, as returned by \code{\link{doMaxCliques}}.}
#' 
#' \item{\code{matrix}} {input binary alteration matrix.}
#' 
#' \item{\code{groupPvalue}} {input threshold for the corrected p-value, lower
#' than which cliques are significant.}
#' }
#' 
#' @source Produced with the function \code{\link{TiMEx}}, on the binary matrix
#' which can be accessed via \code{data(gbm)} in the R package \code{muex}.
#' @name gbmMuexOutput
#' @aliases gbmMuexOutput
NULL



############################################################################### 
#' Stability of mutually exclusive groups in the glioblastoma datased used 
#' by muex
#'
#' Dataset containing the stability of the mutually exclusive groups identified
#' by TiMEx in the glioblastoma dataset used by Szczurek \emph{et. al} in 
#' 'Modeling mutual exclusivity of cancer mutations' (Research in 
#' Computational Molecular Biology, 2014), after subsampling the set
#' of patients at frequencies of \code{30\%}, \code{50\%}, and \code{80\%}, for
#' 100 times. 
#'
#' @format \code{gbmMuexSubsampling} is a list with as many elements as 
#' subsampling frequencies provided. Each element is further a list with as 
#' many elemenets as number of sizes of the significantly mutually exclusive 
#' groups identified. Aditionally, \code{bonf} and \code{fdr} are two lists 
#' corresponding to each of these elements, representing different multiple 
#' correction methods. Finally, each element is a vector of subsampling 
#' frequencies of the significant mutually exclusive groups identified. For 
#' example, \code{gbmMuexSubsampling[[1]][[3]]} represents the relative counts 
#' of the identified mutually exclusive groups of size 3 for a subsampling 
#' frequency of 30\%, for both  \code{fdr} and \code{bonf} (bonferroni) 
#' multiple correction methods.
#' 
#' @source Produced with the function \code{\link{subsampleAnalysis}}, ran with
#' the inputs \code{subsampl<-c(0.3,0.5,0.8)}, \code{noReps<-100}, and the
#' mutually exclusive groups from \code{\link{gbmMuexOutput}}.
#' @name gbmMuexSubsampling
#' @aliases gbmMuexSubsampling
NULL



############################################################################### 
#' Ovarian cancer dataset
#'
#' Dataset containing a binary alteration pattern for the ovarian cancer 
#' dataset downloaded from cBioPortal (TCGA) in July 2014, and preprocessed
#' as explained in 'TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations', by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2016). Rows represent patients, and columns represent 
#' alterations.  
#'
#' @format A binary matrix with 316 rows and 312 columns.
#' @source \url{http://www.cbioportal.org/study.do?cancer_study_id=ov_tcga_pub}
#' @name ovarian
#' @aliases ovarian
NULL



############################################################################### 
#' Metagroups of genes in ovarian cancer
#'
#' Dataset containing the genes with identical alteration patterns in the 
#' ovarian cancer dataset \code{\link{ovarian}} (before preprocessing). It 
#' is represented as a list of metagenes, with as many elements as input 
#' genes which had an identical alteration pattern with at least one other 
#' input gene.
#'
#' @format \code{ovarianGroups} is a list with 263 elements, where each 
#' element is a vector of genes with identical alteration patterns as the 
#' current gene. The numbers indicate the positions of the genes in the input
#' matrix.
#' 
#' @source Produced with the function \code{\link{doMetagene}}.
#' @name ovarianGroups
#' @aliases ovarianGroups
NULL



############################################################################### 
#' Mutually exclusive groups in ovarian cancer
#'
#' Dataset containing the groups identified as significantly 
#' mutually exclusive by TiMEx in ovarian cancer, together with their 
#' intensities of mutual exclusivity, corrected p-values, and other 
#' information.
#'
#' @format \code{ovarianOutput} is a list consisting of:
#' \itemize{
#' \item{\code{genesSignif}} {list of significantly mutually exclusive groups, 
#' as gene names, sorted by corrected p-value. The list contains as many 
#' elements as identified lengths of groups. For example, 
#' \code{genesSignif[[2]]} is a list containing the gene names of the 
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent gene names of significantly mutually exclusive groups.}
#' 
#' \item{\code{idxSignif}} {list of significantly mutually exclusive groups, as
#' indices in the input matrix, sorted by corrected p-value. The list 
#' contains as many elements as identified lengths of groups. For example,
#' \code{idxSignif[[2]]} is a list containing the indices of the  
#' significant groups of size 2. Each list of this type further has two 
#' elements, \code{fdr} and \code{bonf}, corresponding to different multiple 
#' testing correction methods. Each element is a matrix, in which rows 
#' represent indices of significantly mutually exclusive groups.}
#' 
#' \item{\code{pvals}} {list of corrected significant p-values corresponding to
#' the tested cliques, ordered ascendingly. The list contains as many elements 
#' as identified lengths of significant groups. For example, \code{pvals[[2]]}
#' is a list containing the p-values of the significant maximal cliques of 
#' size 2. Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a vector, of length the number of significant 
#' maximal cliques of a given size.} 
#' 
#' \item{\code{posSignif}} {list of positions of the significant groups in the 
#' input list of maximal cliques, ordered ascendingly by corrected p-value. 
#' The list contains as many elements as identified lengths of significant 
#' groups.For example, \code{posSignif[[2]]} is a list containing the 
#' positions of the significant groups of size 2.  Each list of this type 
#' further has two elements, \code{fdr} and \code{bonf}, corresponding to 
#' different multiple correction methods.  Each element is a vector, of 
#' length the number of significant maximal cliques of a given size.}
#' 
#' \item{\code{MusGroup}} {list of inferred mu values corresponding to
#' the tested cliques, ordered ascendingly by the corresponding corrected 
#' p-value. The list contains as many elements as identified lengths of 
#' significant groups. For example, \code{MusGroup[[2]]} is a list containing 
#' the mu values of the significant maximal cliques of size 2. Each list of 
#' this type further has two elements, \code{fdr} and \code{bonf}, 
#' corresponding to different multiple testing correction methods. Each 
#' element is a vector, of length the number of significant maximal cliques 
#' of a given size.} 
#' 
#' \item{\code{mcStruct}} {input structure of maximal cliques to be tested 
#' for mutual exclusivity, as returned by \code{\link{doMaxCliques}}.}
#' 
#' \item{\code{matrix}} {input binary alteration matrix.}
#' 
#' \item{\code{groupPvalue}} {input threshold for the corrected p-value, lower
#' than which cliques are significant.}
#' }
#' 
#' @source Produced with the function \code{\link{TiMEx}}, on the binary matrix
#' in the input dataset \code{\link{ovarian}}.
#' 
#' @name ovarianOutput
#' @aliases ovarianOutput
NULL



############################################################################### 
#' Stability of mutually exclusive groups in ovarian cancer
#'
#' Dataset containing the stability of the mutually exclusive groups identified
#' by TiMEx in the ovarian cancer dataset \code{\link{ovarian}}, after 
#' subsampling the set of patients at frequencies of \code{30\%}, 
#' \code{50\%}, and \code{80\%}, for 100 times.
#'
#' @format \code{ovarianSubsampling} is  a list with as many elements as 
#' subsampling frequencies provided (3 in this case). Each element is further 
#' a list with as many elemenets as number of sizes of the significantly 
#' mutually exclusive groups identified. Aditionally, \code{bonf} and 
#' \code{fdr} are two lists corresponding to each of these elements, 
#' representing different multiple correction methods. Finally, each element 
#' is a vector of relative counts of the significantly mutually exclusive 
#' groups identified. For example, \code{ovarianSubsampling[[1]][[3]]} 
#' represents the relative counts of the identified mutually exclusive groups 
#' of size 3 for a subsampling frequency of 30\%, for both  \code{fdr} and 
#' \code{bonf} (bonferroni) multiple correction methods.
#' 
#' @source Produced with the function \code{\link{subsampleAnalysis}}, ran with
#' the inputs \code{subsampl<-c(0.3,0.5,0.8)}, \code{noReps<-100}, and the
#' mutually exclusive groups from \code{\link{ovarianOutput}}.
#' 
#' @name ovarianSubsampling
#' @aliases ovarianSubsampling
NULL
