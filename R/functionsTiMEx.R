# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
#' @title Finds all mutually exclusive pairs in a dataset
#' 
#' @description \code{analyzePairs} performs step 1 of the TiMEx procedure:
#' tests all gene pairs for mutual exclusivity. It returns a complex list, 
#' including parameter estimates, pairwise p-values and intensities of mutual 
#' exclusivity, log likelihoods, and others.  
#' 
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' 
#' @details In the first step, the TiMEx procedure for identifying mutually
#' exclusive groups of alterations tests all pairs for mutual exclusivity. The 
#' data corresponding to each pair of genes is therefore fitted to both the 
#' Null (Conditional Independence) and the Mutual Exclusivity models. 
#' Parameters under the two models are estimated, and, since they are nested, 
#' a likelihood ratio test is performed between the corresponding log 
#' likelihoods, in order to test whether mu (the intensity of mutual 
#' exclusivity) is different from 0. For more details on the TiMEx procedure, 
#' as well as on the underlying mathematical model, see "TiMEx: A Waiting Time 
#' Model For Mutually Exclusive Cancer Alterations", by Constantinescu 
#' \emph{et al.} (Bioinformatics, 2015).
#' 
#' The output list contains exhaustive information on the pairwise testing of 
#' genes, including parameter estimates, pairwise p-values and intensities of 
#' mutual exclusivity, log likelihoods, and others. 
#' 
#' This function displays progress messages indicating the gene which is 
#' currently being tested against the remaining genes.
#' 
#' @return List consisting of a set of matrices, all of dimension \code{n x n} 
#' (\code{n} being the number of genes). Element \code{(i,j)} of each matrix 
#' corresponds to the pairwise test between genes \emph{i} and \emph{j}:
#' \itemize{
#' \item{\code{lamiEstNull}} {matrix with rate estimates (lambda) of the 
#' waiting time of one gene (gene \code{i}) under the Null model.}
#' \item{\code{lamjEstNull}} {matrix with rate estimates (lambda) of the 
#' waiting time of the other gene (gene \code{j}) under the Null model.}
#' \item{\code{lamiEstME}} {matrix with rate estimates (lambda) of the waiting
#' time of one gene (gene \code{i}) under the Mutual Exclusivity model.}
#' \item{\code{lamjEstME}} {matrix with rate estimates (lambda) of the 
#' waiting time of the other gene (gene \code{j}) under the Mutual Exclusivity 
#' model.}
#' \item{\code{likeNull}} {matrix with log likelihoods under the Null model.}
#' \item{\code{likeAlt}} {matrix with log likelihoods under the Mutual
#' Exclusivity model.}
#' \item{\code{LRT}} {matrix with log likelihood ratio test (LRT) statistics.}
#' \item{\code{pvalueLRTCorrectSym}} {list of the following symmetric matrices:
#' \itemize{
#' \item{\code{bonferroni}} {bonferroni corrected p-values corresponding to 
#' the LRT statistic.}
#' \item{\code{bonferroni}} {fdr corrected p-values corresponding to the 
#' LRT statistic.}
#' \item{\code{uncorrected}} {uncorrected  p-values corresponding to the 
#' LRT statistic.}
#' }
#' }
#' \item{\code{muEstSym}} {symmetric matrix with estimated pairwise mu 
#' (intensity of mutual exclusivity)}
#'}
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015)
#' 
#' @seealso \code{\link{doMaxCliques}} for step 2 of the TiMEx procedure;
#' \code{\link{findSignifCliques}} for step 3 of the TiMEx procedure;
#' the wrapper function \code{\link{TiMEx}} for combining these three steps, 
#' and identifying mutually exclusive groups in a  binary dataset with the 
#' TiMEx model.
#' 
#' @examples
#' # Test all pairs from the ovarian dataset for mutual exclusivity (takes 
#' # approximately 5 minutes)
#' data(ovarian)
#' \donttest{ovarianPairs<-analyzePairs(ovarian)}
#'
#' @export

analyzePairs<-function(mat)
{
    # check inputs
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("input needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    noGenes<-dim(mat)[2]
    
    if (is.null(colnames(mat)))
        colnames(mat)<-paste("gene",c(1:noGenes),sep="")
    
    muEst<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(muEst)<-colnames(mat)
    rownames(muEst)<-colnames(mat)
    
    lamiEstME<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(lamiEstME)<-colnames(mat)
    rownames(lamiEstME)<-colnames(mat)
    
    lamjEstME<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(lamjEstME)<-colnames(mat)
    rownames(lamjEstME)<-colnames(mat)
    
    lamiEstNull<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(lamiEstNull)<-colnames(mat)
    rownames(lamiEstNull)<-colnames(mat)
    
    lamjEstNull<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(lamjEstNull)<-colnames(mat)
    rownames(lamjEstNull)<-colnames(mat)
    
    likeNull<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(likeNull)<-colnames(mat)
    rownames(likeNull)<-colnames(mat)
    
    likeAlt<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(likeAlt)<-colnames(mat)
    rownames(likeAlt)<-colnames(mat)
    
    LRT<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(LRT)<-colnames(mat)
    rownames(LRT)<-colnames(mat)
    
    pvalueLRT<-matrix(NA,nrow=noGenes,ncol=noGenes)
    colnames(pvalueLRT)<-colnames(mat)
    rownames(pvalueLRT)<-colnames(mat)
    
    for (gene1 in 1:(noGenes-1))
    {
        print(paste("gene", gene1))
        for (gene2 in (gene1+1):noGenes)
        {
            genes<-mat[,c(gene1,gene2)] #the two genes to be tested
            ns<-computeContTable(genes) #their contingency table
            nsVect<-c(ns$n00,ns$n01,ns$n10,ns$n11)
            
        # estimate parameters under ME
        #try(opMu<-optim(par=c(0.01,0.01,0.01),fn=computeParmsMu,gr=NULL,
        # n=nsVect,method="L-BFGS-B",lower=c(1e-15,1e-15,0),
        # upper=c(Inf,Inf,1-(1e-15)),control=list(fnscale=-20,maxit=10000)))
        #while (!exists("opMu"))
        #{
        #  try(opMu<-optim(par=c(0.01+abs(rnorm(1)),0.01+abs(rnorm(1)),
        #  min(0.0001+abs(rnorm(1)),0.999)),fn=computeParmsMu,gr=NULL,
        #  n=nsVect,method="L-BFGS-B",lower=c(1e-15,1e-15,0),
        #  upper=c(Inf,Inf,1-(1e-15)),control=list(fnscale=-20,maxit=10000),
        #  silent=TRUE))
        #}
            opMu<-optim(par=c(0.01,0.01,0.01),fn=computeParmsMu,gr=NULL,
                        n=nsVect,
                        method="L-BFGS-B",lower=c(1e-15,1e-15,0),
                        upper=c(Inf,Inf,1-(1e-15)),
                        control=list(fnscale=-20,maxit=10000))
            
            
        # estimate parameters under Null
        #try(opNull<-optim(par=c(0.01,0.01),fn=computeParmsNull,gr=NULL,
        #n=nsVect,method="L-BFGS-B",lower=c(1e-15,1e-15),
        #upper=c(Inf,Inf),control=list(fnscale=-20,maxit=10000)))
        #while (!exists("opNull"))
        #{
        #try(opNull<-optim(par=c(0.01+abs(rnorm(1)),0.01+
        #abs(rnorm(1))),fn=computeParmsNull,gr=NULL,n=nsVect,
        #method="L-BFGS-B",lower=c(1e-15,1e-15),
        #upper=c(Inf,Inf),control=list(fnscale=-20,maxit=10000),
        #silent=TRUE))
        #}
            opNull<-optim(par=c(0.01,0.01),fn=computeParmsNull,gr=NULL,
                        n=nsVect,method="L-BFGS-B",lower=c(1e-15,1e-15),
                        upper=c(Inf,Inf),
                        control=list(fnscale=-20,maxit=10000))
            
            lamiEstME[gene1,gene2]<-opMu$par[1]
            lamiEstME[gene2,gene1]<-opMu$par[2]
            
            lamjEstME[gene1,gene2]<-opMu$par[2]
            lamjEstME[gene2,gene1]<-opMu$par[1]
            
            muEst[gene1,gene2]<-opMu$par[3]
            
            lamiEstNull[gene1,gene2]<-opNull$par[1]
            lamiEstNull[gene2,gene1]<-opNull$par[2]
            
            lamjEstNull[gene1,gene2]<-opNull$par[2]
            lamjEstNull[gene2,gene1]<-opNull$par[1]
            
            # compute likelihoods
            likeNull[gene1,gene2]<-computeParmsNull(c(lamiEstNull[gene1,gene2],
                                                    lamjEstNull[gene1,gene2]),
                                                    nsVect)
            likeAlt[gene1,gene2]<-computeParmsMu(c(lamiEstME[gene1,gene2],
                                                    lamjEstME[gene1,gene2], 
                                                    muEst[gene1,gene2]),nsVect)
            
            # perform the LRT
            LRT[gene1,gene2]<-2*likeAlt[gene1,gene2]-2*likeNull[gene1,gene2]
            pvalueLRT[gene1,gene2]<-pchisq(LRT[gene1,gene2], df = 1, ncp = 0, 
                                        lower.tail = FALSE)
        }
    }
    
    # correct the pvalues
    pvalueLRTCorrect<-list()
    pvalueLRTCorrect$bonferroni<-matrix(p.adjust(as.vector(pvalueLRT),
                                method="bonferroni"),nrow=noGenes,ncol=noGenes)
    colnames(pvalueLRTCorrect$bonferroni)<-colnames(pvalueLRT)
    rownames(pvalueLRTCorrect$bonferroni)<-rownames(pvalueLRT)
    pvalueLRTCorrect$fdr<-matrix(p.adjust(as.vector(pvalueLRT),method="fdr"),
                                nrow=noGenes,ncol=noGenes)
    colnames(pvalueLRTCorrect$fdr)<-colnames(pvalueLRT)
    rownames(pvalueLRTCorrect$fdr)<-rownames(pvalueLRT)
    pvalueLRTCorrect$uncorrected<-pvalueLRT
    
    # transform to symmetric matrices
    muEstSym<-makeSym(muEst)
    likeNull<-makeSym(likeNull)
    likeAlt<-makeSym(likeAlt)
    LRT<-makeSym(LRT)
    pvalueLRTSym<-makeSym(pvalueLRT)
    pvalueLRTCorrectSym<-list()
    pvalueLRTCorrectSym$bonferroni<-makeSym(pvalueLRTCorrect$bonferroni)
    pvalueLRTCorrectSym$fdr<-makeSym(pvalueLRTCorrect$fdr)
    pvalueLRTCorrectSym$uncorrected<-makeSym(pvalueLRTCorrect$uncorrected)
    
    results<-list("lamiEstNull"=lamiEstNull,"lamjEstNull"=lamjEstNull,
                "lamiEstME"=lamiEstME,"lamjEstME"=lamjEstME,
                "likeNull"=likeNull,"likeAlt"=likeAlt,"LRT"=LRT,
                "pvalueLRTCorrectSym"=pvalueLRTCorrectSym,"muEstSym"=muEstSym)
    
    return(results) 
}



###############################################################################
#' @title Identifies maximal cliques from pairwise testing information
#'
#' @description \code{doMaxCliques} performs step 2 of the TiMEx procedure: 
#' identifies maximal cliques using information from pairwise testing. The 
#' maximal clique detection routine only uses the connections between gene 
#' pairs which satisfy the thresholds on mu (\code{pairMu}) and pvalue 
#' (\code{pairPvalue}).
#' 
#' @param pairs list resulting after pairwise testing, as returned by
#' \code{\link{analyzePairs}}
#' @param pairMu pair-level threshold on mu (real number between 0 and 1). 
#' Default is 0.5.
#' @param pairPvalue pair-level threshold on p-value (real number between 0 and
#' 1). Default is 0.01.
#' 
#' @details In the second step, the TiMEx procedure for identifying mutually 
#' exclusive groups of alterations detects maximal cliques using pairwise 
#' testing information from step 1. A graph is constructed, in which genes
#' are vertices, and an edge is drawn between any pair \emph{(i,j)} if both the
#' estimated intensity of mutual exclusivity and the computed p-value satisfy 
#' the chosen thresholds \code{pairMu} and \code{pairPvalue}. Maximal
#' cliques are detected on this graph. 
#' 
#' The two thresholds can be set by the 
#' user, and are recommended to be chosen based on the sensitivity and 
#' specificity levels to which they correspond, as assessed in simulated data. 
#' For details, see "TiMEx: A Waiting Time Model For Mutually Exclusive Cancer 
#' Alterations", by Constantinescu \emph{et al.} (2015). The default values are
#' 0.5 for \code{pairMu} and 0.01 for \code{pairPvalue}.
#' 
#' This function needs functions from the packages \emph{RBGL} and 
#' \emph{igraph} to run.
#' 
#' @return list consisting of:
#' \itemize{
#' \item{\code{detectedLengths}} {vector of lengths of the identified maximal 
#' cliques.}
#' \item{\code{idxInCliques}} {list with as many elements as lengths of 
#' the identified maximal cliques. Each element of the list is a matrix, in 
#' which each row represents the indices of genes in an identified maximal 
#' clique.}
#' \item{\code{genesInCliques}} {list with as many elements as lengths of 
#' the identified maximal cliques. Each element of the list is a matrix, in 
#' which each row represents the names of genes in an identified maximal 
#' clique.}
#' \item{\code{noMaxCliques}} {vector of numbers of identified maximal cliques 
#' corresponding to each length present in the field \code{detectedLengths}.}
#' \item{\code{Mus}} {list of two elements: \code{OrderedGenesInCliques} and 
#' \code{OrderedIdxInCliques}, which have the same structure as the elements
#' \code{genesInCliques} and \code{idxInCliques}. The only difference is that
#' the identified maximal cliques are now ordered by their average pairwise
#' mu.}
#' \item \code{pairMu} input pair-level threshold on mu.
#' \item \code{pairPvalue} input pair-level threshold on p-value.
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso \code{\link{analyzePairs}} for step 1 of the TiMEx procedure;
#' \code{\link{findSignifCliques}} for step 3 of the TiMEx procedure;
#' the wrapper function \code{\link{TiMEx}} for combining these three steps, 
#' and identifying mutually exclusive groups in a  binary dataset with the 
#' TiMEx model.
#' 
#' @examples
#' # First, test all pairs from the ovarian cancer dataset for mutual 
#' # exclusivity (take approximately 5 minutes) 
#' data(ovarian)
#' \donttest{ovarianPairs<-analyzePairs(ovarian)}
#' 
#' # Then, identify all maximal cliques using the default thresholds
#' \donttest{ovarianMaxCliques<-doMaxCliques(ovarianPairs)}
#' 
#' @importFrom igraph graph.adjacency
#' @importFrom igraph igraph.to.graphNEL
#' @import RBGL
#' 
#' @export
doMaxCliques<-function(pairs,pairMu,pairPvalue)
{
    
    # check inputs
    if (missing(pairs))
        stop("a list containing pairwise testing information needs to be 
            specified")
    if (!class(pairs)=="list")
        stop("a list containing pairwise testing information needs to be 
            specified")
    if (!length(setdiff(c("muEstSym","pvalueLRTCorrectSym"),names(pairs)))==0)
        stop("the input list 'pairs' should have contained a field with 
            pairwise information on mu and a field with pairwise information 
            on p-value")
    
    if (class(pairs$pvalueLRTCorrectSym)!="list")
        stop("the field 'pvalueLRTCorrectSym' of the input list should have 
            an element named 'uncorrected'")
    if (length(setdiff(c("uncorrected"),names(pairs$pvalueLRTCorrectSym))==0))
        stop("the field 'pvalueLRTCorrectSym' of the input list should have 
            an element named 'uncorrected'")
    if(!length(unique(c(dim(pairs$muEstSym),
                        dim(pairs$pvalueLRTCorrectSym$uncorrected))))==1)
        stop("the matrices with pvalue and mu information should have the 
            same dimension and be squared")
    if (is.null(colnames(pairs$pvalueLRTCorrectSym$uncorrected)))
        colnames(pairs$pvalueLRTCorrectSym$uncorrected)<-
        paste("gene",c(1:dim(pairs$pvalueLRTCorrectSym$uncorrected)[1]),sep="")
    if (is.null(rownames(pairs$pvalueLRTCorrectSym$uncorrected)))
        rownames(pairs$pvalueLRTCorrectSym$uncorrected)<-
        paste("gene",c(1:dim(pairs$pvalueLRTCorrectSym$uncorrected)[1]),sep="")
    if (is.null(colnames(pairs$muEstSym)))
        colnames(pairs$muEstSym)<-
        paste("gene",c(1:dim(pairs$muEstSym)[1]),sep="")
    if (is.null(rownames(pairs$muEstSym)))
        rownames(pairs$muEstSym)<-
        paste("gene",c(1:dim(pairs$muEstSym)[1]),sep="")
    
    
    if (missing(pairMu))
        pairMu<-0.5
    if (!(class(pairMu)=="numeric"))
        stop("the threshold on mu needs to be a real number between 0 and 1")
    if (length(pairMu)>1)
        stop("the threshold on mu needs to be a real number between 0 and 1")
    if (!(0<=pairMu && pairMu<=1))
        stop("the threshold on mu needs to be a real number between 0 and 1")
    
    if (missing(pairPvalue))
        pairPvalue<-0.01
    if (!(class(pairPvalue)=="numeric"))
        stop("the threshold on p-value needs to be a real number between 0 
            and 1")
    if (length(pairPvalue)>1)
        stop("the threshold on p-value needs to be a real number between 0 
            and 1")
    if (!(0<=pairPvalue && pairPvalue<=1))
        stop("the threshold on p-value needs to be a real number between 0 
            and 1")
    
    cliques<-doClique(pairs$muEstSym,pairs$pvalueLRTCorrectSym$uncorrected,
                    pairMu,pairPvalue)
    mcStruct<-analyzeMaxCliquesBySize(cliques)
    mcStruct$Mus<-getMusAllMaxCliques(mcStruct,pairs$muEstSym)
    # do not return these fields anymore (redundant or unnecessary)
    mcStruct$Mus$detectedLengths<-NULL
    mcStruct$Mus$idxInCliques<-NULL
    mcStruct$Mus$genesInCliques<-NULL
    mcStruct$Mus$noMaxCliques<-NULL
    mcStruct$pairMu<-pairMu
    mcStruct$pairPvalue<-pairPvalue
    return(mcStruct)
}



###############################################################################
#' @title Tests whether a given group is mutually exclusive
#'
#' @description \code{testCliqueAsGroup} tests whether a group, given as gene
#' indices, is mutually exclusive.
#' 
#' @param geneIdx vector of indices in the input matrix of the genes to be 
#' tested
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' @param lo rate of observation time. Default is 1.
#' 
#' @details For deciding whether a group is mutually exclusive, the group of  
#' genes is fitted to both the Null (Conditional Independence) and the Mutual 
#' Exclusivity models. Parameters under the two models are estimated, and, 
#' since they are nested, a likelihood ratio test is performed between the 
#' corresponding log likelihoods, in order to test whether mu (the intensity of
#' mutual exclusivity) is different from 0. For computing the likelihood of 
#' the data under both models, an exhaustive enumeration of all possible 
#' orders of the input alterations needs to be performed. Therefore, the 
#' complexity of the test is exponential in the number of genes to be tested, 
#' which makes it unfeasible for large number of genes (usually more than 10).
#' 
#' \code{lo} (the rate of observation time) is by default set to 1, as both 
#' models are otherwise unindentifiable. We recommend leaving the value of this
#' parameter unchanged, as otherwise the estimated waiting time rates of the
#' genes require additional interpretation.
#' 
#' For more details on the TiMEx procedure, as well as on the underlying 
#' mathematical model, see "TiMEx: A Waiting Time Model For Mutually Exclusive 
#' Cancer Alterations", by Constantinescu \emph{et al.} (Bioinformatics, 2015).
#' 
#' @return List consisting of:
#' \itemize{
#' \item{\code{opMu}} {list as returned by \code{\link[stats]{optim}}. The 
#' field \code{par} is a vector of \code{n+1} positions (\code{n} being the 
#' number of genes) containing the estimates of the waiting time rates 
#' (lambda) for the \code{n} genes under the mutual exclusivity model, 
#' followed by the estimate for mu.}
#' \item{\code{opNull}} {list as returned by \code{\link[stats]{optim}}. The 
#' field \code{par} is a vector of \code{n} positions (\code{n} being the 
#' number of genes) containing the estimates of the waiting time rates 
#' (lambda) for the \code{n} genes under the null model.}
#' \item{\code{countsVec}} {the contingency table of the input genes, as a 
#' vector. The first element is the count of the samples where no gene was 
#' altered, the next \code{n} elements are the counts of the samples where 
#' exactly one gene was altered (starting with gene 1 and ending with gene 
#' \code{n}), the next \code{n} elements are the counts of the samples where
#' exactly two genes were altered (starting with genes 1 and 2, continuing with
#' genes 1 and 3, and ending with genes \code{n-1} and \code{n}), and so on.}
#' \item{\code{genes}} {subset of the input binary matrix, corresponding to the
#' genes to be tested}
#' \item{\code{LRT}} {log likelihood ratio (LRT)}
#' \item{\code{pvalueLRT}} {the p-value corresponding to the LRT}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso the wrapper function \code{\link{TiMEx}} for identifying
#' mutually exclusive groups in a  binary dataset with the TiMEx model.
#' 
#' @examples
#' # Tests for mutual exclusivity the group of genes with indices 13, 204, and 
#' # 310 in the ovarian cancer dataset
#' data(ovarian)
#' testGroup<-testCliqueAsGroup(c(13,204,310),ovarian)
#' 
#' @export
testCliqueAsGroup<-function(geneIdx,mat,lo)
{
    
    # check inputs
    if (missing(geneIdx))
        stop("need to specify a vector of gene indices to be tested for 
            mutual exclusivity")
    if (!(class(geneIdx)%in%c("numeric","integer")))
        stop("the gene indices need to be a numeric vector")
    if (!sum(geneIdx%%1)==0)
        stop("the gene indices need to be positive integers")
    if (any(geneIdx<=0))
        stop("the gene indices need to be positive integers")
    if (length(unique(geneIdx))!=length(geneIdx))
        stop("the gene indices need to be different")
    
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("the input 'mat' needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("the input 'mat' needs to be a binary matrix")
    
    if (any(geneIdx>dim(mat)[2]))
        stop("the gene indices need to be lower than the maximum number of 
            input genes")
    
    if (missing(lo))
        lo<-1
    if (lo!=1)
        warning("you are changing the default value of lo")
    
    
    genes<-mat[,geneIdx]
    n<-length(geneIdx)
    l<-list()
    for (i in 1:n)
    {
        l[[i]]<-seq(0,1)
    }
    allGenos<-as.matrix(expand.grid(l))
    
    tabNonZero<-table(unlist(sapply(c(1:dim(genes)[1]), 
        function(y,genes,allGenos) {which(apply(allGenos,1,function(x,y)
        {w<-genes[y,]; names(w)<-NULL; all.equal(as.vector(x),w)},y=y)==
        "TRUE")},genes=genes, allGenos=allGenos)))
    
    countsVec<-rep(0,dim(allGenos)[1])
    countsVec[as.numeric(names(tabNonZero))]<-as.vector(tabNonZero)
    
    opMu<-optim(par=c(rep(0.1,(n+1))),fn=optimizeParamsMuGroup,gr=NULL,
            countsVec=countsVec,n=n,lo=lo,model="ME",method="L-BFGS-B",
            lower=c(rep(1e-10,n),0),upper=c(rep(Inf,n),(1-(1e-15))),
            control=list(fnscale=-20))
    
    opNull<-optim(par=c(rep(0.1,n)),fn=optimizeParamsMuGroup,gr=NULL,
            countsVec=countsVec,n=n,lo=lo,model="Null",method="L-BFGS-B",
            lower=rep(1e-10,n),upper=rep(Inf,n),control=list(fnscale=-20))
    
    likeMu<-opMu$value
    likeNull<-opNull$value
    LRT<-2*likeMu-2*likeNull
    pvalueLRT<-pchisq(LRT, df = 1, ncp = 0, lower.tail = FALSE)
    
    return(l=list("opMu"=opMu,"opNull"=opNull,"countsVec"=countsVec,
                "genes"=genes,"LRT"=LRT,"pvalueLRT"=pvalueLRT))
}



###############################################################################
#' @title Tests all maximal cliques for mutual exclusivity
#'
#' @description \code{findSignifCliques} performs step 3 of the TiMEx 
#' procedure, namely tests all candidate maximal cliques for mutual 
#' exclusivity and reports the significant ones.
#' 
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' @param mcStruct list containing maximal cliques, as returned by
#' \code{\link{doMaxCliques}}
#' @param groupPvalue threshold for the corrected p-value of the groups, lower 
#' than which cliques are significant (real number between 0 and 1). Default 
#' is 0.1.
#' 
#' @details This function displays progress messages, namely the size of the 
#' clique currently being tested, and the number of cliques to test. 
#' 
#' Note that sequentially performing steps 1, 2, and 3 of the TiMEx procedure
#' (functions \code{\link{analyzePairs}}, \code{\link{doMaxCliques}}, and 
#' \code{\link{findSignifCliques}}) is equivalent to simply running the 
#' function \code{\link{TiMEx}}.
#' 
#' @return list consisting of:
#' \itemize{
#' \item{\code{genesSignif}} {list of significantly mutually exclusive groups, 
#' as gene names, sorted by corrected p-value. The list contains as many 
#' elements as identified lengths of groups. For example, 
#' \code{genesSignif[[2]]}
#' is a list containing the gene names of the significant groups of size 2. 
#' Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a matrix, in which rows represent gene names of 
#' significantly mutually exclusive groups.}
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
#' for mutual exclusivity, as returned by \code{\link{doMaxCliques}}}
#' 
#' \item{\code{matrix}} {input binary alteration matrix}
#' 
#' \item{\code{groupPvalue}} {input threshold for the corrected p-value, lower
#' than which cliques are significant}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso \code{\link{analyzePairs}} for step 1 of the TiMEx procedure;
#' \code{\link{doMaxCliques}} for step 2 of the TiMEx procedure;
#' the wrapper function \code{\link{TiMEx}} for combining these three steps, 
#' and identifying mutually exclusive groups in a  binary dataset with the 
#' TiMEx model. The data structures \code{\link{ovarianOutput.rda}},
#' \code{\link{breastOutput.rda}}, \code{\link{gbmDendrixOutput.rda}}, and
#' \code{\link{gbmMuexOutput.rda}} are examples of structures resulting after 
#' running TiMEx on large cancer datasets.
#' 
#' @examples
#' # First, test all pairs from the ovarian dataset for mutual exclusivity 
#' # (takes approximately 5 minutes).
#' data(ovarian)
#' \donttest{ovarianPairs<-analyzePairs(ovarian)}
#' 
#' # Second, identify all maximal cliques using the default thresholds
#' \donttest{ovarianMaxCliques<-doMaxCliques(ovarianPairs)}
#' 
#' # Then, test all maximal cliques for mutual exclusivity and report the 
#' # significant ones, based on a corrected p-value threshold of 0.1 (default).
#' \donttest{ovarianMEgroups<-findSignifCliques(ovarian,ovarianMaxCliques)}
#' 
#' @export
findSignifCliques<-function(mat,mcStruct,groupPvalue)
{
    
    # check inputs
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("need to specify a binary matrix as input")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    if (missing(mcStruct))
        stop("need to specify a list of identified maximal cliques as input")
    if (length(setdiff(c("detectedLengths","genesInCliques","idxInCliques",
                        "Mus","noMaxCliques","pairMu","pairPvalue"),
                        names(mcStruct)))!=0)
        stop("check the names of the elements of the input list")
    
    if (!(class(mcStruct$detectedLengths)=="numeric" || 
            class(mcStruct$detectedLengths)=="integer"))
        stop("the elements of the field 'detectedLengths' of the input list 
            need to be numeric")
    if (!sum(mcStruct$detectedLengths%%1)==0)
        stop("the elements of the field 'detectedLengths' of the input list 
            need to be positive integers")
    if (any(mcStruct$detectedLengths<=0))
        stop("the elements of the field 'detectedLengths' of the input list 
            need to be positive integers")
    
    if (!(length(mcStruct$idxInCliques)==length(mcStruct$detectedLengths)))
        stop("the fields of the input list detectedLengths and idxInCliques 
            do not match up in size")
    if (sum(is.na(match(unlist(mcStruct$idxInCliques),c(1:dim(mat)[2])))))
        stop("at least one index provided in the element 'idxInCliques' of 
            the input list is not part of the indices of the input genes")
    
    if (!(length(mcStruct$genesInCliques)==length(mcStruct$detectedLengths)))
        stop("the fields of the input list detectedLengths and genesInCliques 
            do not match up in size")
    if (sum(is.na(match(unlist(mcStruct$genesInCliques),colnames(mat))))>0)
        stop("at least one gene name provided in the element 'genesInCliques' 
            of the input list is not part of the input genes")
    
    if (length(setdiff(c("OrderedIdxInCliques","OrderedGenesInCliques"),
                    names(mcStruct$Mus)))!=0)
        stop("at least one of the fields 'OrderedIdxInCliques' 
            and 'OrderedGenesInCliques' of the element 'Mus' of the input 
            list is missing")
    if (sum(is.na(match(unlist(mcStruct$Mus$OrderedIdxInCliques),
                        c(1:dim(mat)[2])))))
        stop("at least one index provided in the element 'idxInCliques' 
            of the element 'Mus' of the input list is not part of the indices 
            of the input genes")
    if (!(length(mcStruct$Mus$OrderedIdxInCliques)==
            length(mcStruct$detectedLengths)))
        stop("the fields of the input list detectedLengths and 
            Mus$OrderedIdxInCliques do not match up in size")
    if (sum(is.na(match(unlist(mcStruct$Mus$OrderedGenesInCliques),
                        colnames(mat))))>0)
        stop("at least one gene name provided in the element 
            'OrderedGenesInCliques' of the element 'Mus' of the input list 
            is not part of the input genes")
    if (!(length(mcStruct$Mus$OrderedGenesInCliques)==
            length(mcStruct$detectedLengths)))
        stop("the fields of the input list detectedLengths and 
            Mus$OrderedGenesInCliques do not match up in size")
    
    if (!(class(mcStruct$noMaxCliques)=="numeric" || 
            class(mcStruct$noMaxCliques)=="integer"))
        stop("the field 'noMaxCliques' of the input list needs to be numeric")
    if (!sum(mcStruct$noMaxCliques%%1)==0)
        stop("the elements of the field 'noMaxCliques' of the input 
            list need to be positive integers")
    if (any(mcStruct$noMaxCliques<=0))
        stop("the elements of the field 'noMaxCliques' of the input 
            list need to be positive integers")
    
    if (!(class(mcStruct$pairMu)=="numeric"))
        stop("the field 'pairMu' of the input list should have been a real 
            number between 0 and 1")
    if (length(mcStruct$pairMu)>1)
        stop("the field 'pairMu' of the input list should have been a real 
            number between 0 and 1")
    if (!(0<=mcStruct$pairMu && mcStruct$pairMu<=1))
        stop("the field 'pairMu' of the input list should have been a real 
            number between 0 and 1")
    
    if (!(class(mcStruct$pairPvalue)=="numeric"))
        stop("the field 'pairPvalue' of the input list should have been a real 
            number between 0 and 1")
    if (length(mcStruct$pairPvalue)>1)
        stop("the field 'pairPvalue' of the input list should have been a real 
            number between 0 and 1")
    if (!(0<=mcStruct$pairPvalue && mcStruct$pairPvalue<=1))
        stop("the field 'pairPvalue' should have been a real 
            number between 0 and 1")
    
    if (missing(groupPvalue))
        groupPvalue<-0.1
    if (!(class(groupPvalue)=="numeric"))
        stop("the level of the corrected p value has to be a real number 
            between 0 and 1")
    if (length(groupPvalue)>1)
        stop("the level of the corrected p value has to be a real number 
            between 0 and 1")
    if (!(0<=groupPvalue && groupPvalue<=1))
        stop("the level of the corrected p value needs to be a real number 
            between 0 and 1")
    
    
    testedCand<-testAllCandidateGroups(mcStruct,mat)
    signifCliques<-filterSignifCliques(mcStruct,testedCand,groupPvalue)
    signifCliques$mcStruct<-mcStruct
    signifCliques$matrix<-mat
    signifCliques$groupPvalue<-groupPvalue
    return(signifCliques)  
}



###############################################################################
#' @title Finds mutually exclusive groups
#'
#' @description \code{TiMEx} is the main function of this package. It 
#' identifies all groups of mutually exclusive cancer alterations in a large 
#' binary input dataset.
#' 
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' @param pairMu pair-level threshold on mu (real number between 0 and 1). 
#' Default is 0.5.
#' @param pairPvalue pair-level threshold on p-value (real number between 0 
#' and 1). Default is 0.01.
#' @param groupPvalue threshold for the corrected p-value, lower than which 
#' cliques are significant. Default to 0.1
#' 
#' @details Dependning on the size of the dataset (both in terms of samples 
#' and alterations), TiMEx can require a reasonable time to run. For example, 
#' the approximate running time is 10 minutes for the 
#' \code{\link{ovarian.rda}} cancer dataset, and 45 minutes for the 
#' \code{\link{breast.rda}} cancer  dataset included in this package, on a 
#' personal computer.
#' 
#' \code{TiMEx} displays progress messages. In a fist step, it indicates
#' the gene which is currently being tested against the remaining genes. In a 
#' later step, it indicates the size of the clique currently being tested, 
#' and the number of cliques to test.
#' 
#' @return list consisting of:
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
#' \item{\code{idxSignif}} {list of significantly mutually exclusive groups, 
#' as indices in the input matrix, sorted by corrected p-value. The list 
#' contains  as many elements as identified lengths of groups. For example,
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
#' element is a vector, of length the number of significant maximal cliques of 
#' a given size.} 
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
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso \code{\link{analyzePairs}} for step 1 of the TiMEx procedure;
#' \code{\link{doMaxCliques}} for step 2 of the TiMEx procedure, and 
#' \code{\link{findSignifCliques}} for step 3 of the TiMEx procedure. The data
#' structures \code{\link{ovarianOutput.rda}}, \code{\link{breastOutput.rda}},
#' \code{\link{gbmDendrixOutput.rda}}, and \code{\link{gbmMuexOutput.rda}} 
#' are examples of structures resulting after running TiMEx on large cancer 
#' datasets.
#' 
#' @examples
#' # Run TiMEx on the ovarian cancer dataset with default parameters 
#' # (takes approximately 10 minutes)
#' data(ovarian)
#' \donttest{ovarianMEGroups<-TiMEx(ovarian)}
#' 
#' @export
TiMEx<-function(mat,pairMu,pairPvalue,groupPvalue)
{
    # check inputs
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("input needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    if (missing(pairMu))
        pairMu<-0.5
    if (!(class(pairMu)=="numeric"))
        stop("the threshold on mu needs to be a real number between 0 and 1")
    if (length(pairMu)>1)
        stop("the threshold on mu needs to be a real number between 0 and 1")
    if (!(0<=pairMu && pairMu<=1))
        stop("the threshold on mu needs to be a real number between 0 and 1")
    
    if (missing(pairPvalue))
        pairPvalue<-0.01
    if (!(class(pairPvalue)=="numeric"))
        stop("the threshold on p-value needs to be a real number between 
            0 and 1")
    if (length(pairPvalue)>1)
        stop("the threshold on p-value needs to be a real number between 
            0 and 1")
    if (!(0<=pairPvalue && pairPvalue<=1))
        stop("the threshold on p-value needs to be a real number between 
            0 and 1")
    
    if (missing(groupPvalue))
        groupPvalue<-0.1
    if (length(groupPvalue)>1)
        stop("the level of the corrected p v-alue has to be a real 
            number between 0 and 1")
    if (!(class(groupPvalue)=="numeric"))
        stop("the level of the corrected p v-alue has to be a real 
            number between 0 and 1")
    if (!(0<=groupPvalue && groupPvalue<=1))
        stop("the level of the corrected p v-alue needs to be a real 
            number between 0 and 1")
    
    pairs<-analyzePairs(mat)
    mcStruct<-doMaxCliques(pairs,pairMu,pairPvalue)
    signifCliques<-findSignifCliques(mat,mcStruct,groupPvalue)
    return(signifCliques)
}



###############################################################################
#' @title Creates metagroups of genes
#' 
#' @description \code{doMetagene} collapses genes with identical alteration
#' patterns across patients into metagroups. It returns a new matrix with
#' the collapsed genes, as well as the members of the metagroups.
#' 
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' 
#' @details It is recommended to run this function on the input binary matrix 
#' before applying TiMEx, because genes with identical alteration patterns 
#' across patients will otherwise be indistinguishable.
#' 
#' Note that in the datasets provided in this package, the genes 
#' have already been collapsed into metagroups.
#' 
#' @return List consisting of:
#' \itemize{
#' \item {\code{newMat}} {the collapsed input binary matrix, with metagenes
#' instead of genes.}
#' \item {\code{groups}} {list of metagenes, with as many elements as input 
#' genes which had an identical alteration pattern with at least one other 
#' input gene.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso \code{\link{ovarianGroups.rda}}, \code{\link{breastGroups.rda}} 
#' for examples of metagroups in large cancer datasets.
#' 
#' @examples
#' # Simulate genes and extract groups
#' simGenes<-simulateGenes(c(0.5,1,0.3),0.8,4000)
#' genes<-cbind(simGenes$genes,simGenes$genes)
#' genesNew<-doMetagene(genes)
#'
#' # In the datasets provided in this package, the genes have already been 
#' # collapsed into metagroups, hence the new matrix will be identical to the 
#' # old one.
#' data(ovarian)
#' \donttest{ovarianNew<-doMetagene(ovarian)}
#' 
#' @export
#'
doMetagene<-function(mat)
{
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("input needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    table.combos <- mat
    matcRow<-function(mine,table.combos)
    {
        row.is.a.match <- apply(table.combos, 2, identical, mine)
        match.idx <- which(row.is.a.match)
        return(match.idx)
    }
    
    ll<-apply(mat,2,matcRow,table.combos=mat)
    
    if (class(ll)=="matrix")
    {
        ll<-split(ll, rep(1:ncol(ll), each = nrow(ll)))  
    }
    
    #list with the groups of metagenes 
    groups<-ll[which(unlist(lapply(ll,length))!=1)]
    #the new alteration matrix with the identical genes collapsed into 
    #metagroups
    newMat<-mat[,which(duplicated(t(mat))==0)] 
    
    return(result=list("newMat"=newMat,"groups"=groups))
}



###############################################################################
#' @title Removes alterations based on frequency
#' 
#' @description \code{removeLowFreqs} returns a binary matrix from which genes 
#' altered with lower frequency than the input level are removed.
#' 
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' @param level frequency level under which the genes are to be removed 
#' (real number between 0 and 1). Default is 0.03.
#' 
#' @details It is only recommended to run this function if the user is not at 
#' all interested in the role of low-frequently altered genes.
#' 
#' @return the binary input matrix without the low-frequently altered genes.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso \code{\link{ovarian.rda}}, \code{\link{breast.rda}}, 
#' and \code{\link{gbmDendrix.rda}}  for examples of biological large cancer
#' datasets.
#' 
#' @examples
#' # Remove genes altered in less than 3% (the default level) of the samples  
#' # in the ovarian cancer dataset
#' data(ovarian)
#' ovarianNew<-removeLowFreqs(ovarian)
#' 
#' @export
#'
removeLowFreqs<-function(mat,level)
{
    # check inputs
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("input needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    if (missing(level))
        level<-0.03
    if (!(class(level)=="numeric"))
        stop("the level parameter needs to be a real number between 0 and 1")
    if (length(level)>1)
        stop("the level parameter needs to be a real number between 0 and 1")
    if (!(0<=level && level<=1))
        stop("the level parameter needs to be a real number between 0 and 1")
    
    nrpats<-dim(mat)[1]
    freqs<-apply(mat,2,sum)/nrpats
    small<-which(freqs<level)
    if (length(small)>0)
    {
        mat<-mat[,-small]
    }
    return(mat)
}



###############################################################################
#' @title Produces tables with groups
#'
#' @description \code{produceTablesSignifGroups} produces tables with the 
#' significant groups. These tables include the names of the genes part of the 
#' groups, their respective frequency in the dataset, and the mu and corrected
#' pvalue corresponding to each group. 
#' 
#' @param signifGroups result structure with the significant groups, as 
#' returned by either \code{\link{TiMEx}} or \code{\link{findSignifCliques}}
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' @param noToShow maximum number of groups to include in the table. Default 
#' is 30.
#' 
#' @details This function summarizes information on the mutually exclusive
#' groups identified by TiMEx in a dataset, as tables.
#' 
#' @return list with as many elements as lengths of the identified mutually
#' exclusive groups, containing tables with the significant groups for each 
#' size. 
#' Each list of this type further has two elements, \code{fdr} and 
#' \code{bonf}, corresponding to different multiple testing correction 
#' methods. Each element is a matrix, in which rows represent significantly 
#' mutually exclusive groups.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso the wrapper function \code{\link{TiMEx}} for identifying
#' mutually exclusive groups in a  binary dataset with the TiMEx model.
#' 
#' @examples
#' # Produce tables on the output of TiMEx on the ovarian cancer dataset
#' data(ovarian)
#' data(ovarianOutput)
#' ovarianTables<-produceTablesSignifGroups(ovarianOutput,ovarian)
#' 
#' @export
produceTablesSignifGroups<-function(signifGroups,mat,noToShow)
{
    # check inputs
    if (missing(signifGroups))
        stop("need to specify a result structure as returned by function 
            TiMEx as input")
    if (length(setdiff(c("idxSignif","mcStruct","genesSignif",
            "MusGroup","pvals"),names(signifGroups)))!=0)
        stop("check that the input list contains the following elements: 
            'idxSignif', 'mcStruct', 'genesSignif', 'MusGroup', and 'pvals' ")
    
    if (class(signifGroups$mcStruct)!="list")
        stop("the element 'mcStruct' of the input list should have been a 
            list")    
    if (is.null(signifGroups$mcStruct$detectedLengths))
        stop("the element 'mcStruct' of the input list should have contained 
            a field, 'detectedLengths'")
    if (!(class(signifGroups$mcStruct$detectedLengths)=="numeric" 
            || class(signifGroups$mcStruct$detectedLengths)=="integer"))
        stop("the elements of the 'detectedLengths' field of 
            'mcStruct' in the input list should have been positive integers")
    if (!sum(signifGroups$mcStruct$detectedLengths%%1)==0)
        stop("the elements of the 'detectedLengths' field of 
            'mcStruct' in the input list should have been positive integers")
    if (any(signifGroups$mcStruct$detectedLengths<=0))
        stop("the elements of the 'detectedLengths' field of 
            'mcStruct' in the input list should have been positive integers")
    
    if (is.null(signifGroups$mcStruct$noMaxCliques))
        stop("the element 'mcStruct' of the input list should 
            have contained a field, 'noMaxCliques'")
    if (!(class(signifGroups$mcStruct$noMaxCliques)=="numeric" 
            || class(signifGroups$mcStruct$noMaxCliques)=="integer"))
        stop("the elements of the 'noMaxCliques' field of 
            'mcStruct' in the input list should have been positive integers")
    if (!sum(signifGroups$mcStruct$noMaxCliques%%1)==0)
        stop("the elements of the 'noMaxCliques' field of 
            'mcStruct' in the input list should have been positive integers")
    if (any(signifGroups$mcStruct$noMaxCliques<=0))
        stop("the elements of the 'noMaxCliques' field of 
            'mcStruct' in the input list should have been positive integers")
    
    if (length(signifGroups$mcStruct$noMaxCliques)!=
            length(signifGroups$mcStruct$detectedLengths))
        stop("the lengths of the fields 'detectedLengths' and 
            'noMaxCliques' of 'mcStruct' in the input list do not coincide")
    
    if (class(signifGroups$idxSignif)!="list")
        stop("the field 'idxSignif' of the input list 
            should have been a non-empty list")
    if (length(signifGroups$idxSignif)==0)
        stop("the field 'idxSignif' of the input list 
            should have been a non-empty list")
    for (i in 1:length(signifGroups$idxSignif))
        if (!is.null(signifGroups$idxSignif[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$idxSignif[[i]])))!=0)
                stop("each element of the 'idxSignif' field of the input list 
                    should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
            
    if (class(signifGroups$genesSignif)!="list")
        stop("the field 'genesSignif' of the input list should have been 
            a non-empty list")
    if (length(signifGroups$genesSignif)==0)
        stop("the field 'genesSignif' of the input list should have been 
            a non-empty list")
    for (i in 1:length(signifGroups$genesSignif))
        if (!is.null(signifGroups$genesSignif[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$genesSignif[[i]])))!=0)
                stop("each element of the 'genesSignif' field of the 
                    input list should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
        
    
    if (class(signifGroups$pvals)!="list")
        stop("the field 'pvals' of the input list should have been 
            a non-empty list")
    if (length(signifGroups$pvals)==0)
        stop("the field 'pvals' of the input list should have been 
            a non-empty list")
    for (i in 1:length(signifGroups$pvals))
        if (!is.null(signifGroups$pvals[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$pvals[[i]])))!=0)
                stop("each element of the 'pvals' field of the input list 
                    should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
        
    
    if (class(signifGroups$MusGroup)!="list")
        stop("the field 'MusGroup' of the input list should have been 
            a non-empty list")
    if (length(signifGroups$MusGroup)==0)
        stop("the field 'MusGroup' of the input list should have been 
            a non-empty list")
    for (i in 1:length(signifGroups$MusGroup))
        if (!is.null(signifGroups$MusGroup[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$MusGroup[[i]])))!=0)
                stop("each element of the 'MusGroup' field of the input list 
                    should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
            
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("input needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    if (missing(noToShow))
        noToShow<-30
    if (length(noToShow)>1)
        stop("the number of groups to show needs to be a positive integer")
    if (!(class(noToShow)=="numeric"))
        stop("the number of groups to show needs to be a positive integer")
    if (!(noToShow%%1==0))
        stop("the number of groups to show needs to be a positive integer")
    if (noToShow<=0)
        stop("the number of groups to show needs to be a positive integer")
    
    
    ff<-getFreqsForCliquesAfterTest(signifGroups$idxSignif,mat)
    
    namesFreqs<-list()
    
    if (length(signifGroups$mcStruct$noMaxCliques)>=2)
    {
        for (i in 2:length(signifGroups$mcStruct$noMaxCliques))
        {
            namesFreqs[[i]]<-list()
            
            namesFreqs[[i]]$fdr<-combineNameWithFreq(
                signifGroups$genesSignif[[i]]$fdr,ff[[i]]$fdr,
                signifGroups$MusGroup[[i]]$fdr,signifGroups$pvals[[i]]$fdr)
            if (!is.null(dim(namesFreqs[[i]]$fdr)))
            {
                if (dim(namesFreqs[[i]]$fdr)[1]>noToShow)
                    namesFreqs[[i]]$fdr<-namesFreqs[[i]]$fdr[c(1:noToShow),]  
            }
            namesFreqs[[i]]$fdr<-as.data.frame(namesFreqs[[i]]$fdr)
            
            namesFreqs[[i]]$bonf<-combineNameWithFreq(
                signifGroups$genesSignif[[i]]$bonf,ff[[i]]$bonf,
                signifGroups$MusGroup[[i]]$bonf,signifGroups$pvals[[i]]$bonf)
            if (!is.null(dim(namesFreqs[[i]]$bonf)))
            {
                if (dim(namesFreqs[[i]]$bonf)[1]>noToShow)
                    namesFreqs[[i]]$bonf<-namesFreqs[[i]]$bonf[c(1:noToShow),]
            }
            namesFreqs[[i]]$bonf<-as.data.frame(namesFreqs[[i]]$bonf)
        }
    } else
    {
        i<-signifGroups$mcStruct$detectedLengths
        
        namesFreqs[[i]]<-list()
        
        namesFreqs[[i]]$fdr<-combineNameWithFreq(
            signifGroups$genesSignif[[i]]$fdr,ff[[i]]$fdr,
            signifGroups$MusGroup[[i]]$fdr,signifGroups$pvals[[i]]$fdr)
        if (!is.null(dim(namesFreqs[[i]]$fdr)))
        {
            if (dim(namesFreqs[[i]]$fdr)[1]>noToShow)
                namesFreqs[[i]]$fdr<-namesFreqs[[i]]$fdr[c(1:noToShow),]  
        }
        namesFreqs[[i]]$fdr<-as.data.frame(namesFreqs[[i]]$fdr)
        
        namesFreqs[[i]]$bonf<-combineNameWithFreq(
            signifGroups$genesSignif[[i]]$bonf,ff[[i]]$bonf,
            signifGroups$MusGroup[[i]]$bonf,signifGroups$pvals[[i]]$bonf)
        if (!is.null(dim(namesFreqs[[i]]$bonf)))
        {
            if (dim(namesFreqs[[i]]$bonf)[1]>noToShow)
                namesFreqs[[i]]$bonf<-namesFreqs[[i]]$bonf[c(1:noToShow),]  
        }
        namesFreqs[[i]]$bonf<-as.data.frame(namesFreqs[[i]]$bonf)
    }
    
    return(namesFreqs)
}


###############################################################################
#' @title Assesses the stability of groups by subsampling
#'
#' @description \code{subsampleAnalysis} subsamples the set of patients and 
#' assess the stability of the identified mutually exclusive groups.
#' 
#' @param subsampl a vector with subsampling frequencies
#' @param noReps number of repetitions of subsampling. Default is 100.
#' @param signifGroups result structure with the significant groups, as 
#' returned by either \code{\link{TiMEx}} or \code{\link{findSignifCliques}}
#' 
#' @details As this function runs TiMEx many times sequentially, it is 
#' computationally very intensive. For a version of this function which can be
#' directly ran on a cluster, please contact me (see e-mail below).
#' 
#' For example, after loading \code{data(ovarianOutput)} and running 
#' 
#' \code{counts<-subsampleAnalysis(subsampl=c(0.3,0.5,0.8),
#' signifGroups=signifGroups)}
#' 
#' \code{counts[[1]][[3]]} will represent
#' the relative counts of the identified mutually exclusive groups of size 3 
#' for a subsampling frequency of 30\%, for both  \code{fdr} and \code{bonf} 
#' (bonferroni) multiple correction methods.
#' 
#' @return list with as many elements as subsampling frequencies provided. Each
#' element is further a list with as many elemenets as number of sizes of the
#' significantly mutually exclusive groups identified. Aditionally, \code{bonf}
#' and \code{fdr} are two lists corresponding to each of these elements, 
#' representing different multiple correction methods. Finally, each element 
#' is a vector of subsampling frequencies of the significantly mutually 
#' exclusive groups identified. For an example, see \emph{Details} above. 
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso the wrapper function \code{\link{TiMEx}} for identifying mutually 
#' exclusive groups in a  binary dataset with the TiMEx model,
#' \code{\link{ovarianSubsampling.RData}}, 
#' \code{\link{breastSubsampling.RData}}, 
#' \code{\link{gbmDendrixSubsampling.RData}}, and 
#' \code{\link{gbmMuexSubsampling.RData}} for examples of outputs after
#' performing the subsampling analysis.
#' 
#' @examples
#' # running this function is time-intensive
#' data(ovarianOutput)
#' \dontrun{subsampleOvarian<-subsampleAnalysis(c(0.3,0.5,0.8),ovarianOutput)}
#' 
#' @export
subsampleAnalysis<-function(subsampl,noReps,signifGroups)
{

    # check inputs
    if (missing(subsampl))
        stop("a vector of subsampling frequencies is needed as input")
    if (!(class(subsampl)=="numeric"))
        stop("the subsampling frequencies need to be real numbers between 
            0 and 1")
    if (sum(sapply(subsampl,function(x){if (!(x>0 && x<=1)) p<-1 else p<-0; 
            return(p)}))>0)
        stop("the subsampling frequencies need to be real numbers between 
            0 and 1")
    
    if (missing(noReps))
        noReps<-100
    if (!(class(noReps)=="numeric"))
        stop("the number of repetitions needs to be positive integer")
    if (!(noReps%%1==0))
        stop("the number of repetitions needs to be positive integer")
    if (!noReps>0)
        stop("the number of repetitions needs to be positive integer")
    
    if (missing(signifGroups))
        stop("need to specify a result structure as returned by 
            function TiMEx as input")
    if (length(setdiff(c("idxSignif","mcStruct","genesSignif",
            "MusGroup","pvals"),names(signifGroups)))!=0)
        stop("check that the input list contains the following elements: 
            'idxSignif', 'mcStruct', 'genesSignif', 'MusGroup', and 'pvals' ")
    
    if (class(signifGroups$mcStruct)!="list")
        stop("the element 'mcStruct' of the input list should have been a 
            list")    
    if (is.null(signifGroups$mcStruct$detectedLengths))
        stop("the element 'mcStruct' of the input list should have contained 
            a field, 'detectedLengths'")
    if (!(class(signifGroups$mcStruct$detectedLengths)=="numeric" 
        || class(signifGroups$mcStruct$detectedLengths)=="integer"))
        stop("the elements of the 'detectedLengths' field of 
            'mcStruct' in the input list should have been positive integers")
    if (!sum(signifGroups$mcStruct$detectedLengths%%1)==0)
        stop("the elements of the 'detectedLengths' field of 
            'mcStruct' in the input list should have been positive integers")
    if (any(signifGroups$mcStruct$detectedLengths<=0))
        stop("the elements of the 'detectedLengths' field of 
            'mcStruct' in the input list should have been positive integers")
    
    if (is.null(signifGroups$mcStruct$noMaxCliques))
        stop("the element 'mcStruct' of the input list 
            should have contained a field, 'noMaxCliques'")
    if (!(class(signifGroups$mcStruct$noMaxCliques)=="numeric" 
        || class(signifGroups$mcStruct$noMaxCliques)=="integer"))
        stop("the elements of the 'noMaxCliques' field of 
            'mcStruct' in the input list should have been positive integers")
    if (!sum(signifGroups$mcStruct$noMaxCliques%%1)==0)
        stop("the elements of the 'noMaxCliques' field of 
            'mcStruct' in the input list should have been positive integers")
    if (any(signifGroups$mcStruct$noMaxCliques<=0))
        stop("the elements of the 'noMaxCliques' field of 
            'mcStruct' in the input list should have been positive integers")
    
    if (length(signifGroups$mcStruct$noMaxCliques)!=
            length(signifGroups$mcStruct$detectedLengths))
        stop("the lengths of the fields 'detectedLengths' and 
            'noMaxCliques' of 'mcStruct' in the input list do not coincide")
    
    if (class(signifGroups$idxSignif)!="list")
        stop("the field 'idxSignif' of the input list should have been a 
            non-empty list")
    if (length(signifGroups$idxSignif)==0)
        stop("the field 'idxSignif' of the input list should have been a 
            non-empty list")
    for (i in 1:length(signifGroups$idxSignif))
        if (!is.null(signifGroups$idxSignif[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),names(signifGroups$
                                                idxSignif[[i]])))!=0)
                stop("each element of the 'idxSignif' 
            field of the input list should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
    
    if (class(signifGroups$genesSignif)!="list")
        stop("the field 'genesSignif' of the input list should have been 
            a non-empty list")
    if (length(signifGroups$genesSignif)==0)
        stop("the field 'genesSignif' of the input list should have been 
            a non-empty list")
    for (i in 1:length(signifGroups$genesSignif))
        if (!is.null(signifGroups$genesSignif[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$genesSignif[[i]])))!=0)
                stop("each element of the 'genesSignif' field of the input 
                    list should have been a list with two elements, 'fdr' and 
                    'bonf'")
        }
    
    
    if (class(signifGroups$pvals)!="list")
        stop("the field 'pvals' of the input list should have been a non-empty 
            list")
    if (length(signifGroups$pvals)==0)
        stop("the field 'pvals' of the input list should have been a non-empty 
            list")
    for (i in 1:length(signifGroups$pvals))
        if (!is.null(signifGroups$pvals[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$pvals[[i]])))!=0)
                stop("each element of the 'pvals' field of the input list 
                    should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
    
    
    if (class(signifGroups$MusGroup)!="list")
        stop("the field 'MusGroup' of the input list should have 
            been a non-empty list")
    if (length(signifGroups$MusGroup)==0)
        stop("the field 'MusGroup' of the input list should have 
            been a non-empty list")
    for (i in 1:length(signifGroups$MusGroup))
        if (!is.null(signifGroups$MusGroup[[i]]))
        {
            if (length(setdiff(c("fdr","bonf"),
                            names(signifGroups$MusGroup[[i]])))!=0)
                stop("each element of the 'MusGroup' field of the input list 
                    should have been a list with two elements, 
                    'fdr' and 'bonf'")
        }
    
    if (length(setdiff(c("pairMu","pairPvalue"),
                    names(signifGroups$mcStruct)))>0)
        stop("the 'mcStruct' field of the input list should have included 
            the parameters pairMu and pairPvalue")
    if (!(class(signifGroups$mcStruct$pairPvalue)=="numeric"))
        stop("the field 'pairPvalue' of 'mcStruct' of the input list should 
            have been a real number between 0 and 1")
    if (!(0<=signifGroups$mcStruct$pairPvalue && 
            signifGroups$mcStruct$pairPvalue<=1))
        stop("the field 'pairPvalue' of 'mcStruct' of the input list should 
            have been a real number between 0 and 1")
    if (length(signifGroups$mcStruct$pairPvalue)!=1)
        stop("the field 'pairPvalue' of 'mcStruct' of the input list should 
            have been a real number between 0 and 1")
    if (!(class(signifGroups$mcStruct$pairMu)=="numeric"))
        stop("the field 'pairMu' of 'mcStruct' of the input list should 
            have been a real number between 0 and 1")
    if (!(0<=signifGroups$mcStruct$pairMu && signifGroups$mcStruct$pairMu<=1))
        stop("the field 'pairMu' of 'mcStruct' of the input list should 
            have been a real number between 0 and 1")
    if (length(signifGroups$mcStruct$pairMu)!=1)
        stop("the field 'pairMu' of 'mcStruct' of the input list should 
            have been a real number between 0 and 1")
    
    if (!(class(signifGroups$groupPvalue)=="numeric"))
        stop("the field 'groupPvalue' of the input list should 
            have been a real number between 0 and 1")
    if (!(0<=signifGroups$groupPvalue && signifGroups$groupPvalue<=1))
        stop("the field 'groupPvalue' of the input list should have been 
            a real number between 0 and 1")
    if (length(signifGroups$groupPvalue)!=1)
        stop("the field 'groupPvalue' of the input list should have been 
            a real number between 0 and 1")

    if (is.null(signifGroups$matrix))
        stop("the input list should have had a binary matrix, 'matrix', 
            as one of the fields")
    if (is.null(dim(signifGroups$matrix)))
        stop("the input list should have had a binary matrix, 'matrix', 
            as one of the fields")
    if ((all.equal(c(0,1),sort(unique(as.vector(signifGroups$matrix))))!=1))
        stop("the input list should have had a binary matrix, 'matrix', 
            as one of the fields")
    
    if (sum(is.na(match(unlist(signifGroups$genesSignif),
                        colnames(signifGroups$matrix))))>0)
        stop("at least one gene name provided in the element 'genesInCliques' 
            of the input list is not part of the input genes")
    
    
    pairMu<-signifGroups$mcStruct$pairMu
    pairPvalue<-signifGroups$mcStruct$pairPvalue
    groupPvalue<-signifGroups$groupPvalue
    
    countsAll<-list()
    
    if (length(signifGroups$genesSignif)>0)
    {
        for (idxSubs in 1:length(subsampl))
        {
            countsFreq<-list()
            
            for (io in 1:length(signifGroups$mcStruct$detectedLengths))
            {
                countsFreq[[signifGroups$mcStruct$detectedLengths[io]]]<-list()
                countsFreq[[signifGroups$mcStruct$detectedLengths[io]]]$fdr<-
                    rep(0,max(1,dim(signifGroups$genesSignif
                        [[signifGroups$mcStruct$detectedLengths[io]]]$fdr)[1]))
                countsFreq[[signifGroups$mcStruct$
                    detectedLengths[io]]]$bonf<-
                    rep(0,max(1,dim(signifGroups$genesSignif
                    [[signifGroups$mcStruct$detectedLengths[io]]]$bonf)[1]))
            }
            
            for (idxR in 1:noReps)
            {
                matrixSubsample<-
                signifGroups$matrix[sample(c(1:dim(signifGroups$matrix)[1]),
                    ceiling(subsampl[idxSubs]*dim(signifGroups$matrix)[1])),]
                signifGroupsSubsample<-TiMEx(matrixSubsample,pairMu,pairPvalue,
                                            groupPvalue)
                
                genesSignifCombined<-list()
                
                for (io in 1:length(signifGroups$mcStruct$detectedLengths))
                {
                    genesSignifCombined[[signifGroups$mcStruct$
                                            detectedLengths[io]]]<-list()
                    
                    if (signifGroups$mcStruct$detectedLengths[io]>1)  
                    {
                        # if groups of the tested length were detected at all
                        if (length(signifGroupsSubsample$genesSignif)>=
                                signifGroups$mcStruct$detectedLengths[io]) 
                        {
genesSignifCombined[[signifGroups$mcStruct$detectedLengths[io]]]$bonf<-
    rbind(signifGroups$genesSignif[[signifGroups$mcStruct$
        detectedLengths[io]]]$bonf,signifGroupsSubsample$
            genesSignif[[signifGroups$mcStruct$detectedLengths[io]]]$bonf)
dph<-duplicated(genesSignifCombined[[signifGroups$mcStruct$
    detectedLengths[io]]]$bonf) | duplicated(genesSignifCombined[[signifGroups$
        mcStruct$detectedLengths[io]]]$bonf, fromLast = TRUE)
# if more than one group was detected initially
if (length(dim(signifGroups$genesSignif[[signifGroups$mcStruct$
        detectedLengths[io]]]$bonf)[1])>0)
    countsFreq[[signifGroups$mcStruct$
        detectedLengths[io]]]$bonf<-countsFreq[[signifGroups$mcStruct$
            detectedLengths[io]]]$bonf+dph[1:dim(signifGroups$
                genesSignif[[signifGroups$mcStruct$
                    detectedLengths[io]]]$bonf)[1]]+0
else
    countsFreq[[signifGroups$mcStruct$
        detectedLengths[io]]]$bonf<-countsFreq[[signifGroups$mcStruct$
            detectedLengths[io]]]$bonf+1
                            
genesSignifCombined[[signifGroups$mcStruct$detectedLengths[io]]]$
    fdr<-rbind(signifGroups$genesSignif[[signifGroups$mcStruct$
        detectedLengths[io]]]$fdr,signifGroupsSubsample$
            genesSignif[[signifGroups$mcStruct$detectedLengths[io]]]$fdr)      
dph<-duplicated(genesSignifCombined[[signifGroups$mcStruct$
    detectedLengths[io]]]$fdr) | duplicated(genesSignifCombined[[signifGroups$
        mcStruct$detectedLengths[io]]]$fdr, fromLast = TRUE)
# if more than one group was detected initially
if (length(dim(signifGroups$genesSignif[[signifGroups$mcStruct$
        detectedLengths[io]]]$fdr)[1])>0)
    countsFreq[[signifGroups$mcStruct$
        detectedLengths[io]]]$fdr<-countsFreq[[signifGroups$mcStruct$
            detectedLengths[io]]]$fdr+dph[1:dim(signifGroups$
                genesSignif[[signifGroups$mcStruct$detectedLengths[io]]]$
                    fdr)[1]]+0
else
    countsFreq[[signifGroups$mcStruct$
        detectedLengths[io]]]$fdr<-countsFreq[[signifGroups$mcStruct$
            detectedLengths[io]]]$fdr+1          
                        }  
                    } 
                } 
            }
            
            for (io in 1:length(signifGroups$mcStruct$detectedLengths))
            {
countsFreq[[signifGroups$mcStruct$detectedLengths[io]]]$
    bonf<-countsFreq[[signifGroups$mcStruct$detectedLengths[io]]]$bonf/noReps
countsFreq[[signifGroups$mcStruct$detectedLengths[io]]]$
    fdr<-countsFreq[[signifGroups$mcStruct$detectedLengths[io]]]$fdr/noReps
            }
            
            countsAll[[idxSubs]]<-countsFreq
        }
        names(countsAll)<-paste("SubsamplingFreq:",subsampl,sep="")
        for (jo in 1:length(countsAll))
        {
            names(countsAll[[jo]])<-paste("groupSize=",
                                    c(1:length(countsAll[[jo]])),sep="")
        }
    }
    return(countsAll)
}



###############################################################################
#' @title Plots a Mutually Exclusive group
#'
#' @description \code{plotGroupByName} plots a mutually exclusive group 
#' (including the frequencies of the genes), given the names of the genes in 
#' the group.
#' 
#' @param group vector of gene names to be plotted
#' @param mat binary alteration matrix, with rows representing patients and
#' columns representing genes
#' 
#' @details The plotting is done based on the function 
#' \code{\link[graphics]{image}}.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso the wrapper function \code{\link{TiMEx}} for identifying
#' mutually exclusive groups in a  binary dataset with the TiMEx model.
#' 
#' @examples
#' # Plot the group consisting of the copy number aberrations of MIEN1 and 
#' # CDKN1B, and the point mutations of CDH1, GATA3, and MAP3K1, in breast 
#' # cancer.
#' data(breast)
#' group<-c("MIEN1-CNA","CDH1-Mut","GATA3-Mut","MAP3K1-Mut","CDKN1B-CNA")
#' plotGroupByName(group,breast)
#' 
#' @export
plotGroupByName<-function(group,mat)
{
    # check inputs
    if (missing(group))
        stop("a group of genes needs to be specified")
    
    if (missing(mat))
        stop("need to specify a binary matrix as input")
    if (is.null(dim(mat)))
        stop("input needs to be a binary matrix")
    if ((all.equal(c(0,1),sort(unique(as.vector(mat))))!=1))
        stop("input needs to be a binary matrix")
    
    if ( sum(is.na(match(group,colnames(mat))))>0)
        stop("at least one gene name is not part of the input genes")
    
    l<-length(group)
    idx<-c()
    for (i in 1:l)
    {
        idx[i]<-which(colnames(mat)==group[i])
    }
    submatrix<-mat[,idx]
    freqs<-round(apply(submatrix,2,sum)*100/dim(submatrix)[1],2)
    ord<-order(freqs,decreasing=TRUE)
    submatrix<-submatrix[,ord]
    
    mat2<-submatrix
    for (i in l:1)
        mat2<-mat2[order(mat2[,i]),]
    
    image(t(mat2),col=c("white","black"),yaxt="n",xaxt="n",ylab="samples")
    axis(1,at=seq(from=0,to=1,length.out=l),tck=-0.03,cex.axis=0.7,
        labels=paste(group[ord],"\n ",freqs[ord]," %",sep=""))
}



###############################################################################
#' @title Recovers members of the metagroups
#'
#' @description \code{recoverAllNameGroups}  recovers, from a metagroup, the 
#' names of the genes part of an identified mutually exclusive group.
#' 
#' @param groupsMeta: list containing groups of equivalent genes, as returned 
#' by the field \code{groups} of \code{\link{doMetagene}}
#' @param clGenes: matrix of mutually exclusive groups of same size, as gene 
#' names. This type of matrix is returned by either \code{\link{TiMEx}} or
#' \code{\link{findSignifCliques}}, as one of the 
#' matrix elements of the \code{genesSignif} field. 
#' 
#' @details This function can be used if the input binary matrix contains 
#' identical events that need to be merged into metagenes using 
#' \code{\link{doMetagene}}. Running \code{recoverAllNamesGroups} provides the
#' set of identical alterations which are part of the identified mutually 
#' exclusive groups. 
#' 
#' In order to run this function on all the identified mutually exclusive 
#' groups as returned by \code{\link{TiMEx}} or 
#' \code{\link{findSignifCliques}}, it is necessary to run it separately on 
#' each matrix element (corresponding to different group sizes and different 
#' correction methods) of the \code{genesSignif} field in the structure 
#' returned by either \code{\link{TiMEx}} or \code{\link{findSignifCliques}}.
#' 
#' For example, after loading \code{data(ovarianGroups)} and 
#' \code{data(ovarianOutput)}, and running 
#' 
#' \code{rGroups<-recoverAllNamesGroups(ovarianGroups,
#' signifGroups$genesSignif[[3]]$bonf)}
#' 
#' \code{rGroups[[14]]} has 3 elements (as many as genes part
#' of the identified mutually exclusive groups). Each element is the metagroup 
#' of each of the genes part of the 14th mutually exclusive group in the input
#' matrix 
#' \code{signifGroups$genesSignif[[3]]$bonf}. Namely, \emph{BRD4-CNA} and 
#' \emph{MYC-CNA} have unique alteration patterns among samples, and are alone 
#' in their metagroup, while \emph{CASC1-CNA} has an identical alteration
#' pattern with \emph{KRAS-CNA} and \emph{LYRM5-CNA}. The numbers
#' below the gene names are the indices of the genes in the initial input 
#' binary matrix of patients.
#' 
#' @return list with as many elements as number of identified mutually 
#' exclusive groups, \emph{i.e.} number of rows in the input matrix. Each of 
#' its elements is further a list, containing, at each position, the metagroup 
#' of the genes in the initial group at that resepective position. For an 
#' example, see \emph{Details} above.
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso \code{\link{doMetagene}} for collapsing the genes of an input 
#' matrix with identical alteration patterns into metagroups.
#' 
#' @examples
#' data(ovarianGroups)
#' data(ovarianOutput)
#' r<-recoverAllNamesGroups(ovarianGroups,ovarianOutput$genesSignif[[3]]$bonf)
#' 
#' @export
recoverAllNamesGroups<-function(groupsMeta,clGenes)
{
    # check inputs
    if (missing(groupsMeta))
        stop("a list of groups needs to be provided as input")
    if (class(groupsMeta)!="list" && class(groupsMeta)!="character")
        stop("a list of groups needs to be provided as input")
    
    if (missing(clGenes))
        stop("a matrix of identified groups, as gene names, needs to be 
            provided as input")
    if (class(clGenes)!="matrix" && class(clGenes)!="character")
        stop("a matrix of identified groups, as gene names, needs to be 
            provided as input")
    
    newGroups<-list()
    for (i in 1:dim(clGenes)[1])
    {
        newGroups[[i]]<-list()
        k<-1
        for (j in 1:dim(clGenes)[2])
        {
            pos<-which(names(groupsMeta)==clGenes[i,j])
            if (length(pos)>0)
                newGroups[[i]][[k]]<-groupsMeta[pos]
            else
                newGroups[[i]][[k]]<-clGenes[i,j]
            k<-k+1
        }
    }
    return(newGroups)
}



###############################################################################
#' @title Generates data from the TiMEx model
#' 
#' @description \code{simulateGenes} returns a list containing a binary matrix 
#' simulated  from the TiMEx model, for given lambdas (exponential rates), 
#' mu (intensity of mutual exclusivity), and N (sample size).
#'   
#' @param lambdas vector of exponential rates (positive real numbers). 
#' The length of the vector equals the number of simulated genes.
#' @param mu intensity of mutual exclusivity (real number between 
#' 0 and 1). Default is 1 (perfect mutual exclusivity).
#' @param N sample size (positive integer).
#' 
#' @details This function needs \code{\link[gtools]{permutations}} in order to 
#' run.
#' 
#' For details on how the values of the exponential rates correspond to 
#' frequencies, see \emph{References} below.
#' 
#' @return list consisting of:
#' \itemize{
#' \item{\code{genes}} {the simulated dataset, as a binary matrix.}
#' \item{\code{genoMat}} {the matrix with genotype probabilities, from which 
#' the dataset was simulated. This matrix has as many dimensions as number of 
#' genes, i.e. 2x2x...x2. For each dimension, the first position corresponds to
#' the probability of observing a 0 for that gene, and the second position 
#' corresponds to the probability of observing an 1. For example, in the case 
#' of 4 genes, the probability of observing the null genotype \code{(0000)} is 
#' given by \code{genoMat[1,1,1,1]}; the probability of observing the genotype 
#' \code{(1011)} is given by \code{genoMat[2,1,2,2]}; the probability of 
#' observing the genotype \code{(1111)} is given by \code{genoMat[2,2,2,2]}. 
#' The entries of this matrix are nonnegative and sum up to 1.}
#' }
#' 
#' @author Simona Constantinescu, \email{simona.constantinescu@@bsse.ethz.ch}
#' 
#' @references "TiMEx: A Waiting Time Model For Mutually
#' Exclusive Cancer Alterations", by Constantinescu \emph{et al.} 
#' (Bioinformatics, 2015).
#' 
#' @seealso the wrapper function \code{\link{TiMEx}} for identifying
#' mutually exclusive groups in a  binary dataset with the TiMEx model,
#' \code{\link{ovarian.rda}}, \code{\link{breast.rda}}, and
#' \code{\link{gbmDendrix.rda}}  for examples of biological large cancer 
#' datasets.
#' 
#' @examples
#' simGenes<-simulateGenes(c(0.5,1,0.3),0.8,4000)
#' 
#' @importFrom gtools permutations
#' 
#' @export
#' 
simulateGenes<-function(lambdas,mu,N)
{
    # check inputs
    if (missing(lambdas))
        stop("a vector of nonegative exponential rates needs to be 
            provided as input")
    if (!(class(lambdas)=="numeric"))
        stop("a vector of nonegative exponential rates needs to be 
            provided as input")
    if (any(lambdas<=0))
        stop("a vector of nonegative exponential rates needs to be 
            provided as input")
    
    if (missing(mu))
        mu<-1
    if (class(mu)!="numeric")
        stop("mu needs to be a real number between 0 and 1")
    if (!(mu>=0 && mu<=1))
        stop("mu needs to be a real number between 0 and 1")
    
    if (missing(N))
        stop("the sample size N needs to be a positive integer")
    if (!(N%%1==0))
        stop("the sample size N needs to be a positive integer")
    if (N<=0)
        stop("the sample size N needs to be a positive integer")
    
    
    # set the observation lambda lo=1
    lo<-1
    
    n<-length(lambdas)
    
    genoMat<-array(NA,dim=rep(2,n))
    
    for (i in 1:n)
    {
        dimnames(genoMat)[[i]]<-paste("pos",i,"=",c(0,1),sep="")
    }
    
    l<-list()
    for (i in 1:n)
    {
        l[[i]]<-seq(0,1)
    }
    
    allGenos<-as.matrix(expand.grid(l))
    allGenosAtLeastTwo<-allGenos[which(apply(allGenos,1,sum)>1),]
    allGenosOne<-allGenos[which(apply(allGenos,1,sum)==1),]
    allGenosZero<-allGenos[which(apply(allGenos,1,sum)==0),]
    allGenoPositive<-rbind(allGenosOne,allGenosAtLeastTwo)
    
    # for the full0 case
    p1<-compProbFullZeroSims(allGenosZero,lambdas,lo) 
    # for all the genotypes which are not full 0
    p2<-apply(allGenoPositive,1,compAllProbsSims,lambdas=lambdas,lo=lo,mu=mu) 
        
    genoMat[1]<-p1
    genoMat[(allGenoPositive+1)]<-p2
    
    simGenes<-produceGenesGroup(N,genoMat,n)
    genes<-simGenes$genes
    
    return(list("genes"=genes,"genoMat"=genoMat))  
}