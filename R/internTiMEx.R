# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
computeContTable<-function(genes)
{
    # (internal) function to compute a contingency table from a matrix of two 
    # genes

    # inputs: 
        # genes: binary matrix whose columns are the two genes

    # outputs:
        # ns: list consisting of the entries in the contingency table

    gene1<-genes[,1]
    gene2<-genes[,2]
    sum<-gene1+gene2
    diff<-gene1-gene2
    n00<-length(which(sum==0))
    n11<-length(which(sum==2))
    n01<-length(which(diff==-1))
    n10<-length(which(diff==1))
    ns = list("n00"=n00,"n01"=n01,"n10"=n10,"n11"=n11)
    return(ns)
}



###############################################################################
computeParmsMu<-function(params,n)
{
    # (internal) helper numerical optimization function for the specific case 
    # n=2, for the mutual exclusivity model

    # inputs:
        # params: vector of three position, consisting of lambda_i, lambda_j, 
            # and mu
        # n: vector of four positions, consisting of the contingency table 
            # characterizing the two genes (n00,n01,n10,n11)

    # outputs:
        # l: log likelihood corresponding to the input contingency table

    lamio<-params[1]
    lamjo<-params[2]
    mu<-params[3]
    n00<-n[1]
    n01<-n[2]
    n10<-n[3]
    n11<-n[4]
    l<-(n01+n11)*log(lamjo)+(n10+n11)*log(lamio)-(n00+n01+n10+n11)*
        log(lamio+lamjo+1)+n01*log(1+mu*lamio)+n10*log(1+mu*lamjo)+
        n11*log(1-mu)-(n01+n11)*log(1+lamio)-(n10+n11)*log(1+lamjo)+
        n11*log(lamio+lamjo+2);
    return(l)
}



###############################################################################
computeParmsNull<-function(params,n)
{
    # (internal) helper numerical optimization function for the specific case 
    # n=2, for the null model

    # inputs:
        # params: vector of three position, consisting of lambda_i, lambda_j, 
            # and mu
        # n: vector of four positions, consisting of the contingency table 
            # characterizing the two genes (n00,n01,n10,n11)

    # outputs:
        # l: log likelihood corresponding to the input contingency table
    
    lamio<-params[1]
    lamjo<-params[2]
    n00<-n[1]
    n01<-n[2]
    n10<-n[3]
    n11<-n[4]
    l<-(n01+n11)*log(lamjo)+(n10+n11)*log(lamio)-(n00+n01+n10+n11)*
        log(lamio+lamjo+1)-(n01+n11)*log(1+lamio)-(n10+n11)*log(1+lamjo)+
        n11*log(lamio+lamjo+2);
    return(l)
}



###############################################################################
makeSym<-function(mat)
{
    # (internal) function to symmetrize a matrix

    # inputs: 
        # mat: matrix to be symmetrized

    # outputs:
        # mat: symmetrized matrix
    
    ind<-lower.tri(mat)
    mat[ind]<-t(mat)[ind]
    return(mat)
}



###############################################################################
doClique<-function(muEstSym,pvalueLRTCorrectSym,pairMu,pairPvalue)
{
    # (internal) function to compute all maximal cliques

    # inputs: 
        # muEstSym: symetrical matrix of estimates of mu (as returned by 
            # function analyzePairs)
        # pvalueLRTCorrectSym: list of pairwise pvalues as matrices (as 
            # returned by function analyzePairs)
        # pairMu: pair-level threshold on mu (greater equal)
        # pairPvalue: pair-level threshold on pvalue (smaller equal)
    
    # outputs: 
        # result: list consisting of:
            # adj: input binary matrix satisfying the mu and pvalue thresholds
            # cliques: identified cliques on the basis of the matrix 'adj'
            # lengths: lengths of the identified cliques
            # genes: names of the input genes
    
    d<-dim(muEstSym)[1]
    adj<-matrix(0,nrow=d,ncol=d)
    allowed<-setdiff(which(pvalueLRTCorrectSym<=pairPvalue),
                    which(muEstSym<pairMu))
    adj[allowed]<-1
    
    gg<-graph.adjacency(adj,mode="undirected")
    pp<-igraph.to.graphNEL(gg)
    v<-maxClique(pp)
    v$maxCliques<-lapply(v$maxCliques,sort)
    ll<-unlist(lapply(v$maxCliques,length))
    genes<-colnames(muEstSym)
    
    return(result=list("adj"=adj,"cliques"=v,"lengths"=ll,"genes"=genes))
}



###############################################################################
analyzeMaxCliquesBySize<-function(cliqueStruct)
{
    # (internal) function to organize the identified maximal cliques by size, 
    # as well as return cleaner structures

    # inputs: 
        # cliqueStruct: structure containing cliques, as returned by the 
            # function doClique

    # outputs: 
        # result: list consisting of:
            # noMaxCliques: vector of unique lengths of the identified maximal 
                # cliques
            # idxInCliques: list with number of elements the lengths of 
                # detected maximal cliques (starts with 2). each element is a 
                # matrix containing the indices of all members of the maximal 
                # cliques of that respective size
            # genesInCliques: list with same structure as idxInCliques, 
                # containing names of genes instead of indices

    result<-list()

    detectedLengths<-sort(unique(cliqueStruct$lengths))
    result$detectedLengths<-detectedLengths

    idxInCliques<-list()
    genesInCliques<-list()
    for (k in detectedLengths)
    {
        clg<-which(cliqueStruct$lengths==k)
        celem<-matrix(unlist(cliqueStruct$cliques$maxCliques[clg]),ncol=k,
                    byrow=TRUE)
        idxInCliques[[k]]<-celem
        class(idxInCliques[[k]])<-"numeric"
        if (is.null((apply(idxInCliques[[k]],1,function(x)
            {cliqueStruct$genes[x]}))))
            genesInCliques[[k]]<-apply(idxInCliques[[k]],1,function(x)
                {cliqueStruct$genes[x]}) else
                genesInCliques[[k]]<-t(apply(idxInCliques[[k]],1,function(x)
                    {cliqueStruct$genes[x]}))
        
                
    }

    noMaxCliques<-unlist(lapply(genesInCliques,function(x){dim(x)[1]}))

    result$idxInCliques<-idxInCliques
    result$genesInCliques<-genesInCliques
    result$noMaxCliques<-noMaxCliques

    return(result)
}



###############################################################################
getMusAllMaxCliques<-function(maxCliqueStruct,matrixMus)
{
    # (internal) function to compute the total weight of the cliques (sum of 
    # estimated mus on the edges)

    # inputs: 
        # maxCliqueStruct, as returned by function analyzeMaxCliquesBySize
        # matrixMus: matrix of estimated pairwise Mus, as returned by function 
            # analyzePairs

    # outputs: 
        # maxCliqueStruct: list consisting of (in addition to the inputs):
            # OrderedGenesInCliques: same as GenesInCliques (part of intput), 
                # just ordered decreasingly by weight of maximal cliques
            # OrderedIdxInCliques: same as IdxInCliques (part of input), just 
                # ordered decreasingly by weight of maximal cliques

    MusCliques<-list() # list of Mu values for the maximal cliques
    weight<-list() # total weight (sum of Mus) for each clique
    avgWeight<-list() # average weight (average Mu) for each clique
    
    MusCliques[[1]]<-NA
    weight[[1]]<-NA
    avgWeight[[1]]<-NA
    if (length(maxCliqueStruct$genesInCliques)>=2)
    {
        for (k in 2:length(maxCliqueStruct$genesInCliques))
        {
        MusCliques[[k]]<-list()
        weight[[k]]<-c(NA)
        avgWeight[[k]]<-c(NA)
        if (!(is.null(dim(maxCliqueStruct$genesInCliques[[k]])[1])))
        {
            for (i in 1:dim(maxCliqueStruct$genesInCliques[[k]])[1])
            {
    MusCliques[[k]][[i]]<-matrixMus[maxCliqueStruct$idxInCliques[[k]][i,],
                            maxCliqueStruct$idxInCliques[[k]][i,]]
    weight[[k]][i]<-
        sum(MusCliques[[k]][[i]][upper.tri(MusCliques[[k]][[i]])])
    avgWeight[[k]][i]<-weight[[k]][i]/(dim(
        maxCliqueStruct$idxInCliques[[k]])[2]*(dim(
            maxCliqueStruct$idxInCliques[[k]])[2]-1)/2)
            } 
        } 
        }
    }

    ordWeight<-list()
    ordWeight[[1]]<-NA

    if (length(maxCliqueStruct$genesInCliques)>=2)
    {
        for (k in 2:length(maxCliqueStruct$genesInCliques))
        {
        ordWeight[[k]]<-order(weight[[k]],decreasing=TRUE)

        if (!(is.null(dim(maxCliqueStruct$genesInCliques[[k]])[1])))
        {
            for (i in 1:dim(maxCliqueStruct$genesInCliques[[k]])[1])
            {
            maxCliqueStruct$OrderedGenesInCliques[[k]]<-
                maxCliqueStruct$genesInCliques[[k]][unlist(ordWeight[k]),]
            maxCliqueStruct$OrderedIdxInCliques[[k]]<-
                maxCliqueStruct$idxInCliques[[k]][unlist(ordWeight[k]),]
            }
        }
        }
    }
    else
    {
        maxCliqueStruct$OrderedGenesInCliques[[1]]<-
            maxCliqueStruct$genesInCliques[[1]]
        maxCliqueStruct$OrderedIdxInCliques[[1]]<-
            maxCliqueStruct$idxInCliques[[1]]
    }
    
    return(maxCliqueStruct)  
}



###############################################################################
optimizeParamsMuGroup<-function(params,countsVec,n,lo,model)
{
    # (internal) function to be fed into the numerical optimization routine

    # inputs:
        # params: parameter vector
        # countsVec: contingency table of the input genes, as a vector 
        # n: number of genes
        # lo: value of the observation lambda; I set it to 1 for convenience 
        # model: either "Null" or "ME", depending on the model for which 
            # parameters are estimated

    # outputs:
        # s: log likelihood of the sample

    # for each genotype, compute its log likelihood and multiply by the 
    # genotype count
    mm<-sapply(1:length(countsVec),produceLogLikeTerms,params,
            countsVec=countsVec,n=n,lo=lo,model=model)
    # sum them to get the log likelihood of the sample
    s<-sum(mm)
    return(s)
}



###############################################################################
produceLogLikeTerms<-function(idx,params,countsVec,n,lo,model)
{
    # (internal) function to compute the total log probability for a given 
    # genotype with index idx (log probability times count) 

    # inputs:
        # idx: index of the genotype for which the total log probability is to 
            # be computed
        # params: parameter vector
        # countsVec: contingency table of the input genes, as a vector 
        # n: number of genes
        # lo: value of the observation lambda; I set it to 1 for convenience 
        # model: either "Null" or "ME", depending on the model for which 
            # parameters are estimated

    # outputs:
        # term: total log likelihood (log likelihood times count)

    if (model=="ME")
    {
        mu<-params[length(params)]
        lambdas<-params[-length(params)]
    } else if (model=="Null")
    {
        mu<-0
        lambdas<-params
    }

    currentCount<-countsVec[idx]
    obs<-arrayInd(idx,rep(2,n))-1
    if (sum(obs)==0)
        currIndivLike<-compIndivLikeFullZero(lambdas,lo)
    else
        currIndivLike<-compIndivLike(obs,lambdas,lo,mu)
    if (currIndivLike!=0)
        currIndivLogLike<-log(currIndivLike)
    else
        currIndivLogLike<-0
    term<-currIndivLogLike*currentCount
    return(term)
}



###############################################################################
compIndivLikeFullZero<-function(lambdas,lo)
{
    # (internal) function to compute the log likelihood of the null (full 0) 
    # genotype
    
    # inputs:
        # lambdas: values of the lambda parameters for the genes, as vector
        # lo: value of the observation lambda; I set it to 1 for convenience

    # outputs:
        # value: log likelihood of the null (full 0) genotype

    value<-lo/(sum(lambdas)+lo)
    return(value)
}



###############################################################################
compIndivLike<-function(obs,lambdas,lo,mu)
{
    # (internal) function to compute the log likelihood of any genotype 
    # besides the null (full 0) one, for both models

    # inputs:
        # obs: current genotype (as a binary vector)
        # lambdas: values of the lambda parameters for the genes, as vector
        # lo: value of the observation lambda; I set it to 1 for convenience
        # mu: value of the mu parameter

    # outputs:
        # prob: log likelihood of the input genotype
    
    #the indices which are 1
    K<-which(obs==1) 
    #the indices which are 0
    NK<-which(obs==0)
    #the corresponding lambdas to the positions which are 1
    Klambdas<-lambdas[K]
    #the corresponding lambdas to the positions which are 0
    NKlambdas<-lambdas[NK]
    #the sum of the lambdas which correspond to 0 (including lo)
    NKSum<-sum(NKlambdas)+lo
    #the product of the lambdas which correspond to 1
    KProd<-prod(Klambdas)
    #the sum of all lambdas (including lo)
    lamSum<-sum(lambdas)+lo 

    KSize<-length(K) #how many 1 are there

    pk<-KProd/lamSum*lo/NKSum #multiplication term
    
    #all the combinations of lambdas, as in the formula
    options<-permutations(KSize,KSize)
    if (KSize>1)
    {
        #remove half of them, as they are equivalent
        options<-options[-seq(2,(dim(options)[1]),by=2),]
        if (!is.null(nrow(options)))
        terms<-apply(options,1,compProdTerms,lamSum=lamSum,KSize=KSize,
                    Klambdas=Klambdas)
        else 
        terms<-compProdTerms(options,lamSum,KSize,Klambdas) #if only one row
        prob<-(1-mu)*pk*sum(terms) # for the null model, mu=0
    } else #if only one position is 1
    {
        if (mu==0) #null model
        prob<-(1-mu)*pk
        else
        prob<-Klambdas/lamSum*(lo+mu*(lamSum-Klambdas-lo))/(lamSum-Klambdas)
    }

    return(prob)
}



###############################################################################
compProdTerms<-function(perCurr,lamSum,KSize,Klambdas)
{
    # (internal) function to recursively compute the product terms of the log 
    # likelihoods

    # inputs:
        # perCurr: current position
        # lamSum: current sum of lambdas
        # KSize: number of lambdas
        # Klambdas: values of lambdas

    # outputs:
        # prodTerms: product of the individual log likelihoods

    # given a lambdas configuration, the product is of the form 
    # 1/(sum-(lambdas which are 1))
    lamSumUpdate<-lamSum
    prodTerms<-1
    if (KSize>2)
    {
        for (i in 1:(KSize-2))
        {
        prodTerms<-prodTerms*1/(lamSumUpdate-Klambdas[perCurr[i]])
        lamSumUpdate<-lamSumUpdate-Klambdas[perCurr[i]]
        }
    }
    prodTerms<-prodTerms*(1/(lamSumUpdate-Klambdas[perCurr[KSize-1]])+
                            1/(lamSumUpdate-Klambdas[perCurr[KSize]]))
    return(prodTerms)
}



###############################################################################
testAllCandidateGroups<-function(mcStruct,mat)
{
    # (internal) function to test all candidate maximal cliques with TiMEx for 
    # groups

    # inputs:
        # mat: input matrix of binary alterations
        # mcStruct: structure containing maximal cliques, as returned by 
            # function doMaxCliques

    # outputs: 
        # a list consisting of the following lists, each with as many 
        # elements as lengths of the identified cliques:
            # groupTests: each element of this list is a list with as many 
                # elements as maximal cliques of certain length identified. 
                # each element is the result of the TiMEx test on each group, 
                # as returned by function testCliqueAsGroup, including 
                # estimated parameters for the two models, pvalue, likelihoods 
                # etc
            # newMusGroup: each element of this list is a vector of the 
                # estimated mu values for the tested maximal cliques of a 
                # given length
            # pvalueLRTMu: each element of this list is a vector of the 
                # estimated pvalues for the tested maximal cliques of a 
                # given length

    lo<-1 #lambda obs is set to 1 (see the manuscript)
    groupTests<-list()
    newMusGroup<-list()
    pvalueLRTMu<-list()

    for (i in 1:length(mcStruct$detectedLengths))
    {
        if (mcStruct$detectedLengths[i]>1)
        {
        print(paste("clique size =",mcStruct$detectedLengths[i]))
        if (length((dim(mcStruct$Mus$OrderedIdxInCliques
                        [[mcStruct$detectedLengths[i]]])))>0)
        {
        print(paste("number of cliques to test =", 
    dim(mcStruct$Mus$OrderedIdxInCliques[[mcStruct$detectedLengths[i]]])[1]))
        
    groupTests[[mcStruct$detectedLengths[i]]]<-
    apply(mcStruct$Mus$OrderedIdxInCliques[[mcStruct$detectedLengths[i]]],1,
        testCliqueAsGroup,mat=mat,lo=lo)
        newMusGroup[[mcStruct$detectedLengths[i]]]<-
        sapply(groupTests[[mcStruct$detectedLengths[i]]],function(x)
        {x$opMu$par[mcStruct$detectedLengths[i]+1]})
        pvalueLRTMu[[mcStruct$detectedLengths[i]]]<-
        sapply(groupTests[[mcStruct$detectedLengths[i]]],function(x)
        {x$pvalueLRT})
        }
        else
        {
        print(paste("number of cliques to test = 1"))
        groupTests[[mcStruct$detectedLengths[i]]]<-testCliqueAsGroup(
    mcStruct$Mus$OrderedIdxInCliques[[mcStruct$detectedLengths[i]]],mat,lo)
        newMusGroup[[mcStruct$detectedLengths[i]]]<-
groupTests[[mcStruct$detectedLengths[i]]]$opMu$par[mcStruct$detectedLengths[i]
                                            +1]
        pvalueLRTMu[[mcStruct$detectedLengths[i]]]<-
        groupTests[[mcStruct$detectedLengths[i]]]$pvalueLRT
        }
    }
    }
    return(list("groupTests"=groupTests,"newMusGroup"=newMusGroup,
                "pvalueLRTMu"=pvalueLRTMu))
}



###############################################################################
filterSignifCliques<-function(mcStruct,testedCand,groupPvalue)
{
    # (internal) function to filter, among all maximal cliques, the 
    # significant ones (after testing with TiMEx for groups)

    # inputs: 
        # mcStruct: structure containing the maximal cliques, as returned by 
            # function doMaxCliques
        # testedCand: structure containing the results after testing all 
            # candidate maximal cliques with TiMEx for groups, as returned by 
            # function testAllCandidateGroups
        # groupPvalue: threshold for the adjusted pvalue, lower than which 
            # cliques are significant

    # outputs: list consisting of the following lists, with as many elements 
        # as lengths of the identified maximal cliques:
            # genesSignif: names of genes part of the significant groups (for 
                # both fdr and bonferroni)
            # posSignif: positions in the input list of the groups which are 
                # significant (for both fdr and bonferroni)
            # idxSignif: indices (in the input matrix) of genes part of the 
                # significant groups (for both fdr and bonferroni)
        # MusGroup: mu values of the significant groups after testing with 
            # TiMEx for groups (for both fdr and bonferroni)
        # pvals: adjusted p-values (fdr and bonferroni) of the significant 
            # groups

    posSignif<-list()
    genesSignif<-list() 
    idxSignif<-list()
    MusGroup<-list() 
    pvals<-list() 

    for (i in 1:length(mcStruct$detectedLengths))
        {
        if (mcStruct$detectedLengths[i]>1)
        {
            posSignif[[mcStruct$detectedLengths[i]]]<-list()
            MusGroup[[mcStruct$detectedLengths[i]]]<-list()
            pvals[[mcStruct$detectedLengths[i]]]<-list()
        
            genesSignif[[mcStruct$detectedLengths[i]]]<-list()
            idxSignif[[mcStruct$detectedLengths[i]]]<-list()
            
posSignif[[mcStruct$detectedLengths[i]]]$fdr<-which(p.adjust(
    testedCand$pvalueLRTMu[[mcStruct$detectedLengths[i]]],method="fdr")<=
    groupPvalue)
posSignif[[mcStruct$detectedLengths[i]]]$bonf<-which(p.adjust(
    testedCand$pvalueLRTMu[[mcStruct$detectedLengths[i]]],method="bonferroni")
    <=groupPvalue)
        
pvals[[mcStruct$detectedLengths[i]]]$fdr<-(p.adjust(
    testedCand$pvalueLRTMu[[mcStruct$detectedLengths[i]]],method=
    "fdr"))[posSignif[[mcStruct$detectedLengths[i]]]$fdr]
pvals[[mcStruct$detectedLengths[i]]]$bonf<-(p.adjust(
    testedCand$pvalueLRTMu[[mcStruct$detectedLengths[i]]],method=
    "bonferroni"))[posSignif[[mcStruct$detectedLengths[i]]]$bonf]
        
MusGroup[[mcStruct$detectedLengths[i]]]$fdr<-testedCand$
    newMusGroup[[mcStruct$detectedLengths[i]]][posSignif[[mcStruct$
    detectedLengths[i]]]$fdr]
MusGroup[[mcStruct$detectedLengths[i]]]$bonf<-testedCand$
    newMusGroup[[mcStruct$detectedLengths[i]]][posSignif[[mcStruct$
    detectedLengths[i]]]$bonf]
        
            # in case only one group was detected
            if (length(dim(mcStruct$Mus$OrderedGenesInCliques[[mcStruct$
                                        detectedLengths[i]]]))==0)
            {
genesSignif[[mcStruct$detectedLengths[i]]]$fdr<-t(matrix(mcStruct$Mus$
    OrderedGenesInCliques[[mcStruct$detectedLengths[i]]]))[posSignif[[mcStruct$
    detectedLengths[i]]]$fdr,]
genesSignif[[mcStruct$detectedLengths[i]]]$bonf<-t(matrix(mcStruct$Mus$
    OrderedGenesInCliques[[mcStruct$detectedLengths[i]]]))[posSignif[[mcStruct$
    detectedLengths[i]]]$bonf,]

idxSignif[[mcStruct$detectedLengths[i]]]$fdr<-t(matrix(mcStruct$Mus$
    OrderedIdxInCliques[[mcStruct$detectedLengths[i]]]))[posSignif[[mcStruct$
    detectedLengths[i]]]$fdr,]
idxSignif[[mcStruct$detectedLengths[i]]]$bonf<-t(matrix(mcStruct$Mus$
    OrderedIdxInCliques[[mcStruct$detectedLengths[i]]]))[posSignif[[mcStruct$
    detectedLengths[i]]]$bonf,]

        } else{

    genesSignif[[mcStruct$detectedLengths[i]]]$fdr<-mcStruct$Mus$
    OrderedGenesInCliques[[mcStruct$detectedLengths[i]]][posSignif[[mcStruct$
    detectedLengths[i]]]$fdr,]
genesSignif[[mcStruct$detectedLengths[i]]]$bonf<-mcStruct$Mus$
    OrderedGenesInCliques[[mcStruct$detectedLengths[i]]][posSignif[[mcStruct$
    detectedLengths[i]]]$bonf,]

idxSignif[[mcStruct$detectedLengths[i]]]$fdr<-mcStruct$Mus$
    OrderedIdxInCliques[[mcStruct$detectedLengths[i]]][posSignif[[mcStruct$
    detectedLengths[i]]]$fdr,]
idxSignif[[mcStruct$detectedLengths[i]]]$bonf<-mcStruct$Mus$
    OrderedIdxInCliques[[mcStruct$detectedLengths[i]]][posSignif[[mcStruct$
    detectedLengths[i]]]$bonf,]


        # order all the lists by the corrected pvalue
        if (length(order(pvals[[mcStruct$detectedLengths[i]]]$fdr))>1)
        {
posSignif[[mcStruct$detectedLengths[i]]]$fdr<-posSignif[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)]
genesSignif[[mcStruct$detectedLengths[i]]]$fdr<-genesSignif[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr),]
idxSignif[[mcStruct$detectedLengths[i]]]$fdr<-idxSignif[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr),]
MusGroup[[mcStruct$detectedLengths[i]]]$fdr<-MusGroup[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)]
pvals[[mcStruct$detectedLengths[i]]]$fdr<-pvals[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)] 
        }
    else 
    {
posSignif[[mcStruct$detectedLengths[i]]]$fdr<-posSignif[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)]
genesSignif[[mcStruct$detectedLengths[i]]]$fdr<-genesSignif[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)]
idxSignif[[mcStruct$detectedLengths[i]]]$fdr<-idxSignif[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)]
MusGroup[[mcStruct$detectedLengths[i]]]$fdr<-MusGroup[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)]
pvals[[mcStruct$detectedLengths[i]]]$fdr<-pvals[[mcStruct$
    detectedLengths[i]]]$fdr[order(pvals[[mcStruct$detectedLengths[i]]]$fdr)] 
    }

    if (length(order(pvals[[mcStruct$detectedLengths[i]]]$bonf))>1)
        {
posSignif[[mcStruct$detectedLengths[i]]]$bonf<-posSignif[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)]
genesSignif[[mcStruct$detectedLengths[i]]]$bonf<-genesSignif[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf),]
idxSignif[[mcStruct$detectedLengths[i]]]$bonf<-idxSignif[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf),]
MusGroup[[mcStruct$detectedLengths[i]]]$bonf<-MusGroup[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)]
pvals[[mcStruct$detectedLengths[i]]]$bonf<-pvals[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)] 
    }
    else 
    {
posSignif[[mcStruct$detectedLengths[i]]]$bonf<-posSignif[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)]
genesSignif[[mcStruct$detectedLengths[i]]]$bonf<-genesSignif[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)]
idxSignif[[mcStruct$detectedLengths[i]]]$bonf<-idxSignif[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)]
MusGroup[[mcStruct$detectedLengths[i]]]$bonf<-MusGroup[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)]
pvals[[mcStruct$detectedLengths[i]]]$bonf<-pvals[[mcStruct$
    detectedLengths[i]]]$bonf[order(pvals[[mcStruct$detectedLengths[i]]]$bonf)] 
    }
            }
        }
    }
    return(signifCliques=list("genesSignif"=genesSignif,"idxSignif"=idxSignif,
                    "pvals"=pvals,"posSignif"=posSignif,"MusGroup"=MusGroup))
}