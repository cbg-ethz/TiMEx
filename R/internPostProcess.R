# Author: Simona Cristea; scristea@@jimmy.harvard.edu

##############################################################################
getFreqsForCliquesAfterTest<-function(cliqueIdxList,matrix)
{
    # (internal) function to compute the frequencies of the genes part of a 
    # group
    
    # inputs:
        # cliqueIdxList: a list of indices of significant groups, as part of 
            # the signifGroups structure, returned by TiMEx. 'fdr' and 'bonf' 
            # represent the different correction methods used for the 
            # testing of groups
        # matrix: binary input matrix of alterations
    
    # outputs: 
        # freqs: a list with two elements, 'fdr' and 'bonf' (corresponding to 
                # the respective multiple correction method), 
                # each consisting of the vector frequencies of the genes in 
                # the group

    l<-length(cliqueIdxList)
    freqs<-list()
    for (i in 2:l)
    {
        freqs[[i]]<-list()
    
        cliqueIdx.bonf<-(cliqueIdxList[[i]])$bonf
        if (length(dim(cliqueIdx.bonf)[1])>0)
        freqs[[i]]$bonf<-apply(cliqueIdx.bonf,2,function(x)
            {apply(matrix[,x],2,sum)/dim(matrix)[1]*100}) else
            freqs[[i]]$bonf<-sum(matrix[,cliqueIdx.bonf])/(dim(matrix)[1]*100)
    
    
        cliqueIdx.fdr<-(cliqueIdxList[[i]])$fdr
        if (length(dim(cliqueIdx.fdr)[1])>0)
        freqs[[i]]$fdr<-apply(cliqueIdx.fdr,2,function(x)
            {apply(matrix[,x],2,sum)/dim(matrix)[1]*100}) else
            freqs[[i]]$fdr<-sum(matrix[,cliqueIdx.fdr])/(dim(matrix)[1]*100)
    }
    return(freqs)
}



##############################################################################
combineNameWithFreq<-function(genesSignifMat,ffMat,muMat,pvalMat)
{
    # (internal) function to adjacently print the name of the genes in a 
    # group, their frequency, the estimated mu, and the pvalue of the group 
    
    # inputs:
        # genesSignifMat: matrix in which each row consists of the gene names 
            # (as characters) of one identified significant maximal clique of 
            # a given size. there are as many rows as identified maximal 
            # cliques
        # ffMat: matrix with same structure as genesSignifMat, consisting 
                # however of positive real numbers between 0 and 100, 
                # representing the frequency in the input dataset of each gene 
                # present in genesSignifMat
        # muMat: vector with as many elements as identified maximal cliques, 
                # consisting of the estimated mu values for each maximal 
                # clique (real values between 0 and 1)
        # pvalMat: vector with as many elements as identified maximal cliques, 
                # consisting of the adjusted pvalue for each maximal clique 
                # (real values between 0 and 1)
    
    # outputs:
        # newMat: matrix with as many rows as identified maximal clique. 
                # each row consists of the names of the genes part of each 
                # clique (followed, in brackets, by their frequency in the 
                # input dataset), the intensity of mutual exclusivity of the 
                # group, and the corrected p-value of the group
    
    d<-dim(genesSignifMat)[1]
    if (!is.null(d) && !d==0)
    {
        newMat<-array(NA,dim=c(d,(dim(genesSignifMat)[2])+2)) 
        for (j in 1:d)
        {
            newMat[j,]<-c(paste(genesSignifMat[j,]," (",
                    round(ffMat[j,],2),"%)",sep=""),
                    round(muMat[j],2),round(pvalMat[j],6))
        }
    } else
    {
        newMat<-rep(NA,(length(genesSignifMat)+2))
        newMat<-c(paste(genesSignifMat," (",
                    round(ffMat,2),"%)",sep=""),
                    round(muMat,2),round(pvalMat,6))
    }
    
    return(newMat)
}