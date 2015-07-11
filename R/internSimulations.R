# Author: Simona Constantinescu; simona.constantinescu@bsse.ethz.ch

###############################################################################
compAllProbsSims<-function(obs,lambdas,lo,mu)
{
    # (internal) wrapper function to compute the log likelihood of any 
    # genotype besides the null (full 0) one, for given lambdas
    
    # inputs:
        # obs: current genotype (as a binary vector)
        # lambdas: values of the lambda parameters for the genes, as vector
        # lo: value of the observation lambda; I set it to 1 for convenience
        # mu: value of the mu parameter
    
    # outputs:
        # value: log likelihood of the input genotype
    
    value<-compIndivLike(obs,lambdas,lo,mu)
    return(value)
}



###############################################################################
compProbFullZeroSims<-function(obs,lambdas,lo,mu)
{
    # (internal) wrapper function to compute the log likelihood of the null 
    # (full 0) genotype, for given lambdas

    # inputs:
        # obs: current genotype (as a binary vector)
        # lambdas: values of the lambda parameters for the genes, as vector
        # lo: value of the observation lambda; I set it to 1 for convenience
        # mu: value of the mu parameter
    
    # outputs:
        # value: log likelihood of the input genotype
    
    value<-compIndivLikeFullZero(lambdas,lo)
    return(value)
}



###############################################################################
produceGenesGroup<-function(NCur,genoMat,n)
{
    #' (internal) main function to simulate genes corresponding to the 
    #' probability distribution in given as input genoMat (n=number of genes; 
    #' NCur=sample size)
    
    # inputs:
        # NCur: number of observations to simulate
        # genoMat: matrix with genotype probabilities, from which the dataset 
                    # was simulated. This matrix has as many dimensions as 
                    # number of genes, i.e. 2x2x...x2. For each dimension, the 
                    # first position corresponds to the probability of 
                    # observing a 0 on that gene, and the second position 
                    # corresponds to the probability of observing an 1. 
                    # For example, in the case of 4 genes, the probability of 
                    # observing (0000) is given by genoMat[1,1,1,1]; the 
                    # probability of observing (1011) is given by 
                    # genoMat[2,1,2,2]; the probability of observing (1111) is 
                    # given by genoMat[2,2,2,2]. The entries of this matrix
                    # should be nonnegative and sum up to 1.
        # n: number of genes
    
    # outputs: a list consisting of:
        # genes: binary matrix containing the simulated genes according to the 
                # probability distribution in genoMat
        # indivLikes: vector of likelihoods (one for each observation), of 
                # length equals number of observations
        # indivLogLikes: vector of log likelihoods (one for each observation), 
                # of length equals number of observations
        # counts: matrix of same dimensions as genoMat. each position 
                # represents the counts of each genotype in the simulated 
                # dataset. For more information on the format of this matrix,
                # see the input parameter genoMat.
    
    p<-produceGenesLikesGroup(NCur,genoMat,n)
    likes<-unlist(p[,2])
    logLikes<-log(likes)
    genes<-matrix(unlist(p[,1]),ncol=n,byrow=TRUE)
    counts<-Reduce('+',p[,3])
    return(l=list("genes"=genes,"indivLikes"=likes,"indivLogLikes"=logLikes,
                "counts"=counts))
}



##############################################################################
produceGenesLikesGroup<-function(NCur,genoMat,n)
{
    # (internal) wrapper function to generate NCur observations, for given 
    # distribution of genotypes
    
    # inputs:
        # NCur: number of observations to simulate
        # genoMat: matrix with genotype probabilities, from which the dataset 
                # was simulated. This matrix has as many dimensions as 
                # number of genes, i.e. 2x2x...x2. For each dimension, the 
                # first position corresponds to the probability of 
                # observing a 0 on that gene, and the second position 
                # corresponds to the probability of observing an 1. 
                # For example, in the case of 4 genes, the probability of 
                # observing (0000) is given by genoMat[1,1,1,1]; the 
                # probability of observing (1011) is given by 
                # genoMat[2,1,2,2]; the probability of observing (1111) is 
                # given by genoMat[2,2,2,2]. The entries of this matrix
                # should be nonnegative and sum up to 1.
    # n: number of genes
    
    # outputs:
        # obs: matrix with as many rows as observations; each row consists of 
            # the observation, its probability, and a matrix of the same 
            # dimensions as genoMat, which has an 1 for the genotype which was 
            # sampled, and the rest of the entries 0.
    
    obs<-t(sapply(1:NCur,produceSingleObsGroup,genoMat=genoMat,n=n))
    return(obs)
}



##############################################################################
produceSingleObsGroup<-function(i,genoMat,n)
{
    # (internal) function to generate a single observation, for given 
    # distribution of genotypes 
    
    # inputs:
        # i: number of observations to simulate
        # genoMat: matrix with genotype probabilities, from which the dataset 
                # was simulated. This matrix has as many dimensions as 
                # number of genes, i.e. 2x2x...x2. For each dimension, the 
                # first position corresponds to the probability of 
                # observing a 0 on that gene, and the second position 
                # corresponds to the probability of observing an 1. 
                # For example, in the case of 4 genes, the probability of 
                # observing (0000) is given by genoMat[1,1,1,1]; the 
                # probability of observing (1011) is given by 
                # genoMat[2,1,2,2]; the probability of observing (1111) is 
                # given by genoMat[2,2,2,2]. The entries of this matrix
                # should be nonnegative and sum up to 1.
    # n: number of genes
    
    # outputs:
        # l: list with the following elements:
            # newObs: binary vector representing the sampled observation
            # newProb: probability of the sampled observation
            # newSample: matrix of the same dimensions as genoMat, which has 
                        # an 1 for the genotype which was sampled, and the 
                        # rest of the entries 0.
    
    # has an 1 for the observation whose probability was sampled
    newSample<-array(rmultinom(1,1,genoMat),dim=rep(2,n))
    newIdx<-which(newSample==1,arr.ind=TRUE)
    newProb<-genoMat[newIdx]
    newObs<-newIdx-1
    return(l=list("newObs"=newObs,"newProb"=newProb,"newSample"=newSample))
}