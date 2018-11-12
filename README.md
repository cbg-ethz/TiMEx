# TiMEx

The TiMEx R package is the implementation of a generative probabilistic graphical model for the *de novo* identification of patterns of various degrees of mutual exclusivity across genetic alterations, which can indicate pathways involved in cancer progression. For more information on the underlying model or on biological applications, see the publication introducing TiMEx: 

*Constantinescu, Simona, et al. "TiMEx: a waiting time model for mutually exclusive cancer alterations." Bioinformatics 32.7 (2015): 968-975*.

## Main idea of the model
We regard tumorigenesis as a dynamic process, and base our model on the temporal interplay between the waiting times to alterations, characteristic for every gene and alteration type (mutation or copy number alteration), and the observation time. Under the assumption of rarity of events over short time intervals, TiMEx models the alteration process for each gene as a Poisson process. The waiting times are therefore modeled as exponentially distributed variables with specific rates, which correspond to the rates of evolution for each alteration. In our modeling framework, the temporal dynamics of each alteration process progresses from the onset of cancer, corresponding to the first alteration related to the growth of a malignant tumor, until the observation time, corresponding to the time of tumor biopsy. The observation time is regarded as a system failure time, and is exponentially distributed with an unknown rate. 

## Computational procedure
The computational procedure implemented in the TiMEx R package is aimed to efficiently identify mutually exclusive groups of alterations of any size in large cancer datasets. The procedure includes three steps:

1. Testing all pairs of alterations and identifying mutually exclusive pairs
2. Finding cliques in the graph defined by the previously identified pairwise interactions
3. Testing all candidate groups identified in step 2, and ranking them by their significance or intensity of mutual exclusivity.

## What is unique about TiMEx
TiMEx is the first method that describes mutual exclusivity as a consequence of a dynamic process in time. Unlike previous *de novo* approaches, TiMEx infers functional relations among genes based on an underlying temporal representation of the process of gene alteration in tumorigenesis. 

Importantly, TiMEx infers the mutual exclusivity intensity of a group as a continuous measure, which can be interpreted as a probability. This feature is biologically justified, since the small, but observable, increase in tumor fitness due to multiple alterations in a group of functionally related genes supports the hypothesis that mutual exclusivity occurs at various degrees, as opposed to a binary classification. 

Unlike most other approaches, TiMEx does not explicitly impose frequency constraints, and detects both high frequent and very low frequent alterations, solely based on the temporal relations between them. It identifies all mutually exclusive gene groups of various, not pre-defined sizes, and performs highly efficiently on large datasets.

## Installation
TiMEx has been tested on Mac and Linux. 

Installing TiMEx requires the R package ```devtools```, and can be easily installed in R as ```devtools::install_github('csimona/TiMEx')``` 

It depends on the R packages ```RBGL```, ```graph``` and ```gtools```, which, if missing, will be installed together with TiMEx, so no extra input should be required from the user. 

## Input
The input to the TiMEx analysis is a binary alteration matrix, with rows representing patients and columns representing alterations. The TiMEx R package already contains preprocessed input datasets that can be used as examples. These datasets can be loaded as, for example, ```data(ovarian)```. 

Additional datasets include ```breast``` (TCGA breast cancer), ```gbmDendrix``` (GBM dataset preprocessed as described in the paper introducing the method *Dendrix*), or ```gbmMuex``` (GBM dataset preprocessed as described in the paper introducing the method *muex*).

## Output
TiMEx returns a list consisting of:

* ```genesSignif```: list of significantly mutually exclusive groups, as gene names, sorted by corrected p-value. The list contains as many elements as identified lengths of groups. For example, ```genesSignif[[2]]``` is a list containing the gene names of the significant groups of size 2. Each list of this type further has two elements, ```fdr``` and ```bonf```, corresponding to different multiple testing correction methods. Each element is a matrix, in which rows represent gene names of significantly mutually exclusive groups.

* ```idxSignif```: list of significantly mutually exclusive groups, as indices in the input matrix, sorted by corrected p-value. The list contains as many elements as identified lengths of groups. For example, ```idxSignif[[2]]``` is a list containing the indices of the significant groups of size 2. Each list of this type further has two elements, ```fdr``` and ```bonf```, corresponding to different multiple testing correction methods. Each element is a matrix, in which rows represent indices of significantly mutually exclusive groups.

* ```pvals```: list of corrected significant p-values corresponding to the tested cliques, ordered ascendingly. The list contains as many elements as identified lengths of significant groups. For example, ```pvals[[2]]``` is a list containing the p-values of the significant maximal cliques of size 2. Each list of this type further has two elements, ```fdr``` and ```bonf```, corresponding to different multiple testing correction methods. Each element is a vector, of length the number of significant maximal cliques of a given size.      
          
* ```posSignif```: list of positions of the significant groups in the input list of maximal cliques, ordered ascendingly by corrected p-value.  The list contains as many elements as identified lengths of significant groups. For example, ```posSignif[[2]]``` is a list containing the positions of the significant groups of size 2.  Each list of this type further has two elements, ```fdr``` and ```bonf```, corresponding to different multiple correction methods.  Each element is a vector, of length the number of significant maximal cliques of a given size.
          
* ```MusGroup```: list of inferred mu values corresponding to the tested cliques, ordered ascendingly by the corresponding corrected p-value. The list contains as many elements as identified lengths of significant groups. For example, ```MusGroup[[2]]``` is a list containing the mu values of the significant maximal cliques of size 2. Each list of this type further has two elements, ```fdr``` and ```bonf```, corresponding to different multiple testing correction methods. Each is a vector, of length the number of significant maximal cliques of a given size.
          
* ```mcStruct```: input structure of maximal cliques to be tested for mutual exclusivity, as returned by ```doMaxCliques```.

* ```matrix```: input binary alteration matrix.

* ```groupPvalue```: input threshold for the corrected p-value, lower than which cliques are significant.


## Example
The main function of this R package is ```TiMEx```. An example run is:

```R
data(ovarian)
ovarianMEGroups <- TiMEx(ovarian)
```

TiMEx will run with default parameters. For additional details on the parameters that the user can set, as well as on the method or the input or output structures, start by typing 

```R
?TiMEx
``` 

at the R console.

Before running TiMEx, it's a good idea to run the function ```doMetagroup```  on your binary matrix, and run TiMEx on the ```$mat``` part of the result. ```doMetagroup``` condenses all genes with identical alteration patterns into one metagene. These genes are equivalent from the point of view of the algorithm, so not condesing them into a single metagene would produce ambigous results.


## Tips
Please remove genes that are all-0 across patients from your data, they are not informative in identifying mutually exclusive patterns.

In cases when the input consists of a very low sample size, errors might occur (estimation fails) when the data contains pairs of genes that have no [00] entries, meaning they are not co-occurring wild type in any patient. If TiMEx fails on your data, please check this first. 

Each function has a detailed and lengthy explanation, together with working examples, which provide an easy-to-follow hans-on guide on how to use TiMEx. Please report any issues you find here; sharing your feedback is greatly appreciated. 