
# GOenrichment

<!-- badges: start -->
<!-- badges: end -->

GOenrichments is a package that facilitates GO set enrichment
analysis and subsequent network visualization. The significance 
of the GO enrichments are estimated using the hypergeometric 
test based on the query set of genes passing the user-defined
fitness score threshold compared to the genes in the background 
set (the gene universe). Because of the high degree of redundancy
between GO terms, an advantageous of the GOenrichment package
is to simplify the output of a GO enrichment analysis. Specifically,
clusters in the enrichment network are defned by GO terms that 
overlap by >= 50%. Node size is proportional to the significane 
of the GO enrichment (FDR score) and edge weight are proportional 
to the degree of overlap between the GO terms. In each cluster
we highlight the most significantly enriched GO term by increasing
the fontsize and setting the fontface to bold. In this manner,
cellular and functional modules emerge in the network, defined by tightly
connected color coordinated clusters connected by multiple edges and
a holistic view of the cellular response to the perturbation ofo 
interest is possible. 

## Installation

You can install the development version of GOenrichment like so:


devtools::install_github("gurinina/GOenrichment")


## Example

This is a basic example which shows you how to solve a common problem:


library(GOenrichment)
# this is a GO enrichment for a yeast fitness dataset, run in presence of a drug called azithromycin

# runGORESP produces the enrichments
goresp = runGORESP(mat=sampleFitdata,coln=1,curr_exp = colnames(sampleFitdata)[1],sig = 1,fdrThresh = 0.2,bp_path = NULL,bp_input = yGOBP.gmt,go_path = NULL,go_input = dfGOBP)

# visSetup prepares the output of runGORESP for downstream plotting
vis = visSetup(goresp$enrichInfo,goresp$edgeMat```

# runNetwork produces the GO enrichment network
runNetwork(vis$nodes,vis$edges)

