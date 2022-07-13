## code to prepare `DATASET` dataset goes here

library(utils)
library(fgsea)

dfGOBP = read.delim("data-raw/2022_January3_GOID_GOBP_SGD.txt",stringsAsFactors = F,check.names = FALSE)

yGOBP.gmt   = gmtPathways(gmt.file = "data-raw/2022_January_3_GO_BP.gmt")

hGOBP.gmt = gmtPathways(gmt.file = "data-raw/STable3_hGOBP_no_iea_July012017_gene.gmt")

sampleFitdata = as.matrix(read.delim("data-raw/sampleFitdata.txt",stringsAsFactors = F,check.names = FALSE))


usethis::use_data(dfGOBP, yGOBP.gmt, hGOBP.gmt, sampleFitdata,overwrite = T)









usethis::use_data(DATASET, overwrite = TRUE)
