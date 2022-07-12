## code to prepare `DATASET` dataset goes here

library(utils)

dfGOBP = read.delim("data-raw/2021_Decenber30_GOID_GOBP_SGD.txt",stringsAsFactors = F,check.names = FALSE)

GOBPgmt   = readRDS("data-raw/2021_Decenber30_GO_BP.RDS")


sampleFitdata = read.delim("data-raw/sampleFitdata.txt",stringsAsFactors = F,check.names = FALSE)


usethis::use_data(GOBPgmt,dfGOBP,sampleFitdata,overwrite = T)









usethis::use_data(DATASET, overwrite = TRUE)
