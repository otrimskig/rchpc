library(tidyverse)
library(dtplyr)
library(GSVA)


#read in your data. 
#to perform analysis, you will need your data in the form of
#a matrix of RPKMs, with rownames corresponding to the gene name/code/symbol, in the same
#format as in the gene signatures list, and the same species. 

#in this case our signatures list is human gene symbols, so that's what we use. 

mat<-readRDS("rebecca/rpkms_example.rds")


#this file contains several lists of pre-made gene signature lists.
load("rebecca/allSignatures.rda")

#use this to add additional gene signatures lists if desired.
#you can get additional lists as .gmt files.
#all they are is a list of signatures and the genes they contain.
#you can make your own or get a series of them from 
#https://www.gsea-msigdb.org/gsea/index.jsp or other source.
#gsea-mgsigdb has many signatures lists for both human and mouse,
#and for different gene id types (ensembl id, or gene symbol, for example.)

#additional_signatures<-qusage::read.gmt("x.gmt")


#concatenate all your desired signature lists together.
signature_lists <- c(Hallmark, kegg, TIMEx, Immune_sig)


#now perform the analysis.
#make parameter list.
ssgeaP<-ssgseaParam(mat, signature_lists)
#call GSVA on the parameter list.
#use supressWarnings to continue analysis if there are 
#any gene lists with only 1 gene. 
ssgsvaG<- suppressWarnings(gsva(ssgeaP))

#get unscaled (non-normalized) matrix of scores for each pathway.
#scores are calculated per sample. So any sample exclusions or 
#group comparisons are not relevant at this point.
gsva_u <- t(t(ssgsvaG))


saveRDS(gsva_u, "rebecca/gsva_u.rds")