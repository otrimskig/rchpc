setwd("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig")
library(tidyverse)
# https://github.com/katharineyu/TCGA_CCLE_paper/blob/master/README.md

install.packages("synapser", repos = c("http://ran.synapse.org", "http://cran.fhcrc.org"))

install.packages("synapserutils")
library(synapser) 
library(synapserutils) 

synLogin(authToken="") 
files <- synapserutils::syncFromSynapse('syn18685536') 