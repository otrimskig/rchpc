source("libs.R")
library(tidyverse)
library(dtplyr)
library(EnhancedVolcano)


df<-readRDS("k19mf/ds/immune2-acral-v-subq-stats.rds")


voldf<-df%>%filter(tumor_type=="acral")%>%
  column_to_rownames("gene_name_ms")%>%
  filter(is.finite(log2_fold_change))



EnhancedVolcano(voldf,
                lab = rownames(voldf),
                pCutoff = 10e-3,
                
                x = 'log2_fold_change',
                y = 'p_unadjusted',
               
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                
                
                
                col = c('#757575', "#6188ff", 'purple', 'red3'))
