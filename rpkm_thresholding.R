library(tidyverse)
mutate(max_rpkm = sum(rpkm))%>%
  filter(max_rpkm>=1)