library(diffHic)
library(Biostrings)
library(BSgenome)
library(dplyr)
library(ggplot2)




fls=list.files(path="/uufs/chpc.utah.edu/common/home/holmen-group1/star_ms/", pattern = "*m39.dna.primary_assembly.fa$", full.names = T)




#cutsM39_cctagg <- cutGenome(fls, "cctagg", overhang=4) #GATC AspA2I AvrII BlnI XmaJI
#cutsM39_cttaag <- cutGenome(fls, "cttaag", overhang=4) #AATT AFIII
cutsM39_AGATCT <- cutGenome(fls, "AGATCT", overhang=4) #CTAG BgIII Bsp1407I BsrGI BstAUI SspBI

#cutsM39_cctagg_tibble<-as_tibble(cutsM39_cctagg)
cutsM39_AGATCT_tibble<-as_tibble(cutsM39_AGATCT)
#cutsM39_cctagg_tibble<-as_tibble(cutsM39_cctagg)







totaln<-cutsM39_AGATCT_tibble%>%
  summarise(wtotal=sum(width))



under_x<-30000
over_x<-50

cutsM39_AGATCT_tibble%>%
  mutate(under_x=if_else(width<=under_x, "un", "ov"))%>%
  group_by(under_x)%>%
  summarise(totalw=sum(width), 
            perc=totalw/totaln)




# ggplot(cutsM39_cctagg_tibble, aes(y=width))+
#   geom_boxplot(outlier.size=2,outlier.colour="green")+
#   scale_y_continuous(trans='log2')+
#   ggtitle("Mouse Genome Cuts by AspA2I, AvrII, BlnI, XmaJI (c*ctagg)")
logbase<-2


gdata<-cutsM39_AGATCT_tibble%>%
  rename(fragl=width)%>%
  
  #filter(fragl<=under_x)%>%
  
  mutate(logfragl=log(fragl, base = logbase))


ladder<- c(0, 100, 200, 300, 
           400, 500, 650, 850, 
           1500,
           2000,
           3000,
           4000,
           5000, 6000, 7000, 8000, 10000, 15000)


logladder<-log(ladder, base=logbase)


# ggplot(gdata, aes(x=fragl))+
#   
#   
#   geom_freqpoly(binwidth=10)+
#   
#   scale_x_continuous(breaks = ladder)+
#   coord_flip()
# 
# 





library(ggplot2)
library(dplyr)

top_cutoff<-16000
bottom_cutoff<-75

# Parameters
start <- bottom_cutoff-1        # starting value (e.g., 1)
stop <- top_cutoff+1       # maximum edge
base <- 1.001         # exponential growth base (e.g., 2x each step)

# Generate exponentially spaced values
bin_edges <- unique(c(
  0,                         # always include zero if desired
  start * base^(0:10000) # change 20 to control how many bins
))
# Filter to not exceed the stop value
bin_edges <- bin_edges[bin_edges <= stop]

# Optionally add the stop as a final bin
bin_edges <- c(bin_edges, stop)

bin_edges<-unique(c(bin_edges, top_cutoff+1))




# Create bin midpoints
bin_midpoints <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2

# Assign bins
gdata_binned <- gdata %>%
  mutate(fragl=if_else(fragl>top_cutoff, top_cutoff, fragl))%>%
 filter(fragl>=bottom_cutoff)%>%
  
  mutate(bin = cut(fragl, breaks = bin_edges, include.lowest = TRUE)) %>%
  group_by(bin) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(midpoint = bin_midpoints[as.integer(bin)])

# # Plot
# ggplot(gdata_binned, aes(x = midpoint, y = count)) +
#   geom_line() +
#   # geom_point() 
#   scale_x_continuous(breaks = ladder)+
#   coord_flip()





ladder<- c(0, 100, 200, 300, 
           400, 500, 650, 850, 
           1500,
           2000,
           3000,
           4000,
           5000, 6000, 7000, 8000, 10000, 15000)






library(ggplot2)
library(scales)

# Plot with log-transformed x-axis
ggplot(gdata_binned, aes(x = midpoint, y = count)) +
  geom_line() +
  scale_x_continuous(
    trans = 'log10',
    breaks = ladder,
    labels = comma_format()
  ) +
  coord_flip()

