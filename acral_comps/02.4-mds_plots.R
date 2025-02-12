source("libs.R")

library(tidyverse)
library(ggrepel)
library(ggnewscale)


#get all sample metadata
all_info<-readRDS("acral_comps/ds/v00-sample_info.rds")%>%
  mutate(exp=if_else(is.na(exp)&as.numeric(mouse_num)>=38533, "paired_ac_sub", exp))


#get coordinates from plot-mds analysis, and dimension length.

mds_files<-list.files("acral_comps/ds", pattern = "mds_dge[0-9]", full.names = T)


for (i in 1:length(mds_files)){


file_base<-sub(".rds", "", basename(mds_files[i]))



mds<-readRDS(mds_files[i])

mds_coords<-tibble(x=mds$x,
       y=mds$y,
       sample_id=attributes(mds[["distance.matrix.squared"]])[["dimnames"]][[1]])
x_dim<-round(mds[["var.explained"]][1]*100, 1)
y_dim<-round(mds[["var.explained"]][2]*100, 1)

dim_ratio<-round(y_dim/x_dim, 2)

#join mds coordinates to sample info.
dfp<-mds_coords%>%left_join(all_info)%>%
  filter(!is.na(sample_type))



plot_title<-"Experiment and Injection Type"
category1<-"exp"
category2<-"sample_type"
plot_to_save<-ggplot(dfp, aes(x, y))+

    ggtitle(plot_title)+
    #labs(subtitle = plot_subtitle)+



    #layer1
    #sets bubble inside colors
    geom_point(aes(color=!!sym(category2)), shape=21, size=6, stroke=5, alpha=.6)+
    

    ####
    #scale_color_manual(values=colors[[category]])+

    #set up legend prior to resetting color scale.



    guides(color = guide_legend(override.aes = list(size = 7),
                                title=""))+

    new_scale_color()+

    #layer2
    #sets bubble outline colors
    geom_point(aes(color=!!sym(category1)), size=9, shape=20, alpha=.6)+
    #sets scale based on color mapping df, and darkens.



    ####
    #scale_color_manual(values=sapply(colors[[category]],function(x) colorspace::darken(x, 0.2)))+

    #set up legend, prior to resetting color scale. Must match above legend otherwise
    #legend will split.
    guides(color = guide_legend(title=""))+

    #new_scale_color()+

    #set size scale for bubbles
    #scale_size(range = c(2,10))+





    #add text labels for points, by mouse number.
    geom_text_repel(aes(label = mouse_num),
                    min.segment.length = 0,
                    segment.color = "grey80",
                    force=30,
                    max.overlaps = 30,
                    point.padding = 1,
                    show.legend = FALSE)+
    #set color of text to be slightly darker than that of bubble.



    ########
  #scale_color_manual(values=sapply(colors[[category]],function(x) colorspace::darken(x, 0.3)))+

  #remove legend for size
  # guides(size="none")+

    #change y and x axes labels to include percent difference in MDS dimension.
    labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
    labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+

    #set aspect ratio - (probably either 1:1 or scale to %diff).
    coord_fixed(ratio=set_dim_ratio)+
    theme_classic()

  # plot_to_save


  object<-plot_to_save

  ggsave(paste0("acral_comps/plots/",fs::path_sanitize(paste0("mds-", file_base, category1, category2, "unscaled", ".", set_dim_ratio, "ratio", "-m.pdf"))),
         plot=object,

         title=paste0("src: ",
                      
                      rstudioapi::getSourceEditorContext()$path%>%
                        sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/","",.)%>%
                        sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                      
                      " at ", 
                      
                      lubridate::round_date(Sys.time(), "second")
         ),
         
         
         limitsize = FALSE,
         
         
         height=8,
         width=12,
         scale = 1.7,
         dpi=600,
         
         
         
  )
  
  
}
