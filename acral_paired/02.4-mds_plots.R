source("libs.R")
library(tidyverse)
library(ggrepel)
library(ggnewscale)


#get all sample metadata
all_info<-readRDS("acral_paired/ds/v00-sample_info.rds")



#get coordinates from plot-mds analysis, and dimension length.
mds<-readRDS("acral_paired/ds/mds_dge.rds")
mds_coords<-tibble(x=mds$x,
       y=mds$y,
       sample_id=attributes(mds[["distance.matrix.squared"]])[["dimnames"]][[1]])


#join mds coordinates to sample info.
dfp<-mds_coords%>%left_join(all_info)





#get dimensions for plot based on ratio of dimensions.
x_dim<-round(mds[["var.explained"]][1]*100, 1)
y_dim<-round(mds[["var.explained"]][2]*100, 1)
dim_ratio<-round(y_dim/x_dim, 2)





bubble_plot<-function(plot_title, plot_subtitle, category, set_dim_ratio, pdf_width, pdf_height){
  plot_to_save<-ggplot(dfp, aes(x, y))+
    
    
    geom_line(aes(group=factor(mouse_num)), color="grey", alpha=.6) +
    
    
    
    
    new_scale_color()+
    ggtitle(plot_title)+
    labs(subtitle = plot_subtitle)+
    
    
    
    #layer1
    #sets bubble inside colors
    geom_point(aes(color=!!sym(category), size=7), shape=16, alpha=.5)+
    #sets scale based on color mapping df. 
    
    ####
    #scale_color_manual(values=colors[[category]])+
    
    #set up legend prior to resetting color scale.
    
    
    
    guides(color = guide_legend(override.aes = list(size = 7),
                                title=""))+
    
    new_scale_color()+
    
    #layer2
    #sets bubble outline colors
    geom_point(aes(color=!!sym(category), size=7), shape=1)+
    #sets scale based on color mapping df, and darkens. 
    
    
    
    ####
    #scale_color_manual(values=sapply(colors[[category]],function(x) colorspace::darken(x, 0.2)))+
    
    #set up legend, prior to resetting color scale. Must match above legend otherwise
    #legend will split. 
    guides(color = guide_legend(title=""))+
    
    new_scale_color()+
    
    #set size scale for bubbles
    #scale_size(range = c(2,10))+
    
    
    
    
    
    #add text labels for points, by mouse number.
    geom_text_repel(aes(label = mouse_num, color=!!sym(category)),
                    min.segment.length = 0,
                    segment.color = "grey80",
                    force=20,
                    point.padding=15,
                    show.legend = FALSE) +
    #set color of text to be slightly darker than that of bubble. 
    
    
    
    ########
  #scale_color_manual(values=sapply(colors[[category]],function(x) colorspace::darken(x, 0.3)))+
  
  #remove legend for size
  guides(size="none")+
    
    #change y and x axes labels to include percent difference in MDS dimension.
    labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
    labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
    
    #set aspect ratio - (probably either 1:1 or scale to %diff).
    coord_fixed(ratio=set_dim_ratio)+
    theme_classic()
  
  # plot_to_save
  
  
  object<-plot_to_save
  
  ggsave(paste0("acral_paired/plots/",fs::path_sanitize(paste0("mds-", category, plot_title, "unscaled", ".", set_dim_ratio, "ratio", "-m.pdf"))),
         
         
         title=paste0("src: ",
                      
                      rstudioapi::getSourceEditorContext()$path%>%
                        sub("/uufs/chpc.utah.edu/common/home/holmen-group1/otrimskig/rchpc/","",.)%>%
                        sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/rchpc/","",.)%>%
                        sub("C:/Users/u1413890/OneDrive - University of Utah/garrett hl-onedrive/R/","",.),
                      
                      " at ", 
                      
                      lubridate::round_date(Sys.time(), "second")
         ),
         
         
         plot=object,
         
         scale = .75,
         dpi=600,
         width = pdf_width,
         height = pdf_height,
         unit="in")
  
  
  
  
 
  
  
}


########################
#set up recurring plot variables here..
plot_subtitle<-""
set_dim_ratio<-dim_ratio
######################################

bubble_plot(plot_title = "By Tumorinjection site2", 
            category = "sample_type",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            # 
            # 
            pdf_width = 10,
            pdf_height = 7)







