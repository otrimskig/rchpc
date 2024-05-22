library(tidyverse)
library(ggrepel)
library(ggnewscale)


#get all sample metadata
all_info<-readRDS("ds/v10-per_sample_updated.rds")


#get coordinates from plot-mds analysis, and dimension length.
mds<-readRDS("ds/mds_dge.rds")
mds_coords<-tibble(x=mds$x,
       y=mds$y,
       sample_id=attributes(mds[["distance.matrix.squared"]])[["dimnames"]][[1]])
x_dim<-round(mds[["var.explained"]][1]*100, 1)
y_dim<-round(mds[["var.explained"]][2]*100, 1)

dim_ratio<-round(y_dim/x_dim, 2)

#join mds coordinates to sample info.
dfp<-mds_coords%>%left_join(all_info)



#source("colors_input_from_gsheets.R")
colors<-readRDS("ds/colors_list.rds")


#make colors darker (depends on colorspace)
# sapply(colors$patho_cat_name, function(x) colorspace::darken(x, 0.5))



bubble_plot<-function(plot_title, plot_subtitle, category, set_dim_ratio, pdf_width, pdf_height){
plot_to_save<-ggplot(dfp, aes(x, y))+
  
  ggtitle(plot_title)+
  labs(subtitle = plot_subtitle)+
  
  
  
  #layer1
  #sets bubble inside colors
  geom_point(aes(color=!!sym(category), size=(150-aod+1)), shape=16, alpha=.5)+
  #sets scale based on color mapping df. 
  scale_color_manual(values=colors[[category]])+
  
  #set up legend prior to resetting color scale.
  guides(color = guide_legend(override.aes = list(size = 5),
                              title=""))+
  
  new_scale_color()+
  
  #layer2
  #sets bubble outline colors
  geom_point(aes(color=!!sym(category), size=(150-aod+1)), shape=1)+
  #sets scale based on color mapping df, and darkens. 
  scale_color_manual(values=sapply(colors[[category]],function(x) colorspace::darken(x, 0.2)))+
  
  #set up legend, prior to resetting color scale. Must match above legend otherwise
  #legend will split. 
  guides(color = guide_legend(title=""))+
 
  new_scale_color()+

  #set size scale for bubbles
  scale_size(range = c(2,10))+
 
  
  
  

  #add text labels for points, by mouse number.
  geom_text_repel(aes(label = mouse_num, color=!!sym(category)),
                  min.segment.length = 0,
                  segment.color = "grey80",
                  force=20,
                  point.padding=15,
                  show.legend = FALSE) +
  #set color of text to be slightly darker than that of bubble. 
  scale_color_manual(values=sapply(colors[[category]],function(x) colorspace::darken(x, 0.3)))+

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

ggsave(paste0("plots/",fs::path_sanitize(paste0("mds-", category, set_dim_ratio, "ratio", "-m.pdf"))),
       plot=object,
       
       scale = .75,
       dpi=600,
       width = pdf_width,
       height = pdf_height,
       unit="in")


}


########################
#set up recurring plot variables here..
plot_subtitle<-"size of point scaled as days remaining to experimental endpoint upon death"
set_dim_ratio<-dim_ratio
######################################
bubble_plot(plot_title = "Pathologist Major Categorization", 
            category = "patho_cat_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)

bubble_plot(plot_title = "Pathologist Detailed Categorization", 
            category = "patho_cat_det_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)


bubble_plot(plot_title = "Pathologist Categorization and Grade", 
            category = "patho_cat2_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)


bubble_plot(plot_title = "Pathologist Grade", 
            category = "patho_grade",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)

bubble_plot(plot_title = "Resultant Genotype", 
            category = "resultant_geno",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)

bubble_plot(plot_title = "Treatment Cohorts", 
            category = "resultant_geno_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)


#########################################
#re-run same as above, but with different ratio.
set_dim_ratio<-1
pdf_width=15
######################################
bubble_plot(plot_title = "Pathologist Major Categorization", 
            category = "patho_cat_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = pdf_width,
            pdf_height = 7)

bubble_plot(plot_title = "Pathologist Detailed Categorization", 
            category = "patho_cat_det_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = pdf_width,
            pdf_height = 7)


bubble_plot(plot_title = "Pathologist Categorization and Grade", 
            category = "patho_cat2_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = pdf_width,
            pdf_height = 7)


bubble_plot(plot_title = "Pathologist Grade", 
            category = "patho_grade",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = pdf_width,
            pdf_height = 7)

bubble_plot(plot_title = "Resultant Genotype", 
            category = "resultant_geno",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 30,
            pdf_height = 7)

bubble_plot(plot_title = "Treatment Cohorts", 
            category = "resultant_geno_name",
            
            plot_subtitle = plot_subtitle,
            set_dim_ratio = set_dim_ratio,
            
            pdf_width = 40,
            pdf_height = 7)








#previous code not in use.
######################################################################



# 
# #plot 2: pathologist sub category, scaled to aod.
# plot_patho_subcat<-ggplot(dfp, aes(x, y, color=patho_cat_det_name, fill=patho_cat_det))+
#   ggtitle("Pathologist Sub-categorization")+
#   
#   geom_text_repel(aes(label = mouse_num),
#                   min.segment.length = 0,
#                   segment.color = "grey80",
#                   force=20,
#                   point.padding=15,
#                   show.legend = FALSE) +
#   
#   geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
#   geom_point(aes(size=(150-aod+1)), shape=1)+
#   
#   scale_size(range = c(2,10))+
#   
#   
#   labs(subtitle = "size of point scaled as days remaining to experimental endpoint upon death")+
#   guides(color = guide_legend(override.aes = list(size = 5)))+
#   guides(fill="none")+
#   guides(size="none")+
#   
#   labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
#   labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
#   
#   labs(color="Category")+
#   coord_fixed(ratio=dim_ratio)+
#   theme_classic()
# 
# #plot_patho_subcat
# 
# 
# #plot 3: resultant "genotype" scaled to aod.
# plot_resultant_geno<-ggplot(dfp, aes(x, y, color=cohort_proper_name, fill=resultant_geno))+
#   ggtitle("Resultant Genotype")+
#   
#   geom_text_repel(aes(label = mouse_num),
#                   min.segment.length = 0,
#                   segment.color = "grey80",
#                   force=20,
#                   point.padding=15,
#                   show.legend = FALSE) +
#   
#   geom_point(aes(size=(150-aod+1)), shape=16, alpha=.5)+
#   geom_point(aes(size=(150-aod+1)), shape=1)+
#   
#   scale_size(range = c(2,10))+
#   
#   
#   labs(subtitle = "size of point scaled as days remaining to experimental endpoint upon death")+
#   guides(color = guide_legend(override.aes = list(size = 5)))+
#   guides(fill="none")+
#   guides(size="none")+
#   
#   labs(x=paste0("Dim X (", x_dim, "% of difference)"))+
#   labs(y=paste0("Dim Y (", y_dim, "% of difference)"))+
#   
#   labs(color="Category")+
#   coord_fixed(ratio=dim_ratio)+
#   theme_classic()
# 
# 
# 

# 
# 
# to_export<-grep("^plot", ls(), value = TRUE)
# 
# for (i in 1:length(to_export)){
# 
# name<-to_export[i]
#   
# object<-get(sym(name))
# 
# ggsave(paste0("plots/",fs::path_sanitize(paste0("mds-", name, "-m.pdf"))),
#        plot=object,
#        
#        scale = .75,
#        dpi=600,
#        width = 40,
#        height =7,
#        unit="in"
#        
# )
# }
# 
