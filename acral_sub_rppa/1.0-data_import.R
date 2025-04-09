source("libs.R")
library(readxl)


#read in all relevant info.
data_file<-"acral_sub_rppa/ds/16_Vashisht_Nanda_A__Kaley_Loftin_A_A.xlsx"

sheet_names<-excel_sheets(data_file)



sample_info <- read_excel("acral_sub_rppa/ds/RPPA Sample List_11.20.25_Holmen.xlsx", 
                                               sheet = "SampleList")



sheet_dfs<-list()
for(s in 1:length(sheet_names)){
  
  sheet_dfs[[sheet_names[s]]]<-read_excel(data_file, sheet=sheet_names[s])
  
}




#I: get ab metadata.
#write ab info from Lx sheets, which are identical across sheets.
ab<-1

original_sheet_name<-names(sheet_dfs)[ab]

output_name<-paste0("acral_sub_rppa/ds/", "ab-info.rds")

ab_df<-sheet_dfs[[ab]]%>%
  slice(1:7)%>%
  select(where(~ !all(is.na(.))))%>%
  
  t()%>%
  as_tibble(rownames = "ab_name", .name_repair = "unique")%>%
  janitor::row_to_names(1)%>%
  janitor::clean_names()




ab_qc<-sheet_dfs[["Antibody QC Scores"]]%>%
  janitor::clean_names()%>%
  rename(antibody_full_name=antibody_name)



ab_df2<-ab_df%>%
  left_join(ab_qc)


saveRDS(ab_df2, output_name)









#II: get all RPPA values per sample. 
#write rppa values from all Lx sheets.
rppa_data<-list()



for (la in 1:length(sheet_dfs)){


original_sheet_name<-names(sheet_dfs)[la]

cat(original_sheet_name, "\n")


if (grepl("^L\\d", original_sheet_name)) {
  
rppa_data[[original_sheet_name]]$original_df<-sheet_dfs[[la]]%>%
  slice(-c(1:7))%>%
  as_tibble(.name_repair = "unique")%>%
  janitor::row_to_names(1)



rppa_data[[original_sheet_name]]$df_long<-rppa_data[[original_sheet_name]]$original_df%>%
  
  
  select(where(~ !all(is.na(.))))%>%
  select(-c(`Sample Source`, `Category_1`, Sample, `Sample Name`,
            `Sample description`, `Sample Type`))%>%
  
  
  rename(order=Order)%>%
  mutate(order=as.numeric(order))%>%
  pivot_longer(cols = -order, names_to = "antibody_name", values_to = "rppa_value")







cat(blue("wrote", original_sheet_name, "to list\n"))



} else {
  
cat(red("not an ab info sheet\n"))
  
  
}


}



saveRDS(rppa_data, "acral_sub_rppa/ds/rppa_df_list0.rds")













#compile sample info, grouping and metadata.
sample_info0<-sample_info%>%
  janitor::clean_names()%>%
  mutate(mouse_num=substr(sample_name, 1,5))%>%
  rename(sample_num=sample_number)%>%
  mutate(sample_type=if_else(grepl("FP", sample_name), "FP", 
                             if_else(grepl("SQ", sample_name), "SQ", "none")))


rds_files<-list.files("acral_sub_rppa/ds", pattern="rppa-df.rds$", full.names = T)

for(r in 1:length(rds_files)){}

r<-1

sample_infob<-readRDS(rds_files[r])%>%
  select(Order:`Sample Type`)%>%
  mutate(sample_id=as.numeric(str_extract(`Sample Name`, "\\d+$")))%>%
  mutate(mouse_num=substr(`Sample description`, 1,5))%>%
  mutate(a_b=str_sub(`Sample description`, -1))%>%
  select(where(~ !all(is.na(.))))%>%
  
  
  janitor::clean_names()%>%
  rename(exp_sample_name=sample_name,
         sample_desc=sample_description)%>%
  select(-sample_type)%>%
  mutate(sample_type=if_else(grepl("FP", sample_desc), "FP", 
                             if_else(grepl("SQ", sample_desc), "SQ", "none")))


sample_info2<-left_join(sample_info0, sample_infob)


sample_qc<-sheet_dfs[["Sample QC metrics"]]%>%
  tibble()%>%
  janitor::clean_names()%>%
  select(where(~ !all(is.na(.))))%>%
  select(sample_description, starts_with("total"))%>%
  rename(sample_desc=sample_description)


sample_info3<-sample_info2%>%
  left_join(sample_qc)%>%
  mutate(order=as.numeric(order))



saveRDS(sample_info3, "acral_sub_rppa/ds/sample_info0.rds")
