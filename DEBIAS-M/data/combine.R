library(dplyr)
library(biomformat)
library(dada2)
library(tidyverse)
# library(ggplot2)

biom.files <- list.files(pattern ="all.biom$", recursive = TRUE) #"C:/Users/brean/Downloads/masters/Fodor/DEBIAS-M/data/all" contains all the all.biom files


OTU_reads <- read_hdf5_biom(biom.files[1]) 
OTU_reads = biom(OTU_reads)
OTU_table <- as.data.frame(as.matrix(biom_data(OTU_reads)))
combine_otus <- OTU_table


seqs <- row.names(OTU_table)
taxa <- assignTaxonomy(seqs, "C:/Users/brean/Downloads/masters/Taxonomy_Dada2/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)



df <- as.data.frame(taxa) %>% select(Genus) 
combine_otus <- merge(df, OTU_table, by="row.names", all=TRUE)
combine_otus <- drop_na(combine_otus, Genus)
row.names(combine_otus) = combine_otus$Row.names
combine_otus <- subset(combine_otus, select=-c(Row.names))  
# combine_otus <- aggregate(.~Genus, data= combine_otus, FUN=sum)
combine_otus <- aggregate(.~Genus, data = combine_otus, FUN=sum)
#combine_otus <- combine
combine <- combine_otus

for (file in biom.files[2:length(biom.files)]){
  print(file)
  OTU_reads <- read_hdf5_biom(file) 
  OTU_reads = biom(OTU_reads)
  OTU_table <- as.data.frame(as.matrix(biom_data(OTU_reads)))
  
  
  seqs <- row.names(OTU_table)
  taxa <- assignTaxonomy(seqs, "C:/Users/brean/Downloads/masters/Fodor/Taxonomy_Dada2/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
  
  df <- as.data.frame(taxa) %>% select(Genus) %>% drop_na(Genus) %>% merge(OTU_table, by="row.names", all=TRUE)
  
  
  row.names(df) = df$Row.names
  df <- subset(df, select=-c(Row.names))  
  # df <- aggregate(.~Genus, data= df, FUN=sum)
  df <- aggregate(.~Genus, data= df, FUN=sum)
  print(head(df))

  combine_otus <- merge(combine_otus, df, by="Genus", all=TRUE)
  print(head(combine_otus))  
  
  print("done")

} 


#### open reference tables #####



write.csv(combine_otus, "combine_otus.csv", row.names = FALSE)
combined_saved <- combine_otus #dim: 1310, 5360 
combine_otus <- read.csv("./csv_files/combine_otus.csv", check.names=FALSE) 

combine_otus[is.na(combine_otus)] <- 0

ref.files <- list.files(pattern =".txt$", recursive = TRUE)


#### function 1: removed undesired samples (columns) ####


columns_to_filter <- function(df, columns, ref.file_num, condition){
  ref <- read.csv(ref.files[ref.file_num], sep="\t", header=TRUE)
  
  ref <- ref[,columns]

  r2 <- ref[ref[[columns[2]]] != condition,]
  
  print(dim(r2))
  r <- df[,!(names(df) %in% r2[[columns[1]]])]
  return(r)

}

# file 11210: #should have 1312, 5265 after

r <- columns_to_filter(combine_otus, c("sample_name", "env_material"), 3, "feces")

# file 11129:

r <- columns_to_filter(r, c("sample_name", "empo_3"), 2, "Animal distal gut")

# file 10959:

r <- columns_to_filter(r, c("sample_name", "geo_loc_name"), 1, "Ghana")
# r <- columns_to_filter(r, c("sample_name", "empo_3"), 1, "Animal distal gut")

# file 11635
r <- columns_to_filter(r, c("sample_name", "empo_3"), 4, "Animal distal gut")

# file 11712
r <- columns_to_filter(r, c("sample_name", "empo_3"), 6, "Animal distal gut")

# file 13241
r <- columns_to_filter(r, c("sample_name", "empo_3"), 10, "Animal distal gut")
r <- columns_to_filter(r, c("sample_name", "timepoint"), 10, "timepoint 2")


#### 1: (in df r) switch rows and columns ####

rownames(r) <- r$Genus
r <- r %>% subset(select = -c(Genus)) %>% t() %>% as.data.frame()
r <- tibble::rownames_to_column(r, "sample_name")
r$Phenotype <- NA
r$Study_ID <- NA
d.f <- r

##### function 2-3: remove undesirable values and append to df r ####

info <- function(ref.file_num, display_cols=NULL){
  d <- read.csv(ref.files[ref.file_num], sep="\t", header=TRUE, colClasses = c("sample_name" = "character"))
  
  if(!is.null(display_cols)){
    print(colnames(d))
  }
  
  print(ref.files[ref.file_num])  
  
  return(d)
}




# col_names[1] MUST be the column with the phenotypes
info_col <- function(ref.file_num, col_names, study_ID, ref.df, undesirable = NULL){
  # read file
  ref <- read.csv(ref.files[ref.file_num], sep="\t", header=TRUE, colClasses = c("sample_name" = "character"))

  # subset columns need: sample_Id, phenotype
  ref <- ref[,col_names]

  ref <- ref[ref[[col_names[2]]] %in% r[["sample_name"]], ]

  # remove undesirable samples
  if(!is.null(undesirable)){
    temp <- ref[ref[[col_names[1]]] == undesirable,]

    ref <- ref[ref[[col_names[1]]] != undesirable,]
    r <- r[!r[["sample_name"]] %in% temp[[col_names[2]]], ]


  }
  
  #add study ID column

  ref$Study_ID <- study_ID

  #  change reference column to 1s and 0s
  ref <- merge(ref, ref.df, by.x = col_names[1], by.y = "k", all.x=TRUE)

  ref <- ref[, 2:length(ref)]
  

  # merge ref to r
  inters <- intersect(names(r), names(ref))
  r <- merge(r, ref, by = inters, all = TRUE) 
  r <- aggregate(.~sample_name, data= r, FUN=sum, na.rm=TRUE, na.action = NULL)
  # r <- aggregate(.~sample_name, data= r, FUN=sum, na.rm = TRUE, na.action = NULL)
  
  
  return(r)
}




#### files ####
# 10959:
j <- info(1)
count(j, case_control__rtl_case_or_contro)

d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))

r <- info_col(1, c("case_control__rtl_case_or_contro", "sample_name"), 10959, d, "not provided")


# 11129
j <- info(2)
count(j, analysis_disease_stage)
d <- data.frame(k = c("non-t1d + normal weight", "non-t1d + overweight", "t1d + normal weight", "t1d + overweight"),
                Phenotype = c(0, 0, 1, 1))
r <- info_col(2, c("analysis_disease_stage", "sample_name"), 11129, d)


# 11210:

j <- info(3)
count(j, bmi_age46)

# j$bmi <- case_when(
#                    j$bmi_age46 == "not provided" ~ "not provided",
#                    j$bmi_age46 == "not applicab" ~ "not provided",
#                    j$bmi_age46 < 25 ~ "control",
#                    j$bmi_age46 >= 25 ~ "case")
# 
# write.table(j, ref.files[3], sep="\t", col.names = TRUE)

j <- info(3)

d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(3, c("bmi", "sample_name"), 11210, d, "not provided")

# 11635:

j <- info(4)
count(j, nafld)

# # determining whether 1 or 0 is nafld:
# # I used hba1c and bp_diastolic
# # group 1 had higher bp and hba1c
# # nafld == 1 no_nafld == 0
# h <- j %>% select(nafld, hba1c) %>% filter(hba1c != "not applicable" & hba1c != "not provided" & hba1c != "not collected")
# h$hba1c <- as.numeric(h$hba1c)
# h %>% group_by(nafld) %>% summarize(mean = mean(hba1c, na.rm=TRUE))

d <- data.frame(k = c(1, 0),
                Phenotype = c(1, 0))

r <- info_col(4, c("nafld", "sample_name"), 11635, d, "not collected")


# 11710 2 separate .biom files! 146 total

j <- info(5)
count(j, diagnosis)

d <- data.frame(k = c("SZ", "HC"),
                Phenotype = c(1, 0))

r <- info_col(5, c("diagnosis", "sample_name"), 11710, d, "") # 96 blank samples that weren't in the .biom files, so no filtering will occur


# 11712

# j <- info(6)
# count(j, any_formula)

d <- data.frame(k = c("yes", "no"),
                Phenotype = c(1, 0))
r <- info_col(6, c("any_formula", "sample_name"), 11712, d, "not collected")



# 11993
j <- info(7)
count(j, bmi)

# j$bmi <- case_when(
#                    j$bmi_class == "Lean" ~ "control",
#                    j$bmi_class == "Obese" ~ "case",
#                    j$bmi_class == "Overweight" ~ "case")
# 
# write.table(j, ref.files[7], sep="\t", col.names = TRUE)

d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(7, c("bmi", "sample_name"), 11993, d) 


# 13010
j <- info(8)
count(j, comparison_group)
# j$infected <- case_when(
#   j$comparison_group == "Both" ~ "case",
#   j$comparison_group == "GeohelmnintosPositivos" ~ "case",
#   j$comparison_group == "MalariaPositivos" ~ "case",
#   j$comparison_group == "a_Uninfected" ~ "control")
# 
# write.table(j, ref.files[8], sep="\t", col.names = TRUE)

count(j, infected)
d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(8, c("infected", "sample_name"), 13010, d) 

# 13187 
j <- info(9)
count(j, category)

# j$group <- case_when(
#   j$category == 1 ~ "case",
#   j$category == 2 ~ "control",
#   j$category == 3 ~ "control")
# 
# write.table(j, ref.files[9], sep="\t", col.names = TRUE)

d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(9, c("group", "sample_name"), 13187, d)



# 13241
# j <- info(10)
# count(j, diagnosis)
d <- data.frame(k = c("asthma", "control"),
                Phenotype = c(1, 0))
r <- info_col(10, c("diagnosis", "sample_name"), 13241, d, "not applicable") # amount of samples do not change since I filtered not app. samples with columns_to_filter()


# 13338

j <- info(11)
count(j, cohort)

# j$treatment <- case_when(
#   j$cohort == "HIV treated LD" ~ "case",
#   j$cohort == "HIV treated" ~ "case",
#   j$cohort == "HIV negative MSM" ~ "control",
#   j$cohort == "HIV negative MSW" ~ "control",
#   j$cohort == "HIV untreated" ~ "untreated")
# count(j, treatment)
# 
# write.table(j, ref.files[11], sep="\t", col.names = TRUE)

d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(11, c("treatment", "sample_name"), 13338, d, "untreated") 


# 13631

j <- info(12)
count(j, empo_3)
d <- data.frame(k = c("0", "Age_Sex_Match"),
                Phenotype = c(1, 0))
r <- info_col(12, c("control_relation", "sample_name"), 13631, d) 

# 13652
j <- info(13)
count(j, control_relation)
# j$group <- case_when(
#     j$control_relation == "subject" ~ "case",
#     j$control_relation == "Mother" ~ "control",
#     j$control_relation == "Age_and_Gender" ~ "control")
# 
# write.table(j, ref.files[13], sep="\t", col.names = TRUE)
d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(13, c("group", "sample_name"), 13652, d)


# 13695

j <- info(14)
count(j, control_relation)
d <- data.frame(k = c("subject", "sibling"),
                Phenotype = c(1, 0))
r <- info_col(14, c("control_relation", "sample_name"), 13695, d)

# 14130

j <- info(15)
count(j, host_disease)
d <- data.frame(k = c("esophageal squamous cell carcinoma", "healthy"),
                Phenotype = c(1, 0))
r <- info_col(15, c("host_disease", "sample_name"), 14130, d)

# 14669

j <- info(16)
count(j, chemotherapy)
d <- data.frame(k = c("Breast Cancer", "Healthy Control"),
                Phenotype = c(1, 0))
r <- info_col(16, c("chemotherapy", "sample_name"), 14669, d)

# 14812

# j <- info(17)
# count(j, covid_positive)
# 
# # j$group <- case_when(
# #   j$covid_positive == "no" ~ "case",
# #   j$covid_positive == "yes" ~ "control",
# #   j$covid_positive == "recovered" ~ "control")
# # 
# # write.table(j, ref.files[17], sep="\t", col.names = TRUE)
# j <- info(17)
# count(j, group)
# 
d <- data.frame(k = c("case", "control"),
                Phenotype = c(1, 0))
r <- info_col(17, c("group", "sample_name"), 14812, d)


# 15006

j <- info(18)
count(j, description)

d <- data.frame(k = c("celiac", "control"),
                Phenotype = c(1, 0))
r <- info_col(18, c("description", "sample_name"), 15006, d)

#### Remove bacteria with no observations ####
# now that filtering of samples occurred, remove columns with no samples:

save <- r
final_df <- r


filtered_cols = c()
filtered_indicies = c()
final_df = r %>% filter(Study_ID != 13241) %>% filter(Study_ID != 14812) %>% filter(Study_ID != 11712)

for(i in 4:ncol(r)){
  if(sum(as.array(r[[colnames(r)[i]]]), na.rm = TRUE) <= 0){
    filtered_cols <- append(filtered_cols, colnames(r)[i])
    filtered_indicies <- append(filtered_indicies, i)
  } 
  
  
}

dim(r[,colSums(r[,4:length(r)]) <= 0])

filtered_cols <- data.frame("columns" = filtered_cols, "indices" = filtered_indicies)

final_df <- final_df[,-c(filtered_indicies)]




#### write output to file ####
final.df = r
final.df = r %>% filter(Study_ID != 13241) %>% filter(Study_ID != 14812) %>% filter(Study_ID != 11712)

final_df <- final.df %>% select(where(~ is.numeric(.x) && sum(.x) !=0 ))
final_df <- add_column(final_df, sample_name=final.df$sample_name, .before = colnames(final_df)[1])
write.csv(final_df, "./csv_files/DEBIAS-M.csv", row.names = FALSE)



#### for troubleshooting ####

ref <- read.csv(ref.files[18], sep="\t", header=TRUE, colClasses = c("sample_name" = "character"))
rm2 <- read_hdf5_biom(biom.files[19]) 
rm2 = biom(rm2)
rm2 <- as.data.frame(as.matrix(biom_data(rm2)))
rm2 <- rm2 %>% t() %>% as.data.frame() %>% tibble::rownames_to_column("sample_name")
# ref <- filter(ref, empo_3 == "Animal distal gut" & timepoint == "timepoint 2")
sm <- merge(rm2, ref, by = "sample_name")

ref[!(ref$sample_name %in% rm2$sample_name),]

write.csv(ref, "delete2.csv", row.names = FALSE)



