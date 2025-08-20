
###### Initial File Preparation #####
rm(list = ls()) 
gc()
library(readxl)
library(data.table)
library(tidyverse)
library(magrittr)
#--------------------------- Set Parameters --------------------------#
output_path = "./LAVA/IBD_SLE_CN/info_files"
file_dir = "./LAVA/IBD_SLE_CN/Temp/"
code = "./Disease_Codes_[230729].xlsx"
Sheet = "EAS"
code = read_xlsx(code, sheet = Sheet)
code <- as.data.frame(code)
code %<>% filter(Phenotype == "SLE" | token == "Gastric_Diseases")

#--------------------------- Separator Line --------------------------#
setwd(output_path)

for (i in 1:22) {
  code %>% 
    mutate(dir1 = file_dir) %>% 
    mutate(dir = paste0(dir1, File_Prefix, "_chr", i, ".txt.gz")) %>% 
    select(File_Prefix, case, control, dir) %>% 
    rename(phenotype = File_Prefix,
           cases = case,
           controls = control,
           filename = dir) %>% 
    mutate(filename = sprintf('"%s"', filename)) -> c
  
  write.table(c,
              file = paste0("input.info_chr", i, ".txt"),
              sep = " ",
              quote = F,
              row.names = F,
              col.names = T)
}