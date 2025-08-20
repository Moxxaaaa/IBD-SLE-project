rm(list = ls()) 
gc()

library(data.table)
library(magrittr)
library(tidyverse)
library(readxl)
#------------------ Set Parameters ------------------#
code <- read_xlsx("./GWASsummary_cleaned_data/Disease_Codes_[231030].xlsx", sheet = "EAS")
code <- as.data.frame(code)
gcta_dir = "./tools/YANGLAB/gcta/exe/gcta-win-1.94.1.exe"
bfile_dir = "./tools/CTGLAB/Ref_Panel/EAS/g1000_eas/g1000_eas"
gene_list_dir = "./resource/Example/mBAT_combo/glist_symbol_hg19_v40"

traits <- c("SLE", "IBD", "CD", "UC")
for(trait in traits){
  print(trait)
  print(paste0(code[code$Phenotype == trait, "File_path"]))
  data <- fread(paste0(code[code$Phenotype == trait, "File_path"]))
  data <- subset(data, select = c(SNP, A1, A2, EAF, BETA, SE, P, N))
  fwrite(data,
         file = paste0("./202311/1115mBAT/data/", trait, ".ma"), sep = "\t")
  gwas_dir = paste0("./202311/1115mBAT/data/", trait, ".ma")
  out_dir = paste0("./202311/1115mBAT/result/", trait)
  c1 = gcta_dir
  c2 = paste0("--bfile ", bfile_dir)
  c3 = paste0("--mBAT-combo ", gwas_dir)
  c4 = paste0("--mBAT-gene-list ", gene_list_dir)
  c5 = "--mBAT-print-all-p "
  c6 = paste0("--out ", out_dir)
  c7 = "--thread-num 10"
  command = paste(c1, c2, c3, c4, c5, c6, c7, sep = " ")
  system(command)
}
# Adjust p-values
for(trait in traits){
  data <- fread(paste0("./202311/result/", trait, ".gene.assoc.mbat"))
  data$p.adjust <- p.adjust(data$P_mBATcombo, method = "BH", n = length(data$P_mBATcombo))
  fwrite(data,
         file = paste0("./202311/result/", trait, "2.gene.assoc.mbat"), sep = "\t")
}