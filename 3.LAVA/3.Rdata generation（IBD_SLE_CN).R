rm(list = ls()) 
gc()
library(LAVA)
library(tidyverse)
library(magrittr)
#--------------------------- Set Parameters --------------------------#
ref_path = "./Database/Ref_Panel/g1000_eas_chr/g1000_eas_chr" # East Asian
loci_path = "./Database/LAVA/blocks_s2500_m25_f1_w200_eas.locfile" # East Asian

filename = "SLE_IBD"
Rdata_dirout = "./LAVA/IBD_SLE_CN/Rdata"
sample.overlap = NULL
input.info_path = "./LAVA/IBD_SLE_CN/info_files/input.info_chr"
phenotype = c("SLE_CN", "IBD_eas", "UC_eas", "CD_eas")

#--------------------------- Separator Line --------------------------#
setwd(Rdata_dirout)
for (m in 1:22) {
  input.info = paste0(input.info_path, m, ".txt")
  ref = paste0(ref_path, m)
  loci = read.loci(loci_path)
  loci %<>% filter(CHR == m) 
  input = process.input(input.info.file = input.info,
                        sample.overlap.file = sample.overlap,
                        ref.prefix = ref,
                        phenos = phenotype)
  save.image(paste0(filename, "_chr", m, ".RData"))
  print(m)
}