library(coloc)
library(data.table)
library(tidyverse)
library(magrittr)
library(readxl)
library(susieR)
#------------------------------------- Set Parameters ----------------------------------#
code = "./Project/IBD_SLE_CN/coloc/FUMA_tables/Pleiotropic_Loci.csv"
code_gwas_summary = "./Project/IBD_SLE_CN/coloc/code.xlsx"
phe_c = c("IBD_SLE_eas", "CD_SLE_eas", "UC_SLE_eas")
out_path = "./Project/IBD_SLE_CN/coloc"
all_result <- list()
#------------------------------- Separator Line -------------------------------#
setwd(out_path)
code <- read.csv(code)
code_gwas_summary <- read_xlsx(code_gwas_summary)
code <- as.data.frame(code)
code_gwas_summary <- as.data.frame(code_gwas_summary)
code %<>% separate(Locus_boundary, into = c("start", "end"), sep = "-", convert = TRUE)
code <- code[c(1:5, 8:27, 31:48, 53:64),]
for (phe in phe_c) {
  f1 = code_gwas_summary[code_gwas_summary$trait_pair == phe, "dir"]
  f1 <- fread(f1)
  f1 %>% 
    filter(EAF.f1 >= 0.5) %>% 
    mutate(MAF = 1 - EAF.f1) %>% 
    mutate(nbeta.f1 = BETA.f1 / -1) %>% 
    mutate(nbeta.f2 = BETA.f2 / -1) %>%
    mutate(a1 = A2) %>% 
    mutate(a2 = A1) %>% 
    mutate(SE.f1 = BETA.f1 / Z.f1) %>% 
    mutate(SE.f2 = BETA.f2 / Z.f2) %>% 
    select(SNP, CHR, POS, MAF, a1, a2, nbeta.f1, nbeta.f2, P.f1, P.f2, N.f1, N.f2, SE.f1, SE.f2) %>% 
    rename(A1 = a1, A2 = a2,
           BETA.f1 = nbeta.f1, BETA.f2 = nbeta.f2) -> f1.1
  f1 %>% 
    filter(EAF.f1 < 0.5) %>% 
    mutate(MAF = EAF.f1) %>% 
    mutate(nbeta.f1 = BETA.f1) %>% 
    mutate(nbeta.f2 = BETA.f2) %>%
    mutate(a1 = A1) %>% 
    mutate(a2 = A2) %>% 
    mutate(SE.f1 = BETA.f1 / Z.f1) %>% 
    mutate(SE.f2 = BETA.f2 / Z.f2) %>% 
    select(SNP, CHR, POS, MAF, a1, a2, nbeta.f1, nbeta.f2, P.f1, P.f2, N.f1, N.f2, SE.f1, SE.f2, MAF) %>% 
    rename(A1 = a1, A2 = a2,
           BETA.f1 = nbeta.f1, BETA.f2 = nbeta.f2) -> f1.2
  f1 <- rbind(f1.1, f1.2); rm(f1.1); rm(f1.2); gc()
  f1 %<>%
    mutate(varbeta.f1 = SE.f1^2) %>% 
    mutate(varbeta.f2 = SE.f2^2) %>% 
    select(-c("SE.f1", "SE.f2"))
  #------------ Calculate number of pleiotropic loci for a phenotype pair --------------------#
  code %>% filter(trait_pair == phe) -> temp
  nlocus = nrow(temp)
  #------------- Extract phenotype pair information ---------------------#
  f1_type = code_gwas_summary[code_gwas_summary$trait_pair == phe, "trait1_type"]
  f2_type = code_gwas_summary[code_gwas_summary$trait_pair == phe, "trait2_type"]
  ref_panel = code_gwas_summary[code_gwas_summary$trait_pair == phe, "ref_dir"]
  #------- Create list -----------#
  single_coloc_result <- list()
  #------------- Separator Line -----------------------------#
  
  for (i in 1:nlocus) {
    locus_start = temp[i, "start"]
    locus_end = temp[i, "end"]
    chr = temp[i, "CHR"]
    f1 %>% 
      filter(CHR == chr) %>% 
      filter(POS >= locus_start & POS <= locus_end) -> region1
    
    #------------------ Perform single-coloc analysis -------------------#
    dataset1 = list(snp = region1$SNP,
                    chr = region1$CHR,
                    pos = region1$POS,
                    beta = region1$BETA.f1, 
                    varbeta = region1$varbeta.f1,
                    type = f1_type, 
                    N = max(region1$N.f1),
                    MAF = region1$MAF)
    
    dataset2 = list(snp = region1$SNP,
                    chr = region1$CHR,
                    pos = region1$POS,
                    beta = region1$BETA.f2, 
                    varbeta = region1$varbeta.f2,
                    type = f2_type, 
                    N = max(region1$N.f2),
                    MAF = region1$MAF)
    check_dataset(dataset1)
    check_dataset(dataset2)
    singlecoloc_result <- coloc.abf(dataset1, dataset2, MAF = region1$MAF)
    single_coloc_result[[paste(chr, locus_start, locus_end, sep = ":")]] = singlecoloc_result
    print(paste0(phe, ":", i))
  }
  all_result[[phe]] = single_coloc_result
}
rm(list = setdiff(ls(), c("all_result", "code", "code_gwas_summary", "phe_c"))); gc()

save.image("single-coloc_result.RData")