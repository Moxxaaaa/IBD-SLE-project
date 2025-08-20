library(data.table)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(TwoSampleMR)
library(phenoscanner)
library(ieugwasr)
#------------------------- Set Parameters --------------------#
out_path = "./MR/IBD_SLE/eas"
code = "./Disease_Codes_[231114].xlsx"
Sheet = "EAS"
MR_pop = "EAS"
plink_dir = "./Database/Plink/plink.exe"
bfile_dir = "./Database/Ref_Panel/g1000_eas/g1000_eas" ## Corresponding to East Asian
#---------------------- Separator Line -------------------#
setwd(out_path)
code <- read_xlsx(code, sheet = Sheet)
code <- as.data.frame(code)
code %<>% filter(Phenotype == "IBD" | Phenotype == "CD" |
                   Phenotype == "UC" | Phenotype == "SLE")

IV_list <- list() # Create instrument variable list
#-------------------- 1. SLE to IBD ---------------------
sumstats1_for = c("SLE")
sumstats2_for = c("IBD", "UC", "CD")
for (m in 1:length(sumstats1_for)) {
  sumstats1 = sumstats1_for[m]
  f1 = code[code$Phenotype == sumstats1, "File_path"] #### Get file path
  f1 <- fread(f1)
  f1$phenotype = sumstats1
  f1_exposure <- format_data(f1, type = "exposure", snps = NULL, header = TRUE,
                             phenotype_col = "phenotype",
                             snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF",
                             effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
                             chr_col = "CHR", pos_col = "POS", samplesize_col = "N",
                             min_pval = 1e-500, log_pval = FALSE)
  f1_exposure1 <- f1_exposure[f1_exposure$pval.exposure < 5*10^-8,]
  f1_exposure1 %<>% mutate(rsid = SNP) %>% mutate(pval = pval.exposure)
  f1_exposure_clump <- ld_clump(dat = f1_exposure1,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                pop = MR_pop,
                                bfile = bfile_dir,
                                plink_bin = plink_dir)
  
  #f1_exposure_clump <- clump_data(f1_exposure1, clump_r2 = 0.001, clump_kb = 10000, pop = "EUR")
  rm(f1); rm(f1_exposure); rm(f1_exposure1); gc()
  for (n in 1:length(sumstats2_for)) {
    sumstats2 = sumstats2_for[n]
    f2 = code[code$Phenotype == sumstats2, "File_path"] #### Get file path
    f2 <- fread(f2)
    f2$phenotype = sumstats2
    f2_outcome <- format_data(f2, type = "outcome", snps = NULL, header = TRUE,
                              phenotype_col = "phenotype",
                              snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF",
                              effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
                              chr_col = "CHR", pos_col = "POS", samplesize_col = "N",
                              min_pval = 1e-500, log_pval = FALSE)
    dat <- harmonise_data(exposure_dat = f1_exposure_clump, outcome_dat = f2_outcome, action = 2)
    rm(f2); rm(f2_outcome); gc()
    #------------ MR steiger ---------#
    dat <- steiger_filtering(dat)
    #--------------- F-value calculation --------------#
    dat %<>%
      mutate(R2 = (2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
      mutate(Fvalue = ((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
    
    IV_list[[paste0(sumstats1, " to ", sumstats2)]] = dat
  }
}

#-------------------- 2. IBD to SLE ---------------------
sumstats2_for = c("SLE")
sumstats1_for = c("IBD", "UC", "CD")
for (m in 1:length(sumstats1_for)) {
  sumstats1 = sumstats1_for[m]
  f1 = code[code$Phenotype == sumstats1, "File_path"] #### Get file path
  f1 <- fread(f1)
  f1$phenotype = sumstats1
  f1_exposure <- format_data(f1, type = "exposure", snps = NULL, header = TRUE,
                             phenotype_col = "phenotype",
                             snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF",
                             effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
                             chr_col = "CHR", pos_col = "POS", samplesize_col = "N",
                             min_pval = 1e-500, log_pval = FALSE)
  f1_exposure1 <- f1_exposure[f1_exposure$pval.exposure < 5*10^-8,]
  f1_exposure1 <- f1_exposure[f1_exposure$pval.exposure < 5*10^-8,]
  f1_exposure1 %<>% mutate(rsid = SNP) %>% mutate(pval = pval.exposure)
  f1_exposure_clump <- ld_clump(dat = f1_exposure1,
                                clump_kb = 10000,
                                clump_r2 = 0.001,
                                pop = MR_pop,
                                bfile = bfile_dir,
                                plink_bin = plink_dir)
  #f1_exposure_clump <- clump_data(f1_exposure1, clump_r2 = 0.001, clump_kb = 10000, pop = "EUR")
  rm(f1); rm(f1_exposure); rm(f1_exposure1); gc()
  for (n in 1:length(sumstats2_for)) {
    sumstats2 = sumstats2_for[n]
    f2 = code[code$Phenotype == sumstats2, "File_path"] #### Get file path
    f2 <- fread(f2)
    f2$phenotype = sumstats2
    f2_outcome <- format_data(f2, type = "outcome", snps = NULL, header = TRUE,
                              phenotype_col = "phenotype",
                              snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "EAF",
                              effect_allele_col = "A1", other_allele_col = "A2", pval_col = "P",
                              chr_col = "CHR", pos_col = "POS", samplesize_col = "N",
                              min_pval = 1e-500, log_pval = FALSE)
    dat <- harmonise_data(exposure_dat = f1_exposure_clump, outcome_dat = f2_outcome, action = 2)
    rm(f2); rm(f2_outcome); gc()
    #------------ MR steiger ---------#
    dat <- steiger_filtering(dat)
    #--------------- F-value calculation --------------#
    dat %<>%
      mutate(R2 = (2*(beta.exposure^2)*eaf.exposure*(1 - eaf.exposure))) %>% 
      mutate(Fvalue = ((samplesize.exposure - 1 - 1)/1)*(R2/(1 - R2)))
    
    IV_list[[paste0(sumstats1, " to ", sumstats2)]] = dat
  }
}
rm(list = setdiff(ls(), c("IV_list")))
save.image("IVs_MR.Rdata")