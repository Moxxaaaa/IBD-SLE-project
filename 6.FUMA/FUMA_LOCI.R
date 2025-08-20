library(magrittr)
library(tidyverse)
library(readxl)
library(data.table)
result <- list()
Pleiotropic_Loci <- list()
Top_loci <- list()
#---------------------------------------- Set Parameters -------------------------#
input_path = "./202310/1014SLE_IBD_EAS/PLACO/FUMA_results/"
out_path = "./results/SLE_IBD_EAS/PLACO_tables/"
code_dir = "./202310/1014SLE_IBD_EAS/PLACO/FUMA_results/Code_Correspondence.xlsx"
Sheet = "Sheet1"
type1 = "snps.txt"
type2 = "leadSNPs.txt"
type3 = "GenomicRiskLoci.txt"
annov_dir = "annov.txt"
#------------------------------- Separator Line -----------------------------#
setwd(out_path)
code <- read_xlsx(code_dir, sheet = Sheet)
code <- as.data.frame(code)
for (i in 1:4) {
  filename = code[i, "filename"]
  phe = code[i, "phe"]
  #------------- Input files -----------#
  f1 = paste0(input_path, filename, "/", type1)
  f2 = paste0(input_path, filename, "/", type2)
  f3 = paste0(input_path, filename, "/", type3)
  annov = paste0(input_path, filename, "/", annov_dir)
  f1 <- fread(f1)
  f2 <- fread(f2)
  f3 <- fread(f3)
  annov <- fread(annov)
  #----------------- IndSig SNP ------------#
  f1 %<>% 
    filter(is.na(gwasP) != T) %>% 
    filter(gwasP < 5e-8) %>% 
    filter(IndSigSNP == rsID) %>% 
    mutate(trait_pair = phe)
  f1 %<>%
    left_join(annov, by = c("uniqID", "chr", "pos")) %>% 
    arrange(dist.y) %>% 
    distinct(uniqID, .keep_all = T) %>% 
    arrange(chr, pos)
  
  #----------------- Lead SNP ------------#
  f2 %<>% 
    select(uniqID, nIndSigSNPs, IndSigSNPs) %>% 
    mutate(Lead_SNP = 1)
  
  f1 %<>% left_join(f2, by = "uniqID") 
  
  #--------------- Top SNP ----------------#
  f3 %<>% 
    select(uniqID, start, end, nGWASSNPs, nLeadSNPs, LeadSNPs) %>% 
    mutate(Top_SNP = 1)
  
  f1 %<>% left_join(f3, by = "uniqID") 
  
  #------------ Extract OR and P from single trait GWAS ------------#
  single_gwas = code[i, "PLACO_file_path"]
  single_gwas = fread(single_gwas)
  single_gwas %<>%
    rename(rsID = SNP,
           P_psy = P.f1,
           P_smk = P.f2) %>% 
    mutate(OR_psy = exp(BETA.f1)) %>% 
    mutate(OR_smk = exp(BETA.f2)) %>% 
    select(rsID, A1, A2, OR_psy, P_psy, OR_smk, P_smk)
  
  f1 %<>% left_join(single_gwas, by = "rsID")
  
  #-------------- Organize ----------------#
  f1 %>% 
    mutate(Locus_boundary = paste0(start, "-", end)) %>% 
    select("trait_pair", "GenomicLocus",
           "uniqID", "rsID", "chr", "pos", "non_effect_allele", "effect_allele", "MAF", "gwasP",
           "nearestGene", "dist.x", "func", "CADD", "RDB",    
           "gene", "symbol", "annot", "dist.y", "exonic_func", "exon", # ANNOVAR annotation results
           "Locus_boundary", "start", "end", "Lead_SNP", "Top_SNP",
           "nGWASSNPs", "nIndSigSNPs", "IndSigSNPs", "nLeadSNPs", "LeadSNPs",
           "A1", "A2", "OR_psy", "P_psy", "OR_smk", "P_smk") -> temp
  
  result[[phe]] = temp
  
  #------------- Organize into Pleiotropic_Loci format ---------#
  temp %<>% 
    filter(Top_SNP == 1) %>% 
    select("trait_pair", "rsID", "chr", "pos", "Locus_boundary", "gwasP",
           "A1", "A2", "nearestGene", "dist.x", "func", "CADD", "RDB", "OR_psy", "OR_smk", "P_psy", "P_smk") %>% 
    rename(Top_SNP = rsID,
           CHR = chr,
           POS = pos,
           P_placo = gwasP,
           dist = dist.x,
           Functional_annotation = func)
  
  n = nrow(temp)
  temp$Sig_psy = NA
  temp$Sig_smk = NA
  for (m in 1:n) {
    p1 = temp[m, "P_psy"]
    p2 = temp[m, "P_smk"]
    if(p1 < 5e-8){
      temp[m, "Sig_psy"] <- 1
    }else{
      temp[m, "Sig_psy"] <- 0
    }
    if(p2 < 5e-8){
      temp[m, "Sig_smk"] <- 1
    }else{
      temp[m, "Sig_smk"] <- 0
    }
  }
  temp$Sig_psy <- as.numeric(temp$Sig_psy)
  temp$Sig_smk <- as.numeric(temp$Sig_smk)
  
  temp %>% 
    select("trait_pair", "Top_SNP", "CHR", "POS", "Locus_boundary", "P_placo",
           "A1", "A2", "nearestGene", "Functional_annotation", "CADD", "RDB", "Sig_psy", "Sig_smk") -> temp1
  
  Pleiotropic_Loci[[phe]] = temp1
  
  #------------------- Top SNP ---------------------#
  temp %>% 
    select("trait_pair", "Top_SNP", "CHR", "Locus_boundary",
           "A1", "A2", "OR_psy", "OR_smk", "P_psy", "P_smk") -> temp2
  
  Top_loci[[phe]] <- temp2
  
  rm(single_gwas)
  print(i)
}

rm(list = setdiff(ls(), c("result", "Pleiotropic_Loci", "Top_loci")))
save.image("FUMA_loci_result.RData")

Pleiotropic_Loci <- bind_rows(Pleiotropic_Loci)
Top_loci <- bind_rows(Top_loci)

write.csv(Pleiotropic_Loci,
          file = "Pleiotropic_Loci.csv",
          row.names = F)

write.csv(Top_loci,
          file = "Top_loci.csv",
          row.names = F)