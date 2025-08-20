library(readr)
library(readxl)
library(dplyr)
library(cause)
library(magrittr)
library(data.table)

#-------------------------------------------- Set Parameters ------------------------#
out_path = "./CAUSE/IBD_SLE_CN/"
code = "./CAUSE/IBD_SLE_CN/CAUSE.xlsx"
Sheet = "eas"
#----- ld pruning ----#
plink_dir = "./Database/Plink/plink"
bfile_dir = "./Database/Ref_Panel/g1000_eas/g1000_eas"
r2_thresh = 0.1
pval_thresh = 1e-3
#----------------------- Separator Line ------------------------#
code <- read_xlsx(code, sheet = Sheet)
code = as.data.frame(code)

for (i in 1:nrow(code)) {
  ex_name = code[i, "exposure"]
  ou_name = code[i, "outcome"]
  output_name = paste0(ex_name, "_to_", ou_name)
  X1 = code[i, "exposure_dir"]; X2 = code[i, "outcome_dir"]
  X1 = fread(X1); X2 = fread(X2)
  #----------- Remove MHC region ------------#
  X1 %<>% filter(!(CHR == 6 & POS > 25000000 & POS < 35000000))
  X2 %<>% filter(!(CHR == 6 & POS > 25000000 & POS < 35000000))
  
  #-------------- Merge --------------#
  X <- gwas_merge(X1, X2, 
                  snp_name_cols = c("SNP", "SNP"), 
                  beta_hat_cols = c("BETA", "BETA"), 
                  se_cols = c("SE", "SE"), 
                  A1_cols = c("A1", "A1"), 
                  A2_cols = c("A2", "A2"), 
                  pval_cols = c("P", "P"))
  
  #--------- 2. Calculate nuisance parameters ----------#
  set.seed(100)
  varlist <- with(X, sample(snp, size = 1000000, replace = FALSE))
  params <- est_cause_params(X, varlist)
  
  class(params)
  names(params)
  params$rho
  head(params$mix_grid)
  
  #---------- LD Pruning ------------#
  X_clump <- X %>%
    rename(rsid = snp,
           pval = p1) %>%
    ieugwasr::ld_clump(dat = .,
                       clump_r2 = r2_thresh,
                       clump_p = pval_thresh,
                       plink_bin = plink_dir, 
                       pop = "eas",
                       bfile = bfile_dir)
  
  top_vars <- X_clump$rsid
  
  #---------------- 4. Fit CAUSE ---------#
  res <- cause(X = X, variants = top_vars, param_ests = params)
  res$loos[[2]]
  loo::pareto_k_table(res$loos[[2]])
  res$loos[[3]]
  loo::pareto_k_table(res$loos[[3]])
  
  #-------------- 5. save result --------------#
  if(!dir.exists(paste0(out_path, output_name))){
    dir.create(paste0(out_path, output_name))}
  setwd(paste0(out_path, output_name))
  
  rm(X); rm(X1); rm(X2); gc()
  save.image(paste0(output_name, ".RData"))
  
  #--------------- Output results ---------------#
  sink(paste0(output_name, ".log"), type = c("output", "message"))
  print(res$elpd)
  print(summary(res, ci_size = 0.95))
  sink() 
  
  pdf(paste0(output_name, ".pdf"), width = 12, height = 6)
  plot(res)
  dev.off()
  print(i)
}