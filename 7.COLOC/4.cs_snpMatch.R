library(coloc)
library(data.table)
library(tidyverse)
library(magrittr)
library(readxl)
library(susieR)
#------------------------------------- Set Parameters ----------------------------------#
code="./IBD_SLE_CN/coloc/cs_result.csv"
code_gwas_summary="./IBD_SLE_CN/coloc/code.xlsx"
phe_c=c("IBD_SLE_eas","CD_SLE_eas","UC_SLE_eas")
out_path="./IBD_SLE_CN/coloc"
all_result<-list()
#------------------------------- Separator Line -------------------------------#
setwd(out_path)
code<-read.csv(code)
code_gwas_summary<-read_xlsx(code_gwas_summary)
code<-as.data.frame(code)
code_gwas_summary<-as.data.frame(code_gwas_summary)

for (phe in phe_c) {
  f1=code_gwas_summary[code_gwas_summary$trait_pair==phe,"dir"]
  f1<-fread(f1)
  code %>% 
    filter(traitpair==phe) %>% 
    select(snp,SNP.PP.H4,CumulativeSum,traitpair) %>% 
    rename(SNP=snp)->tmp
  tmp<-inner_join(tmp,f1,by="SNP")
  all_result[[phe]]<-tmp
}
all_result<-bind_rows(all_result)
write.csv(all_result,
          file = "4.cs_snp_gwas.csv")
fwrite(all_result,
      file = "4.cs_snp_gwas.txt", sep = "\t",col.names = T,row.names = F,
      quote = F)


