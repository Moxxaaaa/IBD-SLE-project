###### Initial File Preparation #####
rm(list = ls()) 
gc()
library(readxl)
library(data.table)
library(tidyverse)
library(magrittr)
#--------------------------- Set Parameters --------------------------#
out_dir="./LAVA/IBD_SLE_CN/Temp" # Output directory prefix
code="./Disease Codes [230729].xlsx"
Sheet="EAS"
code=read_xlsx(code,sheet = Sheet)
code<-as.data.frame(code)
code %<>%filter(Phenotype=="SLE"|token == "Gastric Diseases")

phenotype=code$Phenotype  ### List of phenotype names
write.table(phenotype,paste0(out_dir,"/","filename.txt"),col.names = F,row.names = F,quote = T)

#--------------------------------- Separator Line -----
Nphenotype=length(phenotype)
setwd(out_dir)
for (i in 1:Nphenotype) {
  sumstats=code[i,"File Name Prefix"]
  f1=code[i,"File Path"]
  f1<-fread(f1)
  f1 %<>%
    mutate(Z=BETA/SE) %>% 
    select(SNP,CHR,A1,A2,N,Z)
  for (n in 1:22) {
    f1 %>% 
      filter(CHR==n) %>% 
      select(SNP,A1,A2,N,Z)->part
    fwrite(part,file = paste0(sumstats,"_chr",n,".txt.gz"),
           sep = "\t",col.names = T,row.names = F,quote = F,compress = "gzip")
    print(paste0(i,".",n))
  }
}