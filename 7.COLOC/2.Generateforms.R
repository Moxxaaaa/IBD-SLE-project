library(coloc)
library(data.table)
library(tidyverse)
library(magrittr)
library(readxl)
library(susieR)
#------------------------------------- Set Parameters ----------------------------------#
rdata_dir="./coloc/single-coloc_result.RData"
phe_c=c("IBD_SLE_eas","CD_SLE_eas","UC_SLE_eas")
out_path="./IBD_SLE_CN/coloc"
#------------------------------- Separator Line -------------------------------#
setwd(out_path)
load(rdata_dir)
code %<>% mutate(id=paste(CHR,start,end,sep=":"))
code$PP.H3=NA
code$PP.H4=NA
code$Best_causal_SNP=NA
code$SNP.PP.H4=NA
temp_list<-list()
for (phe in phe_c) {
  f1=all_result[[phe]]
  code %>% filter(trait_pair==phe)->temp
  n=nrow(temp)
  for (i in 1:n) {
    id1=temp[i,"id"]
    cs=f1[[id1]][["results"]]
    cs %<>% arrange(desc(SNP.PP.H4)) 
    #-------------提取想要的数据-------------------------#
    PPH3=as.numeric(f1[[id1]][["summary"]][5])
    PPH4=as.numeric(f1[[id1]][["summary"]][6])
    Best.causal.SNP=cs[1,"snp"]
    SNP_PP_H4=cs[1,"SNP.PP.H4"]
    #------------放入应该放入的地方----------------#
    temp[temp$id==id1,"PP.H3"]=PPH3
    temp[temp$id==id1,"PP.H4"]=PPH4
    temp[temp$id==id1,"Best_causal_SNP"]=Best.causal.SNP
    temp[temp$id==id1,"SNP.PP.H4"]=SNP_PP_H4
  }
  temp_list[[phe]]=temp
}
temp<-bind_rows(temp_list)

temp %<>%
  mutate(Locus_boundary=paste(CHR,start,end,sep = ":")) %>% 
  select("trait_pair","Top_SNP","Locus_boundary","nearestGene",
         "PP.H3","PP.H4","Best_causal_SNP","SNP.PP.H4" )

write.csv(temp,
          file="single-coloc_result.csv")






