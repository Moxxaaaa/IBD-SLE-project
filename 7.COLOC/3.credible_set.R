library(magrittr)
library(tidyverse)

#------------------------------------- Set Parameters ----------------------------------#
rdata_dir="./coloc/single-coloc_result.RData"
coloc_result="./coloc/single-coloc_result.csv"
out_path="./IBD_SLE_CN/coloc"
cs_res<-list()
#------------------------------- Separator Line -------------------------------#
setwd(out_path)
load(rdata_dir)
coloc_result<-read.csv(coloc_result)
coloc_result %<>% 
  filter(PP.H4>=0.9) %>% 
  arrange(desc(PP.H4))

for (i in 1:nrow(coloc_result)) {
  tarit_pair=coloc_result[i,"trait_pair"]
  Locus_boundary=coloc_result[i,"Locus_boundary"]
  temp<-all_result[[tarit_pair]][[Locus_boundary]]$results
  temp %<>%
    arrange(desc(SNP.PP.H4))%>% 
    mutate(CumulativeSum=cumsum(SNP.PP.H4))
  
  temp %>% filter(CumulativeSum<0.99) %>% nrow()->a

  temp<-temp[1:(a+1),]
  temp$traitpair=tarit_pair
  temp$Locus_boundary=Locus_boundary
  cs_res[[tarit_pair]][[Locus_boundary]]=temp

  
}

cs1<-bind_rows(cs_res[[1]])
cs2<-bind_rows(cs_res[[2]])
cs3<-bind_rows(cs_res[[3]])
cs<-rbind(cs1,cs2,cs3)
rm(list = setdiff(ls(),c("cs","cs_res")))
save.image("cs_result.RData")
write.csv(cs,row.names = F,
          file = "cs_result.csv")
