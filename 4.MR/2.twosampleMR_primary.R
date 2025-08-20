library(data.table)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(readxl)
library(TwoSampleMR)

#-------------------------Set Parameters--------------------#
out_path = "./MR/IBD_SLE/eas/TwosampleMR"
rda = "./MR/IB_SLE/eas/IVs_MR.Rdata"
code = "./MR/IBD_SLE/code.xlsx"
MR_pop = "eas"
#----------------------Divider-------------------#
MR_result <- list()
leaveone.out <- list()
code <- as.data.frame(read_xlsx(code, sheet = MR_pop))
setwd(out)
load(rda)

#-----------------------TwoSampleMR-------------#
for (i in 1:nrow(code)) {
  ex_name = code[i, "exposure"]
  ou_name = code[i, "outcome"]
  f_name = paste(ex_name, "to", ou_name)
  dat = IV_list[[f_name]]
  dat %<>%
    filter(ambiguous == "FALSE") %>%
    filter(pval.outcome >= 5e-8)
  #--------------TwoSample MR-------------#
  res<-mr(dat)
  het<-mr_heterogeneity(dat)
  ple<-mr_pleiotropy_test(dat)
  restBindsub<-bind_rows(res,het,ple)
  leaveone_out<-mr_leaveoneout(dat)
  #-------------Save Results-------------#
  MR_result[[f_name]] = restBindsub
  leaveone.out[[f_name]] = leaveone_out
  print(i)
}

rm(list = setdiff(ls(), c("MR_result","leaveone.out")))
save.image("primary_twosample_MR.RData")
MR_result <- bind_rows(MR_result)
leaveone.out <- bind_rows(leaveone.out)
write.csv(MR_result, row.names = F,
          file = "TwosampleMR_primary.csv")
write.csv(leaveone.out,row.names = F,
          file = "leaveone.out_primary.csv")