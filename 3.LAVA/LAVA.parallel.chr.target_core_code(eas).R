
library(LAVA)
library(magrittr)
library(tidyverse)
#---------------chr相关
loci = read.loci(loci_path)
#table(loci$CHR)
nlocus<-scan(text ="1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22 
167 177 149 143 133 149 128 123  99 117 115 117  76  68  67  66  64  64  61  49  32  33 ", what = numeric())
nlocus<- matrix(nlocus, ncol = 22, byrow = TRUE)
nlocus<-as.data.frame(t(nlocus))
nlocus %<>% rename(CHR=V1,Nlocus=V2)
locus_end=nlocus[nlocus$CHR==chr,"Nlocus"]
Rdata_path=paste0(Rdata_dir_prefix,
                  filename,
                  "_chr",
                  chr,
                  ".RData")

#----------------并行计算------------------
library(parallel)
u =b=NULL
cl <- makeCluster(ncore)
clusterExport(cl,c("Rdata_path","u","b","target_phenos"))
clusterEvalQ(cl, {library(LAVA)
  load(Rdata_path)})

start_time <- Sys.time() # 记录初始时间
out <- clusterApply(cl,1:locus_end, function(i) {
  locus <- process.locus(loci[i,], input) # process locus
  if (target_phenos %in% locus$phenos){
    ub<- run.univ.bivar(locus,univ.thresh=0.05,target = target_phenos)
    loc.info <- data.frame(locus = locus$id, chr = locus$chr, start = locus$start, 
                           stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
    u<- cbind(loc.info, ub$univ)
    if (!is.null(ub$bivar)) {
      b <- cbind(loc.info, ub$bivar)
    }
  }
  return(list(univ=u, bivar=b))
})
stopCluster(cl)# stop the cluster of worker processes
end_time <- Sys.time() # 记录终止时间
print(end_time - start_time) # 计算时间差
rm(cl)
#保存
setwd(output_path)
save.image(paste0(filename,"_chr",chr,".RData"))
write.table(do.call(rbind, lapply(out,"[[","univ")), paste0(filename,"_chr",chr,".univ.txt"), row.names=F, quote=F,sep = "\t")
write.table(do.call(rbind, lapply(out, "[[", "bivar")), paste0(filename,"_chr",chr,".bivar.txt"), row.names=F, quote=F,sep = "\t")
rm(out)
gc()
