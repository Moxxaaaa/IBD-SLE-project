rm(list = ls()) 
gc()

#--------------------------- Set Parameters --------------------------#
output_path = "./LAVA/IBD_SLE_CN/result"
loci_path = "./Database/LAVA/blocks_s2500_m25_f1_w200_eas.locfile" # East Asian
filename = "SLE_IBD"
Rdata_dir_prefix = "./LAVA/IBD_SLE_CN/Rdata/"
target_phenos = "SLE_CN"

#-----------------chr1-22---------------------------
ST <- Sys.time() # Record start time
for (chr in 1:22) {
  ncore = 8 ### Set number of cores
  source("./LAVA/IBD_SLE_CN/LAVA.parallel.chr.target_core_code(eas).R")
  print(chr)
}

ET <- Sys.time() # Record end time
print(ST - ET) # Calculate time difference