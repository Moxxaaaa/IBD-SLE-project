rm(list=ls())
library(data.table)
library(magrittr)
library(tidyverse)
library(readxl)
library(corrplot)
library(RColorBrewer)
#------------------------------------ Set Parameters -----------------------------------#
input_path = "./LAVA/result/" # Validation results
output_path = "./initial_results/LAVA"
input.filename.prefix = "SLE_IBD_chr"
output.filename.prefix = "SLE_IBD_eas_ellipse"
code = "./Database/Disease_Codes_[240112].xlsx"
Sheet = "EAS"
#----------------------------------- Read Files -----------------------------------
h2 <- list()
local.gc <- list()
for (i in 1:22) {
  f.h2 = paste0(input_path, input.filename.prefix, i, ".univ.txt")
  f.local.gc = paste0(input_path, input.filename.prefix, i, ".bivar.txt")
  f.h2 <- fread(f.h2)
  f.local.gc <- fread(f.local.gc)
  h2[[i]] <- f.h2
  local.gc[[i]] <- f.local.gc
}

h2 <- bind_rows(h2)
local.gc <- bind_rows(local.gc)

# Add adjusted p-values
h2 %<>% arrange(p); local.gc %<>% arrange(p)
h2$bonferroni <- p.adjust(h2$p, method = "bonferroni")
h2$fdr <- p.adjust(h2$p, method = "BH")
local.gc$bonferroni <- p.adjust(local.gc$p, method = "bonferroni")
local.gc$fdr <- p.adjust(local.gc$p, method = "BH")

# Add locus information
h2 %<>% mutate(chrpos = paste0("chr", chr, ":", start, "-", stop))
local.gc %<>% mutate(chrpos = paste0("chr", chr, ":", start, "-", stop)) 

# Save
setwd(output_path)
h2.dir = paste0(output_path, "/", output.filename.prefix, ".local_h2.csv")
local.gc.dir = paste0(output_path, "/", output.filename.prefix, ".local_gc.csv")
if(!file.exists(h2.dir)){
  write.csv(h2, file = h2.dir, row.names = F)
}
if(!file.exists(local.gc.dir)){
  write.csv(local.gc, file = local.gc.dir, row.names = F)
}

#--------------------------------
local.gc %>% filter(fdr < 0.05) %>% distinct(locus) %>% nrow() #23
local.gc %>% filter(fdr < 0.01) %>% distinct(locus) %>% nrow() #10
local.gc %>% filter(bonferroni < 0.05) %>% distinct(locus) %>% nrow() #7

#-------------1. Use bonferroni < 0.05 for plotting: phenotypes as columns, loci as rows ----
sig = "bonferroni05"
local.gc %>% filter(bonferroni < 0.05) -> local.gc.bonferroni05
# Get column names
local.gc %>% filter(bonferroni < 0.05) %>% 
  distinct(locus, .keep_all = T)  %>% 
  select(locus, chrpos) %>% arrange(locus) -> bonferroni05
bonferroni05_col = bonferroni05$locus
# Get row names
bonferroni05_row = c("IBD_eas", "CD_eas", "UC_eas")
# Create empty data frames
gccor = data.frame(matrix(nrow = length(bonferroni05_row), ncol = length(bonferroni05_col)))
pcor = data.frame(matrix(nrow = length(bonferroni05_row), ncol = length(bonferroni05_col)))
bonferronicor = data.frame(matrix(nrow = length(bonferroni05_row), ncol = length(bonferroni05_col)))

rownames(gccor) <- bonferroni05_row; colnames(gccor) <- bonferroni05_col
rownames(pcor) <- bonferroni05_row; colnames(pcor) <- bonferroni05_col
rownames(bonferronicor) <- bonferroni05_row; colnames(bonferronicor) <- bonferroni05_col


for (i in 1:length(bonferroni05_row)) {
  for (m in 1:length(bonferroni05_col)) {
    phenname <- bonferroni05_row[i]  
    locusname <- bonferroni05_col[m]
    if(length(local.gc[locus == locusname & phen1 == phenname, rho]) > 0){
      gccor[i, m] = local.gc[locus == locusname & phen1 == phenname, rho]
      pcor[i, m] = local.gc[locus == locusname & phen1 == phenname, p]
      bonferronicor[i, m] = local.gc[locus == locusname & phen1 == phenname, bonferroni]
    }else{
      gccor[i, m] = 0
      pcor[i, m] = 1
      bonferronicor[i, m] = 1
      print(paste0(i, ".", m, ":", "locus did not calculate gc"))
    }
    print(paste0(i, ".", m))
  }
}
### Replace NAs, NAs occur when rg calculation is outbound
for (i in 1:length(bonferroni05_row)) {
  for (m in 1:length(bonferroni05_col)) {
    phenname <- as.character(bonferroni05_row[i])
    locusname <- bonferroni05_col[m]
    gccor[i, m]
    pcor[i, m]
    if(is.na(gccor[i, m]) == T){
      gccor[i, m] = 0
      pcor[i, m] = 1
      #bonferronicor[locusname, phenname] = 1
      #boncor[locusname, phenname] = 1
      print(paste0(i, ".", m, " is NA"))
    }
    print(paste0(i, ".", m))
  }
}

# fdr/bon matrix processing, extract x and y information for q < 0.01 in fdr matrix
local.gc.bonferroni05 %>% 
  select(locus, phen1) %>%
  mutate(x3 = locus) %>% 
  mutate(x4 = phen1) -> cor.sig
cor.sig$locus <- as.character(cor.sig$locus)
cor.sig$x3 <- as.character(cor.sig$x3)

#--------------------------------- Plotting -----------------------------------#
setwd(output_path)
a <- as.matrix(gccor)
b <- as.matrix(pcor)

# Change row and column names
colnames(a) <- bonferroni05$chrpos
colnames(b) <- bonferroni05$chrpos
rownames(a) <- c("IBD-SLE_CN", "CD-SLE_CN", "UC-SLE_CN")
rownames(b) <- c("IBD-SLE_CN", "CD-SLE_CN", "UC-SLE_CN")

# fdr/bon matrix processing, extract x and y information for q < 0.01 in fdr matrix
local.gc.bonferroni05 %>%
  separate(phen1, into = c("phena", "phenb")) %>% 
  mutate(phen1 = paste(phena, phen2, sep = "-")) %>% 
  select(chrpos, phen1) %>%
  mutate(x3 = chrpos) %>% 
  mutate(x4 = phen1) -> cor.sig

cor.sig$chrpos <- as.character(cor.sig$chrpos)
cor.sig$x3 <- as.character(cor.sig$x3)

##------------------------------ Color Adjustment: Reduce depth of red and blue -----------------#
up_color <- "#B2182B"
median_color <- "white"
bottom_color <- "#1F63A8"
up_palette <- colorRampPalette(c(up_color, median_color))
bottom_palette <- colorRampPalette(c(median_color, bottom_color))
# Set number of colors in palette
n_colors <- 100
up_palette_colors <- up_palette(n_colors)
bottom_palette_colors <- bottom_palette(n_colors)
palette_colors <- c(up_palette_colors, bottom_palette_colors)
#------------------------------ Generate Plot ------------------------------------#
pdf(paste0(output.filename.prefix, ".localgc.", sig, ".pdf"),
    width = 10, 
    height = 10)
corrplot(a,
         method = 'circle', #'circle','square',"ellipse"
         #addCoef.col = 'black', number.cex = 0.27, #### Color and size of numbers in boxes
         col = palette_colors, #COL2('RdBu', 200),
         #is.corr = F,
         tl.cex = 0.7, tl.srt = 45, cl.cex = 0.7,
         tl.col = 'black',
         mar = c(1, 1, 1, 1),
         p.mat = b,
         sig.level = 0.05,
         insig = 'blank'
)$corrPos -> p1
corrplot(a,
         method = 'circle', #'circle','square',"ellipse"
         #addCoef.col = 'black', number.cex = 0.27, #### Color and size of numbers in boxes
         col = palette_colors, #COL2('RdBu', 200),
         #is.corr = F,
         tl.cex = 0.7, tl.srt = 45, cl.cex = 0.7,
         tl.col = 'black',
         mar = c(1, 1, 1, 1),
         p.mat = b,
         sig.level = 0.05,
         insig = 'blank'
) %>% corrRect(namesMat = cor.sig, col = "#049352FF")
p1$corr = round(p1$corr, 2)

for (i in 1:nrow(p1)) {
  gcvalue = p1[i, "corr"]
  pvalue = p1[i, "p.value"]
  if(gcvalue == 0 & pvalue == 1){
    p1[i, "corr"] = ""
    gcvalue = p1[i, "corr"] = ""
  }
  print(i)
}
text(p1$x, p1$y, p1$corr, cex = 0.6)
dev.off()