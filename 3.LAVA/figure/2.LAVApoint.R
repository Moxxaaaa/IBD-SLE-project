library(ggplot2)
library(readxl)
library(magrittr)
library(tidyverse)
#-------------------------------------- Set Parameters ---------------------------------------#
f1 = "./results/SLE_IBD_EAS/LAVA/SLE_IBD_eas.local_gc.csv"
filename = "SLE-IBD"
output.dir = "./results/SLE_IBD_EAS/scatter_plots"
f1 <- read.csv(f1)
f1 %<>% 
  filter(is.na(rho) == F) %>% 
  mutate(MlogP = -log10(p))

f1 %>% filter(p < 0.05) -> f_nom
setwd(output.dir)


#-------------------------------------------------- Add Group Information ---------------------------------#
# Divide into three groups
f1$group.color <- NA
for (i in 1:nrow(f1)) {
  if(f1[i, "p"] >= 0.05){f1[i, "group.color"] <- 1}
  if(f1[i, "p"] < 0.05 & f1[i, "rho"] > 0){
    f1[i, "group.color"] <- 2}
  if(f1[i, "p"] < 0.05 & f1[i, "rho"] < 0){
    f1[i, "group.color"] <- 3}
}

# Divide into three groups
f1$group <- NA
for (i in 1:nrow(f1)) {
  if(f1[i, "p"] >= 0.05){f1[i, "group"] <- 1}
  if(f1[i, "p"] < 0.05 & f1[i, "rho"] > 0){
    f1[i, "group"] <- 2}
  if(f1[i, "p"] < 0.05 & f1[i, "rho"] < 0){
    f1[i, "group"] <- 3}
}

# Group points by shape based on bonferroni threshold
f1$group.shape <- NA
for (i in 1:nrow(f1)) {
  if(f1[i, "bonferroni"] >= 0.05){
    f1[i, "group.shape"] <- 1}
  if(f1[i, "bonferroni"] < 0.05){
    f1[i, "group.shape"] <- 2}
}

#-------------------------------------------- Plotting ----------------------------------------------#

#------------------------------------ Scatter plot for P < 0.05 -------------------------------------
ggplot(f_nom, aes(x = locus, y = rho)) + 
  geom_point(size = 0.5) + 
  labs(x = "Locus", y = "Local Genetic Correlation")

ggsave(paste0(filename, ".p05.png"), 
       path = output.dir,
       width = 8,
       height = 5,
       dpi = 300)

ggsave(paste0(filename, ".p05.pdf"),
       path = output.dir,
       width = 8,
       height = 5,
       dpi = 300)

#------------------------------------------- All points + color for significant points ----------------------------------
ggplot(f1, aes(x = locus, y = rho, color = factor(group), shape = factor(group))) +
  geom_point(size = 1) +
  scale_color_manual(values = c("gray", "#477FB7", "#C04151"),
                     labels = c("P ≥ 0.05", "P < 0.05 & rg > 0", "P < 0.05 & rg < 0")) +
  scale_shape_manual(values = c(16, 16, 16),
                     labels = c("P ≥ 0.05", "P < 0.05 & rg > 0", "P < 0.05 & rg < 0")) +
  labs(color = "", shape = "") +
  labs(x = "Locus", y = "Local Genetic Correlation")

ggsave(paste0(filename, ".all.png"),
       path = output.dir,
       width = 10,
       height = 5,
       dpi = 300)

ggsave(paste0(filename, ".all.pdf"),
       path = output.dir,
       width = 10,
       height = 5,
       dpi = 300)

#--------------------------------------- Volcano plot of P and rg -----------------------------#####
ggplot(f1, aes(x = rho, y = MlogP, color = factor(group.color), shape = factor(group.shape))) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("gray", "#477FB7", "#C04151"),
                     labels = c("P ≥ 0.05", "rg > 0", "rg < 0")) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("bonferroni < 0.05", "")) +
  labs(color = NULL, shape = NULL) +
  labs(x = "Local Genetic Correlation", y = "-Log10P") +
  guides(shape = guide_legend(override.aes = list(shape = c(17, NA)))) +
  theme(legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.key = element_blank()) +
  theme(legend.spacing.y = unit(0, "cm"))

ggsave(paste0(filename, ".volcano.png"),
       path = output.dir,
       width = 6,
       height = 6,
       dpi = 300)

ggsave(paste0(filename, ".volcano.pdf"),
       path = output.dir,
       width = 6,
       height = 6,
       dpi = 300)

#------------------------------------ Single phenotype volcano plot --------------------------------------
f_c <- c("IBD_eas", "CD_eas", "UC_eas")
for (i in 1:length(f_c)) {
  f.name = f_c[i]
  f1 %>% filter(phen1 == f.name) -> temp
  
  ggplot(temp, aes(x = rho, y = MlogP, color = factor(group.color), 
                   #shape = factor(group.shape)
  )) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("gray", "#477FB7", "#C04151"),
                       labels = c("P ≥ 0.05", "rg > 0", "rg < 0")) +
    #scale_shape_manual(values = c(16, 17), labels = c("bonferroni < 0.05", "")) +
    labs(color = NULL, shape = NULL) +
    labs(x = "Local Genetic Correlation", y = "-Log10P") +
    guides(shape = guide_legend(override.aes = list(shape = c(17, NA)))) +
    theme(legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          legend.key = element_blank()) +
    theme(legend.spacing.y = unit(0, "cm"))
  
  ggsave(paste0(f.name, ".volcano.png"),
         path = output.dir,
         width = 5,
         height = 6,
         dpi = 300)
  
  ggsave(paste0(f.name, ".volcano.pdf"),
         path = output.dir,
         width = 5,
         height = 6,
         dpi = 300)
  print(i)
}