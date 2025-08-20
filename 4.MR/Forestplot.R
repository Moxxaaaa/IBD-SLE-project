library(grid)
library(forestploter)
library(tidyverse)
library(magrittr)
library(readxl)
library(eoffice)
#------------------------ Set Parameters --------------#
MR = "./202407/initial_results/1014SLE_IBD_EAS/MR/MR_forest.xlsx"
Sheet = "Sheet1"
out_path = "./202407/initial_results/1014SLE_IBD_EAS/MR/forest_plots"
phe_c = c("IBD-SLE", "CD-SLE", "UC-SLE")
#---------------------- Separator Line ---------------------#
setwd(out_path)
mr <- read_xlsx(MR, sheet = Sheet)
mr$" " = "                                                "
tm <- forest_theme(base_size = 10,
                   refline_lty = "solid",
                   ci_pch = c(15),
                   ci_col = c("#377eb8"),
                   footnote_col = "blue",
                   legend_name = " ",
                   legend_value = c("IBDs on SLE", "SLE on IBDs"),
                   core = list(bg_params = list(fill = c("white"))))
for (i in 1) {
  phe = phe_c[i]
  mr %>% filter(Trait_pair == phe) -> mr1
  p <- forest(mr1[c(2, 13)],
              est = list(mr1$or, mr1$or1),
              lower = list(mr1$`or_low`, mr1$`or_low1`),
              upper = list(mr1$`or_up`, mr1$`or_up1`),
              ci_column = 2,
              ref_line = 1, # Position of null line
              xlim = c(0.8, 1.4), # Set x-axis limits
              xticks = c(0.8, 1, 1.4),
              ticks_at = c(0.8, 1, 1.2, 1.4), # Set x-axis ticks
              nudge_y = 0.2,
              xlab = "Casual effect (OR)",
              theme = tm)
  topptx(p, filename = paste0(phe, ".pptx"),
         width = 10, height = 10)
}

for (i in 2) {
  phe = phe_c[i]
  mr %>% filter(Trait_pair == phe) -> mr1
  p <- forest(mr1[c(2, 13)],
              est = list(mr1$or, mr1$or1),
              lower = list(mr1$`or_low`, mr1$`or_low1`),
              upper = list(mr1$`or_up`, mr1$`or_up1`),
              ci_column = 2,
              ref_line = 1,
              xlim = c(0.8, 1.4), # Set x-axis limits
              xticks = c(0.8, 1, 1.4),
              ticks_at = c(0.8, 1, 1.2, 1.4), # Set x-axis ticks
              nudge_y = 0.2,
              xlab = "Casual effect (OR)",
              theme = tm)
  topptx(p, filename = paste0(phe, ".pptx"),
         width = 10, height = 10)
}

for (i in 3) {
  phe = phe_c[i]
  mr %>% filter(Trait_pair == phe) -> mr1
  p <- forest(mr1[c(2, 13)],
              est = list(mr1$or, mr1$or1),
              lower = list(mr1$`or_low`, mr1$`or_low1`),
              upper = list(mr1$`or_up`, mr1$`or_up1`),
              ci_column = 2,
              ref_line = 1,
              xlim = c(0.8, 1.4), # Set x-axis limits
              xticks = c(0.8, 1, 1.4),
              ticks_at = c(0.8, 1, 1.2, 1.4), # Set x-axis ticks
              nudge_y = 0.2,
              xlab = "Casual effect (OR)",
              theme = tm)
  topptx(p, filename = paste0(phe, ".pptx"),
         width = 10, height = 10)
}