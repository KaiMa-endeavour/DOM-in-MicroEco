# @Project: DOM-in-MicroEco
# @File: RandomSample_alphaDiversity.R
# @Author: Kai Ma
# Date: 22/05/2022

# This script uses dominance(), rarity(), rarefy(), and alpha_diversity() from the MicroEcoTk R package.
#   install.packages("remotes")
#   remotes::install_github("KaiMa-endeavour/MicroEcoTk")
library(MicroEcoTk)
library(tidyfst)

# --- Project-specific analysis (uses MicroEcoTk::rarefy, MicroEcoTk::dominance, MicroEcoTk::rarity, MicroEcoTk::alpha_diversity) ---

data <- fread('./DOM_scaled_rarefied_int_tab.csv')
alpha_dom <- lapply(c(2**seq(5, 100)[2**seq(5, 100) < max(rowSums(data[, -1]))], max(rowSums(data[, -1]))), function(p) {
  data_sub <- data.table(Samples = data$Samples, rarefy(data[, -1], depth = p))
  alpha_diversity(data_sub[, -1], data_sub$Samples, methods = c('Richness', 'Pielou')) %>%
    cbind(Dominance = apply(data_sub[, -1], 1, dominance),
          Rarity = apply(data_sub[, -1], 1, rarity),
          Sample_size = p,
          S = apply(data_sub[, -1], 1, function(a) length(a[a>0])))
}) %>% rbindlist() %>% mutate_dt(group = 'DOM')

data <- fread('./ASVstable_rarefied_Bacteria.csv')
alpha_bac <- lapply(c(2**seq(5, 100)[2**seq(5, 100) < max(rowSums(data[, -1]))], max(rowSums(data[, -1]))), function(p) {
  data_sub <- data.table(Samples = data$Samples, rarefy(data[, -1], depth = p))
  alpha_diversity(data_sub[, -1], data_sub$Samples, methods = c('Richness', 'Pielou')) %>%
    cbind(Dominance = apply(data_sub[, -1], 1, dominance),
          Rarity = apply(data_sub[, -1], 1, rarity),
          Sample_size = p,
          S = apply(data_sub[, -1], 1, function(a) length(a[a>0])))
}) %>% rbindlist() %>% mutate_dt(group = 'Bacteria')

data <- fread('./ASVstable_rarefied_Fungi.csv')
alpha_fun <- lapply(c(2**seq(5, 100)[2**seq(5, 100) < max(rowSums(data[, -1]))], max(rowSums(data[, -1]))), function(p) {
  data_sub <- data.table(Samples = data$Samples, rarefy(data[, -1], depth = p))
  alpha_diversity(data_sub[, -1], data_sub$Samples, methods = c('Richness', 'Pielou')) %>%
    cbind(Dominance = apply(data_sub[, -1], 1, dominance),
          Rarity = apply(data_sub[, -1], 1, rarity),
          Sample_size = p,
          S = apply(data_sub[, -1], 1, function(a) length(a[a>0])))
}) %>% rbindlist() %>% mutate_dt(group = 'Fungi')

alpha <- rbind(alpha_dom, alpha_bac, alpha_fun)
