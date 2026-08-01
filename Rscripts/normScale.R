# @Project: DOM-in-MicroEco
# @File: normScale.R
# @Author: Kai Ma
# Date: 16/08/2022

# This script uses normScale() and min_error() from the MicroEcoTk R package.
#   install.packages("remotes")
#   remotes::install_github("KaiMa-endeavour/MicroEcoTk")
library(MicroEcoTk)
library(tidyfst)

# --- Project-specific analysis (uses MicroEcoTk::normScale & MicroEcoTk::min_error) ---

DOM_int_tab <- fread('./DOM_int_tab.csv')
error_bray <- min_error(dom = DOM_int_tab[, -1], L_margin = 8e-5, R_margin = 12e-5, step = 1e-6, floor = 200)
DOM_scaled_int_tab <- as.data.table(cbind(Samples = DOM_int_tab$Samples, t(apply(DOM_int_tab[, -1], 1, normScale, floor = 200, multiple = 9.8e-5, intFmt = TRUE))))
