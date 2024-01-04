# @Project: DOM-in-MicroEco
# @File: DDR.R
# @Author: Kai Ma
# Date: 11/04/2022

library(tidyfst)

mole_traits <- fread('./molecular_traits_by_pool.csv')[, c(1:5, 8, 9, 11:14, 18:21, 31, 32, 34)]
dom <- fread('./DOM_scaled_rarefied_int_tab.csv')
mole_traits2 <- mole_traits[mole_traits$`Molecular formula` %in% colnames(dom)[-1], ]

DBE_range <- c(3:9)

df <- data.frame(DBE = DBE_range)
df$labile <- sapply(DBE_range, function(cc) {
  DBEx <- mole_traits2[mole_traits2$DBE <= cc]
  nrow(DBEx[DBEx$H_C >= 1.5])/nrow(DBEx)
})
df$recalcitrant <- sapply(DBE_range, function(cc) {
  DBEd <- mole_traits2[mole_traits2$DBE > cc]
  nrow(DBEd[DBEd$H_C < 1.5])/nrow(DBEd)
})
df$score <- rowSums(df[, -1])
df$DBE[which.max(df$score)]
