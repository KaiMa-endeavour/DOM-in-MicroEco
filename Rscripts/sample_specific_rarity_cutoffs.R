# @Project: DOM-in-MicroEco
# @File: sample_specific_rarity_cutoffs.R
# @Author: Xiu Jia
# @Reviser: Kai Ma
# Date: 09/07/2022
# The R function from https://github.com/Jia-Xiu/rare_biosphere_assembly_2020/blob/master/sample_specific_rarity_cutoffs.R was improved in order to divide the occasional taxa.


library(vegan)


if (nrow(comm) >= 2) {
  comm <- apply(comm, 2, as.numeric)
}

Chao <- as.data.frame(t(estimateR(comm)))
Chao$slope <- Chao$S.obs / Chao$S.chao1
cutoffs <- matrix(NA, nrow(comm), 3)
cutoffs[, 2] <- Chao$slope
cutoffs[, 3] <- Chao$S.obs
colnames(cutoffs) <- c("Rarity.cutoffs", "Slopes", "S.obs")
for (i in 1:nrow(comm)) {
  comm_i <- sort(as.numeric(comm[i, ]), decreasing = T)
  comm_i <- comm_i[comm_i > 0]
  slope <- cutoffs[i, 2]
  cutoffs[i, 1] <- max(which(ifelse(comm_i >= (1:length(comm_i)) * slope, T, F)))
}
df <- comm
for (j in 1:nrow(df)) {
  for (i in 1:ncol(df)) {
    if (df[j, i] > cutoffs[j, 1]) {
      df[j, i] <- NA
    }
  }
}
df[is.na(df)] <- 0
if (nrow(comm) >= 2) {
  df <- df[, colSums(df) != 0]
  df_rare <- comm[, colnames(df)]
} else {
  df <- t(as.matrix(df[, df != 0]))
  df_rare <- comm[, colnames(df)]
  df_rare <- t(as.matrix(df_rare))
}

df <- comm
for (j in 1:nrow(df)) {
  for (i in 1:ncol(df)) {
    if (df[j, i] <= cutoffs[j, 1]) {
      df[j, i] <- NA
    }
  }
}
df[is.na(df)] <- 0
if (nrow(comm) >= 2) {
  df <- df[, colSums(df) != 0]
  df_abundant <- comm[, colnames(df)]
} else {
  df <- t(as.matrix(df[, df != 0]))
  df_abundant <- comm[, colnames(df)]
  df_abundant <- t(as.matrix(df_abundant))
}

occasional <- comm[, intersect(colnames(df_abundant), colnames(df_rare))]
if (!is.null(dim(occasional))) {
  if (nrow(occasional) >= 2) {
    occasional <- occasional[, colSums(occasional) > 0]
    occasional <- occasional[rowSums(occasional) > 0, ]
  } else if (nrow(occasional) == 1) {
    occasional <- t(as.matrix(occasional[occasional != 0]))
  }
}

abundant <- df_abundant[, setdiff(colnames(df_abundant), intersect(colnames(df_abundant), colnames(df_rare)))]
if (nrow(df_abundant) >= 2) {
  abundant <- abundant[, colSums(abundant) > 0]
  abundant <- abundant[rowSums(abundant) > 0, ]
} else {
  abundant <- t(as.matrix(abundant[abundant != 0]))
}

rare <- df_rare[, setdiff(colnames(df_rare), intersect(colnames(df_abundant), colnames(df_rare)))]
if (nrow(df_rare) >= 2) {
  rare <- rare[, colSums(rare) > 0]
  rare <- rare[rowSums(rare) > 0, ]
} else {
  rare <- t(as.matrix(rare[rare != 0]))
}

write.csv(abundant, 'abundant taxa.csv')
write.csv(occasional, 'occasional taxa.csv')
write.csv(rare, 'rare taxa.csv')
