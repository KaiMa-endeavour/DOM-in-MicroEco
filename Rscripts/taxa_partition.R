# @Project: DOM-in-MicroEco
# @File: taxa_partition.R
# @Author: Kai Ma
# Date: 14/12/2020

library(tidyfst)

rel_ab <- function(comm, MARGIN = 1) {
  t(apply(comm, MARGIN = MARGIN, function(x) {
    x / sum(x) * 100
  })) %>% as.data.frame()
}

taxa_partition <- function(comm, mode = "rarity", relab_abundant = 0.1, relab_rare = 0.01, freq_abundant = 0.75, freq_rare = 0.5, out_relabma = F) {
  if (!mode %in% c("relab", "rarity", "frequency")) stop('The parameter mode can only be "relab", "rarity" or "frequency"')
  if (nrow(comm) > 1) {
    comm <- apply(comm, 2, as.numeric)
    if (mode == "relab") {
      comm_relab <- rel_ab(comm, MARGIN = 1)
      abundant <- comm[, colSums(comm_relab > relab_abundant) > nrow(comm_relab) * 0.5]
      rare <- comm[, colSums(comm_relab < relab_rare) > nrow(comm_relab) * 0.5]
      occasional <- comm[, setdiff(colnames(comm_relab), c(colnames(abundant), colnames(rare)))]
    }
    if (mode == "rarity") {
      Chao <- as.data.frame(t(estimateR(comm)))
      Chao$slope <- Chao$S.obs / Chao$S.chao1
      cutoffs <- matrix(NA, nrow(comm), 3)
      cutoffs[, 2] <- Chao$slope
      cutoffs[, 3] <- Chao$S.obs
      colnames(cutoffs) <- c("rarity.cutoffs", "Slopes", "S.obs")
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
      df <- df[, colSums(df) != 0]
      df_rare <- comm[, colnames(df)]

      df <- comm
      for (j in 1:nrow(df)) {
        for (i in 1:ncol(df)) {
          if (df[j, i] <= cutoffs[j, 1]) {
            df[j, i] <- NA
          }
        }
      }
      df[is.na(df)] <- 0
      df <- df[, colSums(df) != 0]
      df_abundant <- comm[, colnames(df)]

      occasional <- comm[, intersect(colnames(df_abundant), colnames(df_rare))]
      abundant <- df_abundant[, setdiff(colnames(df_abundant), intersect(colnames(df_abundant), colnames(df_rare)))]
      rare <- df_rare[, setdiff(colnames(df_rare), intersect(colnames(df_abundant), colnames(df_rare)))]
    }
    if (mode == "frequency") {
      abundant <- comm[, colnames(comm)[colSums(comm > 0) / nrow(comm) >= freq_abundant]]
      rare <- comm[, colnames(comm)[colSums(comm > 0) / nrow(comm) < freq_rare]]
      occasional <- comm[, setdiff(colnames(comm), c(colnames(abundant), colnames(rare)))]
    }
  } else {
    stop("Single sample was not suitable for use with these methods.")
  }

  abundant <- abundant[, colSums(abundant) > 0]
  abundant <- abundant[rowSums(abundant) > 0, ]
  occasional <- occasional[, colSums(occasional) > 0]
  occasional <- occasional[rowSums(occasional) > 0, ]
  rare <- rare[, colSums(rare) > 0]
  rare <- rare[rowSums(rare) > 0, ]

  if (out_relabma) {
    abundant <- rel_ab(abundant, MARGIN = 1)
    occasional <- rel_ab(occasional, MARGIN = 1)
    rare <- rel_ab(rare, MARGIN = 1)
  }
  cat("\nAbundant subcommunity contains ", ncol(abundant), " OTUs (ASVs).", "\n")
  cat("Occasional subcommunity contains ", ncol(occasional), " OTUs (ASVs).", "\n")
  cat("Rare subcommunity contains ", ncol(rare), " OTUs (ASVs).", "\n\n")
  list(abundant_subcomm = as.data.table(abundant), occasional_subcomm = as.data.table(occasional), rare_subcomm = as.data.table(rare))
}
