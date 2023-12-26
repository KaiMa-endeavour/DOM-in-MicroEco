# @Project: RandomSampling
# @File: rarefy_vt.R
# @Author: Kai Ma
# Date: 25/11/2020


rarefy_vt <- function(x, depth, prob = NULL, prob_dec = FALSE, replace = FALSE) {
  depth <- round(as.numeric(depth), 0)
  stopifnot(is.numeric(x))
  if (depth > sum(x) & replace == F) {
    replace <- TRUE
    warning("replace = TRUE")
  }
  if (is.null(prob)) {
    prob <- prob
  } else {
    if (prob_dec) {
      prob <- prob[match(x, sort(x, decreasing = T))]
    }
    prob <- rep(prob, x)
  }
  y <- sample(x = rep(1:length(x), x), size = depth, prob = prob, replace = replace)
  y_tab <- table(y)
  z <- numeric(length(x))
  z[as.numeric(names(y_tab))] <- y_tab
  z
}

rarefy <- function(otu_table, MARGIN = 1, depth = min(rowSums(otu_table)), prob = NULL, replace = F, nworker = 1) {
  otu_table <- apply(otu_table, 2, as.numeric)
  if (nworker == 1) {
    df <- t(apply(X = otu_table, MARGIN = MARGIN, FUN = rarefy_vt, depth = depth, prob = prob, replace = replace))
  } else {
    require(parallel)
    cl <- makeCluster(nworker, type = "PSOCK")
    df <- t(parApply(cl, X = otu_table, MARGIN = MARGIN, FUN = rarefy_vt, depth = depth, prob = prob, replace = replace))
    stopCluster(cl)
  }
  rownames(df) <- rownames(otu_table)
  colnames(df) <- colnames(otu_table)
  df <- df[rowSums(df) > 0, ]
  df <- df[, colSums(df) > 0]
  df <- as.data.frame(df)
}
