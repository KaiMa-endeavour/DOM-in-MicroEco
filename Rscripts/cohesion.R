# @Project: DOM-in-MicroEco
# @File: cohesion.R
# @Author: Kai Ma
# Date: 23/05/2023

rel_ab <- function(comm, MARGIN = 1) {
  as.data.frame(t(apply(comm, MARGIN = MARGIN, function(x) {
    x / sum(x)
  })))
}

cohesion <- function(comm, method = "spearman", true_cor = TRUE, nworker = 4, reps = 999) {
  relab <- rel_ab(comm)
  obs_cor <- stats::cor(relab, method = method)
  if (true_cor) {
    require(parallel)
    require(foreach)
    require(doParallel)

    cl <- makeCluster(nworker, type = "PSOCK")
    registerDoParallel(cl)
    perm_fun <- function(taxon_id) {
      null_taxon <- sapply(1:reps, function(i) {
        shuffle_relab <- sapply(relab, sample)
        shuffle_relab[, taxon_id] <- relab[, taxon_id]
        null_cor <- stats::cor(shuffle_relab, method = method)
        null_cor[, taxon_id]
      })
      expect_cor <- apply(null_taxon, 1, stats::median)
    }
    expected_cor <- foreach(j = 1:ncol(comm), .combine = "cbind") %dopar% perm_fun(j)
    stopCluster(cl)
    diff_cor <- obs_cor - expected_cor
  } else {
    diff_cor <- obs_cor
    diag(diff_cor) <- 0
  }
  # apply(diff_cor, 1, function(x) {(x - mean(x))/sd(x)}) %>% rowSums()
  posi <- apply(diff_cor, 2, function(x) ifelse(length(x[x > 0]) != 0, mean(x[x > 0]), 0))
  nega <- apply(diff_cor, 2, function(x) ifelse(length(x[x < 0]) != 0, mean(x[x < 0]), 0))
  posi_co <- as.matrix(relab) %*% posi
  nega_co <- as.matrix(relab) %*% nega
  data.frame(Posi_cohesion = posi_co, Nega_cohesion = nega_co)
}