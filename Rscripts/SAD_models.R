# @Project: DOM-in-MicroEco
# @File: SAD_models.R
# @Author: Kai Ma
# Date: 03/12/2020

# This script uses Rsquare(), se(), rank_abun(), and fitsad_zipf() from the MicroEcoTk R package.
#   install.packages("remotes")
#   remotes::install_github("KaiMa-endeavour/MicroEcoTk")
library(MicroEcoTk)

# --- Project-specific analysis (uses MicroEcoTk::Rsquare, MicroEcoTk::fitsad_zipf, etc.) ---

SAD_models <- function(data) {
  require(sads)
  require(mobsim)
  require(tidyfst)

  xx <- matrix(0, nrow(data), 4)
  rank_sample <- list()
  for (r in 1:nrow(data)) {
    cat(r, '\n')
    obs <- round(data[r, ], 0) %>% as.integer() %>% sort(decreasing = T)
    obs <- obs[obs > 0]

    pl <- fitsad(obs, "poilog")
    pred_pl <- sim_sad(s_pool = length(pl@data$x), n_sim = sum(pl@data$x), sad_type = "poilog", sad_coef = list("mu" = pl@fullcoef[1], "sig" = pl@fullcoef[2]), fix_s_sim = T)
    Rq_pl <- Rsquare(obs, pred_pl)

    ls <- fitsad(obs, "ls")
    pred_ls <- sim_sad(s_pool = NULL, n_sim = sum(ls@data$x), sad_type = "ls", sad_coef = list("N" = ls@fullcoef[1], "alpha" = ls@fullcoef[2]), fix_s_sim = T)
    Rq_ls <- Rsquare(obs, pred_ls)

    bs <- fitsad(obs, "bs")
    pred_bs <- sim_sad(s_pool = NULL, sum(bs@data$x), sad_type = "bs", sad_coef = list("N" = bs@fullcoef[1], "S" = bs@fullcoef[2]), fix_s_sim = T)
    Rq_bs <- Rsquare(obs, pred_bs)

    zi <- fitsad_zipf(obs)
    pred_zi <- VGAM::rzipf(length(zi@data$x), zi@fullcoef[1], zi@fullcoef[2]) %>% sort(decreasing = T)
    prop <- pred_zi/sum(pred_zi); names(prop) <- paste0("species", 1:length(pred_zi))
    pred_zi <- sample(paste0("species", 1:length(pred_zi)), zi@fullcoef[1], replace = T, prob = prop) %>% factor(levels = names(prop)) %>% table() %>% sort(decreasing = T)
    pred_zi[pred_zi == 0] <- 1
    Rq_zi <- Rsquare(obs, pred_zi)

    rank_sample[[r]] <- sapply(list(Observed = obs, Lognormal = pred_pl, `Broken-stick` = pred_bs, `Log-series` = pred_ls, Zipf = pred_zi), rank_abun) %>% mutate_dt(rank = seq(1, length(obs)), reps = r)

    Rsquare <- c(Rq_pl, Rq_ls, Rq_bs, Rq_zi)
    xx[r, ] <- Rsquare
    colnames(xx) <- c('Lognormal', 'Log-series', 'Broken-stick', 'Zipf')
  }
  M <- colMeans(xx)
  SE <- apply(xx, 2, se)
  prob <- (apply(xx, 1, which.max) %>% table())/nrow(data)
  list(Rsq = xx, mean = M, se = SE, prob = prob, rank_sample = rbindlist(rank_sample))
}
