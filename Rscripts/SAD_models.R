# @Project: DOM-in-MicroEco
# @File: SAD_models.R
# @Author: Kai Ma
# Date: 03/12/2020


Rsquare <- function(obs, pred) {
  1 - (sum((log10(obs) - log10(pred))^2) / sum((log10(obs) - mean(log10(obs)))^2))
}

se <- function(x, na.rm = F) {
  x <- as.numeric(x)
  stats::sd(x, na.rm) / sqrt(length(x) - sum(is.na(x)))
}

rank_abun <- function(abun_vet) {
  rep(names(table(log10(abun_vet))), table(log10(abun_vet))) %>%
    rev() %>%
    as.numeric()
}

fitsad_zipf <- function(x, N, trunc, start.value, upper = 20, ...) {
  dots <- list(...)
  if (any(x <= 0)) stop("All x must be positive integers")
  if (!missing(trunc)) {
    if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(N)) {
    N <- sum(x)
  }
  if (missing(start.value)) {
    p <- x / sum(x)
    lzipf <- function(s, N) -s * log(1:N) - log(sum(1 / (1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    sss <- opt$minimum
  } else {
    sss <- start.value
  }
  if (missing(trunc)) {
    LL <- function(N, s) -sum(dzipf(x, N = N, s = s, log = TRUE))
  } else {
    LL <- function(N, s) -sum(dtrunc("zipf", x = x, coef = list(N = N, s = s), trunc = trunc, log = TRUE))
  }
  result <- do.call("mle2", c(list(LL, start = list(s = sss), data = list(x = x), fixed = list(N = N), method = "Brent", lower = 0, upper = upper), dots))
  if (abs(as.numeric(result@coef) - upper) < 0.001) {
    warning("mle equal to upper bound provided. \n Try increase value for the 'upper' argument")
  }
  new("fitsad", result, sad = "zipf", distr = "discrete", trunc = ifelse(missing(trunc), NaN, trunc))
}

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
    # pred_zi <- VGAM::rzipf(length(zi@data$x), zi@fullcoef[1], zi@fullcoef[2]) %>% rarefy_vt(zi@fullcoef[1]) %>% sort(decreasing = T)
    # pred_zi[pred_zi == 0] <- 1
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
