# @Project: DOM-in-MicroEco
# @File: normScale.R
# @Author: Kai Ma
# Date: 16/08/2022

library(tidyfst)
# Scales each element in the vector by the same proportion.
normScale <- function(v, floor = 100, ceiling = NULL, multiple = NULL, intFmt = FALSE) {
  vmax <- max(v[v > 0])
  vmin <- min(v[v > 0])
  if (is.null(ceiling)) {
    if (is.null(multiple)) {
      xx <- (sqrt(max(v)) - floor) * (v[v > 0] - vmin) / (vmax - vmin) + floor
    } else {
      xx <- (max(v) * multiple - floor) * (v[v > 0] - vmin) / (vmax - vmin) + floor
    }
  } else {
    xx <- (ceiling - floor) * (v[v > 0] - vmin) / (vmax - vmin) + floor
  }
  v[v > 0] <- xx
  if (!intFmt) v else round(v, 0)
}
# Find the multiple with the smallest error in the interval set with L_margin and R_margin.
min_error <- function(dom, L_margin, R_margin, step, method = "bray", floor = 100, ceiling = NULL) {
  beta_dom <- vegan::vegdist(dom, method = method)
  error <- sapply(seq(L_margin, R_margin, step), function(k) {
    dom_scaling <- t(apply(dom, 1, normScale, floor = floor, ceiling = ceiling, multiple = k))
    beta <- vegan::vegdist(dom_scaling, method = method)
    abs(mean(beta_dom) - mean(beta))
  }) %>% as.vector()
  dt <- data.table(multiple = seq(L_margin, R_margin, step), error = error)
  res <- seq(L_margin, R_margin, step)[which.min(dt$error)]
  list(dt, res)
}
DOM_int_tab <- fread('./DOM_int_tab.csv')
error_bray <- min_error(dom = DOM_int_tab[, -1], L_margin = 8e-5, R_margin = 12e-5, step = 1e-6, floor = 200)
DOM_scaled_int_tab <- as.data.table(cbind(Samples = DOM_int_tab$Samples, t(apply(DOM_int_tab[, -1], 1, normScale, floor = 200, multiple = 9.8e-5, intFmt = TRUE))))
