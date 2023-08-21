# @Project: DOM-in-MicroEco
# @File: normScale.R
# @Author: Kai Ma
# Date: 16/08/2022

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

min_error <- function(dom, L_margin, R_margin, step, method = "bray", floor = 100, ceiling = NULL) {
  error <- sapply(seq(L_margin, R_margin, step), function(k) {
    dom_scaling <- t(apply(dom, 1, normScale, floor = floor, ceiling = ceiling, multiple = k))
    beta_dom <- vegan::vegdist(dom, method = method)
    beta <- vegan::vegdist(dom_scaling, method = method)
    abs(mean(beta_dom) - mean(beta))
  }) %>% as.vector()
  dt <- data.table(multiple = seq(L_margin, R_margin, step), error = error)
  res <- seq(L_margin, R_margin, step)[which.min(dt$error)]
  list(dt, res)
}
