# --*-- encoding: utf-8 --*--
'''
@File       :    normScale.R
@Author     :    Kai Ma
@Time       :    2022/08/16 14:04:34
'''

normScale <- function(v, floor = 100, ceiling = NULL, multiple = NULL) {
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
  return(round(v, 0))
}

min_error <- function(dom, L_margin, R_margin, step, method = "bray") {
  error <- sapply(seq(L_margin, R_margin, step), function(k) {
    dom_scaling <- t(apply(dom, 1, normScale, multiple = k))
    beta_dom <- vegan::vegdist(dom, method = method)
    beta <- vegan::vegdist(dom_scaling, method = method)
    abs(mean(beta_dom) - mean(beta))
  }) %>% as.vector()
  dt <- data.table(multiple = seq(L_margin, R_margin, step), error = error)
  res <- seq(L_margin, R_margin, step)[which.min(dt$error)]
  list(dt, res)
}
