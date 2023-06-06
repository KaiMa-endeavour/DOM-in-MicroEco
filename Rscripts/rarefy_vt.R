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
