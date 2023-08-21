# @Project: DOM-in-MicroEco
# @File: keyMole.R
# @Author: Kai Ma
# Date: 29/04/2022


keyMole <- function(dom_t, prefix) {
  require(tidyfst)
  
  pairRI <- function(x1, x2) {
    p0 <- (x1 - x2)/x2*100
    ifelse(is.infinite(p0)|is.na(p0), 0, p0)
  }
  
  site_set <- c('SY', 'BH', 'ZH', 'ST', 'XM', 'WZ', 'NB', 'YC', 'LYG', 'QD', 'DY', 'DD')
  col1 <- 1; pairSite_RI <- matrix(NA, nrow(dom_t), length(site_set)*(length(site_set) - 1) / 2); cnames1 <- rep(NA, length(site_set)*(length(site_set) - 1)/2)
  
  for (s1 in 1:(length(site_set) - 1)) {
    for (s2 in (s1 + 1):length(site_set)) {
      cat(site_set[s1], 'to', site_set[s2], '\n')
      site1_data <- dom_t[, grep(site_set[s1], colnames(dom_t))]
      site2_data <- dom_t[, grep(site_set[s2], colnames(dom_t))]
      
      site_RI <- matrix(NA, nrow(dom_t), ncol(site1_data)*ncol(site2_data) )
      col2 <- 1
      for (i in 1:ncol(site1_data)) {
        for (j in 1:ncol(site2_data)) {
          site_RI[, col2] <- pairRI(site1_data[, i], site2_data[, j])
          col2 <- col2 + 1
        }
      }
      cnames1[col1] <- paste0(site_set[s1], 'to', site_set[s2])
      pairSite_RI[, col1] <- apply(site_RI, 1, median)
      col1 <- col1 + 1
    }
  }
  colnames(pairSite_RI) <- cnames1
  rownames(pairSite_RI) <- rownames(dom_t)
  
  pairSite_RI <- rn_col(as.data.frame(pairSite_RI), var = 'mf')
  rank_df <- apply(pairSite_RI[, -1], 2, function(x) match(x, sort(x, decreasing = T))) 
  
  rank_score_inverse <- 1 - normScale(rowSums(rank_df[, -1]), floor = 0, ceiling = 1, intFmt = FALSE)  
  occu_score_inverse <- normScale(apply(pairSite_RI[, -1], 1, function(x) sum(x > 0)/ncol(pairSite_RI[, -1])), floor = 0, ceiling = 1, intFmt = FALSE)
  LatiAbun_inverse <- pairSite_RI[order(rank_score_inverse+occu_score_inverse, decreasing = TRUE)]
  # LatiAbun_inverse <- arrange_dt(pairSite_RI, -(rank_score_inverse+occu_score_inverse)) 
  
  rank_score_proportional <- normScale(rowSums(rank_df[, -1]), floor = 0, ceiling = 1, intFmt = FALSE) 
  occu_score_proportional <- normScale(apply(pairSite_RI[, -1], 1, function(x) sum(x < 0)/ncol(pairSite_RI[, -1])), floor = 0, ceiling = 1, intFmt = FALSE)
  LatiAbun_proportional <- pairSite_RI[order(rank_score_proportional+occu_score_proportional, decreasing = TRUE)]
  # LatiAbun_proportional <- arrange_dt(pairSite_RI, -(rank_score_proportional+occu_score_proportional))
  
  i <- 1
  while (length(intersect(LatiAbun_inverse[1:i, ]$mf, LatiAbun_proportional[1:i, ]$mf)) == 0) {
    i <- i + 1
    print(i)
  }
  inter_num <- length(intersect(LatiAbun_inverse[1:i, ]$mf, LatiAbun_proportional[1:i, ]$mf))

  fwrite(LatiAbun_inverse[1:(i-inter_num)], paste0(prefix, '-Lati-Abun-inverse.csv'))
  fwrite(LatiAbun_proportional[1:(i-inter_num)], paste0(prefix, '-Lati-Abun-proportional.csv'))
  return(rbind(LatiAbun_inverse[1:(i-inter_num)], LatiAbun_proportional[1:(i-inter_num)]))
}
