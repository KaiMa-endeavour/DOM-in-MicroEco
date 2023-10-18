# @Project: DOM-in-MicroEco
# @File: keyMole.R
# @Author: Kai Ma
# Date: 29/04/2022


keyMole <- function(dom_t, temporal_or_spatial = 'spatial', prefix) {
  stopifnot(temporal_or_spatial %in% c('temporal', 'spatial'))
  require(tidyfst)
  
  pairRI <- function(x1, x2) {
    p0 <- (x1 - x2)/x2*100
    ifelse(is.infinite(p0)|is.na(p0), 0, p0)
  }
  if (temporal_or_spatial == 'spatial') {
    set <- c('SY', 'BH', 'ZH', 'ST', 'XM', 'WZ', 'NB', 'YC', 'LYG', 'QD', 'DY', 'DD')
  }
  if (temporal_or_spatial == 'temporal') {
    set <- c('Aug', 'Oct', 'Dec', 'Feb', 'Apr', 'Jun')
  }
  
  col1 <- 1; dataPair_RI <- matrix(NA, nrow(dom_t), length(set)*(length(set) - 1) / 2); cnames1 <- rep(NA, length(set)*(length(set) - 1)/2)
  
  for (s1 in 1:(length(set) - 1)) {
    for (s2 in (s1 + 1):length(set)) {
      cat(set[s1], 'to', set[s2], '\n')
      sub_data1 <- dom_t[, grep(set[s1], colnames(dom_t))]
      sub_data2 <- dom_t[, grep(set[s2], colnames(dom_t))]
      
      site_RI <- matrix(NA, nrow(dom_t), ncol(sub_data1)*ncol(sub_data2) )
      col2 <- 1
      for (i in 1:ncol(sub_data1)) {
        for (j in 1:ncol(sub_data2)) {
          site_RI[, col2] <- pairRI(sub_data1[, i], sub_data2[, j])
          col2 <- col2 + 1
        }
      }
      cnames1[col1] <- paste0(set[s1], '_2_', set[s2])
      dataPair_RI[, col1] <- apply(site_RI, 1, median)
      col1 <- col1 + 1
    }
  }
  colnames(dataPair_RI) <- cnames1
  rownames(dataPair_RI) <- rownames(dom_t)
  
  dataPair_RI <- rn_col(as.data.frame(dataPair_RI), var = 'mf') #[rowSums(dataPair_RI) > 0, ]
  rank_df <- apply(dataPair_RI[, -1], 2, function(x) match(x, sort(x, decreasing = T)))
  
  rank_score_negative <- 1 - normScale(rowSums(rank_df[, -1]), floor = 0, ceiling = 1, intFmt = FALSE)
  occu_score_negative <- normScale(apply(dataPair_RI[, -1], 1, function(x) sum(x > 0)/ncol(dataPair_RI[, -1])), floor = 0, ceiling = 1, intFmt = FALSE)
  negative_dt <- dataPair_RI[order(rank_score_negative+occu_score_negative, decreasing = TRUE)]
  # negative_dt <- arrange_dt(dataPair_RI, -(rank_score_negative+occu_score_negative)) # [1:10, ]
  
  rank_score_positive <- normScale(rowSums(rank_df[, -1]), floor = 0, ceiling = 1, intFmt = FALSE) 
  occu_score_positive <- normScale(apply(dataPair_RI[, -1], 1, function(x) sum(x < 0)/ncol(dataPair_RI[, -1])), floor = 0, ceiling = 1, intFmt = FALSE)
  positive_dt <- dataPair_RI[order(rank_score_positive+occu_score_positive, decreasing = TRUE)]
  # positive_dt <- arrange_dt(dataPair_RI, -(rank_score_positive+occu_score_positive))
  
  i <- 1
  while (length(intersect(negative_dt[1:i, ]$mf, positive_dt[1:i, ]$mf)) == 0) {
    i <- i + 1
    print(i)
  }
  
  fwrite(negative_dt[1:(i-1)], paste(prefix, temporal_or_spatial, 'negative.csv', sep = '-'))
  fwrite(positive_dt[1:(i-1)], paste(prefix, temporal_or_spatial, 'positive.csv', sep = '-'))
  return(rbind(negative_dt[1:(i-1)], positive_dt[1:(i-1)]))
}
