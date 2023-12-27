library(tidyfst)

dominance <- function(a, rel_D = F) {
  if (!rel_D) {
    return(max(a) / sum(a))
  } else {
    return(max(a / sum(a)))
  }
}

rarity <- function(a) {
  as <- unlist(a[a > 0])
  mean(((as - mean(as)) / stats::sd(as))^3)
}

alphaDiversity <- function(comm, sample_names, methods = c("Richness", "Chao1", "ACE", "Shannon", "Simpson", "Pielou", "goods_coverage"), suffix = NULL, tree = NULL, base = exp(1)) {
  stopifnot(methods %in% c("Richness", "Chao1", "ACE", "Shannon", "Simpson", "Pielou", "goods_coverage"))
  Richness <- Chao1 <- ACE <- Shannon <- Simpson <- Pielou <- goods_coverage <- NA
  if ("Richness" %in% methods) {
    Richness <- vegan::estimateR(comm)[1, ]
  }
  if ("Chao1" %in% methods) {
    Chao1 <- vegan::estimateR(comm)[2, ]
  }
  if ("ACE" %in% methods) {
    ACE <- vegan::estimateR(comm)[4, ]
  }
  if ("Shannon" %in% methods) {
    Shannon <- vegan::diversity(comm, index = "shannon", base = base)
  }
  if ("Simpson" %in% methods) {
    Simpson <- vegan::diversity(comm, index = "simpson") # Gini-Simpson
  }
  if ("Pielou" %in% methods) {
    Shannon <- vegan::diversity(comm, index = "shannon", base = base)
    Richness <- vegan::estimateR(comm)[1, ]
    Pielou <- Shannon / log(Richness, base)
  }
  if ("goods_coverage" %in% methods) {
    goods_coverage <- 1 - rowSums(comm == 1) / rowSums(comm)
  }
  
  res <- data.frame(Samples = sample_names, Richness, Chao1, ACE, Shannon, Simpson, Pielou, goods_coverage)[, c(1, 1 + match(methods, c("Richness", "Chao1", "ACE", "Shannon", "Simpson", "Pielou", "goods_coverage")))]
  
  if (!is.null(tree)) {
    PD_whole_tree <- picante::pd(comm, tree, include.root = FALSE)[1]
    res <- cbind(res, PD = PD_whole_tree)
  }
  if (!is.null(suffix)) {
    colnames(res) <- c("Samples", paste(methods, suffix, sep = "_"))
  }
  as.data.table(res)
}

data <- fread('./DOM_scaled_rarefied_int_tab.csv')
alpha_dom <- lapply(c(2**seq(5, 100)[2**seq(5, 100) < max(rowSums(data[, -1]))], max(rowSums(data[, -1]))), function(p) {
  data_sub <- data.table(Samples = data$Samples, rarefy(data[, -1], depth = p))
  alphaDiversity(data_sub[, -1], data_sub$Samples, method = c('Richness', 'Pielou')) %>% 
    cbind(Dominance = apply(data_sub[, -1], 1, function(a) max(a)), 
          Rarity = apply(data_sub[, -1], 1, rarity), 
          Sample_size = p, 
          S = apply(data_sub[, -1], 1, function(a) length(a[a>0])))
}) %>% rbindlist() %>% mutate_dt(group = 'DOM')

data <- fread('./ASVstable_rarefied_Bacteria.csv')
alpha_bac <- lapply(c(2**seq(5, 100)[2**seq(5, 100) < max(rowSums(data[, -1]))], max(rowSums(data[, -1]))), function(p) {
  data_sub <- data.table(Samples = data$Samples, rarefy(data[, -1], depth = p))
  alphaDiversity(data_sub[, -1], data_sub$Samples, method = c('Richness', 'Pielou')) %>% 
    cbind(Dominance = apply(data_sub[, -1], 1, function(a) max(a)), 
          Rarity = apply(data_sub[, -1], 1, rarity), 
          Sample_size = p, 
          S = apply(data_sub[, -1], 1, function(a) length(a[a>0])))
}) %>% rbindlist() %>% mutate_dt(group = 'Bacteria')

data <- fread('./ASVstable_rarefied_Fungi.csv')
alpha_fun <- lapply(c(2**seq(5, 100)[2**seq(5, 100) < max(rowSums(data[, -1]))], max(rowSums(data[, -1]))), function(p) {
  data_sub <- data.table(Samples = data$Samples, rarefy(data[, -1], depth = p))
  alphaDiversity(data_sub[, -1], data_sub$Samples, method = c('Richness', 'Pielou')) %>% 
    cbind(Dominance = apply(data_sub[, -1], 1, function(a) max(a)), 
          Rarity = apply(data_sub[, -1], 1, rarity), 
          Sample_size = p, 
          S = apply(data_sub[, -1], 1, function(a) length(a[a>0])))
}) %>% rbindlist() %>% mutate_dt(group = 'Fungi')

alpha <- rbind(alpha_dom, alpha_bac, alpha_fun)
