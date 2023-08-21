# @Project: DOM-in-MicroEco
# @File: nullStoc.R
# @Author: Kai Ma
# Date: 23/11/2020


nullModel <- function(comm, null_model, sp_freq, samp_rich, output_nullcomm = F, dist_method = "bray") {
  stopifnot(null_model %in% c("shuffle", "region", "homogeneity"))
  if (null_model == "shuffle") {
  	# Chase, J.M., Kraft, N.J.B., Smith, K.G., Vellend, M., Inouye, B.D., 2011. Using null models to disentangle variation in community dissimilarity from variation in α-diversity. Ecosphere 2, art24. https://doi.org/10.1890/es10-00117.1
    gamma <- ncol(comm)
    comm_b <- comm
    comm_b[comm_b > 0] <- 1
    occur <- colSums(comm_b)
    richness_levels <- rowSums(comm_b)

    null_comm <- t(sapply(1:nrow(comm), function(i) {
      cc <- rep(0, gamma)
      presence_id <- comm[i, ] > 0
      if (length(occur[which(presence_id)]) > 1) {
        shuffle_id <- sample(which(presence_id), richness_levels[i], replace = F, prob = occur[which(presence_id)])
        shuffle_num <- sample(comm[i, presence_id], richness_levels[i])
      } else {
        shuffle_id <- which(presence_id)
        shuffle_num <- comm[i, presence_id]
      }
      cc[shuffle_id] <- shuffle_num
      cc
    }))
  }
  if (null_model == "region") {
  	# Stegen, J.C., Lin, X., Fredrickson, J.K., Chen, X., Kennedy, D.W., Murray, C.J., Rockhold, M.L., Konopka, A., 2013. Quantifying community assembly processes and identifying features that impose them. The ISME Journal 7, 2069–2079. https://doi.org/10.1038/ismej.2013.93
    gamma <- ncol(comm)
    sp_abun <- colSums(comm)
    sample_richness <- rowSums(comm)

    comm_b <- comm
    comm_b[comm_b > 0] <- 1
    richness_levels <- rowSums(comm_b)

    null_comm <- t(sapply(1:nrow(comm), function(i) {
      if (length(which(comm[i, ] > 0)) > 1) {
        c1 <- c2 <- rep(0, gamma)
        presence_id <- comm[i, ] > 0
        c1[which(presence_id)] <- 1
        samp_sub <- sample(which(presence_id), sample_richness[i] - richness_levels[i], replace = T, prob = sp_abun[presence_id])
        c2[as.numeric(names(table(samp_sub)))] <- table(samp_sub)
        cc <- c1 + c2
      } else {
        cc <- comm[i, ]
      }
      cc
    }))
  }
  if (null_model == "homogeneity") {
  	# Zhang, X., Liu, S., Wang, J., Huang, Y., Freedman, Z., Fu, S., Liu, K., Wang, H., Li, X., Yao, M., Liu, X., Schuler, J., 2020. Local community assembly mechanisms shape soil bacterial β diversity patterns along a latitudinal gradient. Nature Communications 11, 5428. https://doi.org/10.1038/s41467-020-19228-4
    if (is.null(rownames(comm))) rownames(comm) <- seq(nrow(comm))
    if (is.null(colnames(comm))) colnames(comm) <- paste0("Species", 1:ncol(comm))

    sp_abun <- colSums(comm)
    sample_richness <- rowSums(comm)
    regional_abun <- sum(sp_abun)

    if (sp_freq == "fix") {
      ind_pool <- rep(colnames(comm), sp_abun)
    }
    if (sp_freq == "shuffle") {
      ind_pool <- c(colnames(comm), sample(colnames(comm), regional_abun - ncol(comm), replace = T))
      ind_pool <- sample(ind_pool, length(ind_pool))
    }
    if (samp_rich == "fix") {
      samp_pool <- sample(rep(rownames(comm), sample_richness), regional_abun)
    }
    if (samp_rich == "shuffle") {
      samp_pool <- sample(rownames(comm), regional_abun, replace = TRUE)
    }
    null_comm <- tapply(rep(1, regional_abun), list(samp_pool, ind_pool), sum)
    null_comm[is.na(null_comm)] <- 0
    null_comm <- null_comm[, match(colnames(comm), colnames(null_comm))]
    null_comm <- null_comm[match(rownames(comm), rownames(null_comm)), ]
  }
  null_dist <- vegan::vegdist(null_comm, method = dist_method)
  if (!output_nullcomm) {
    return(null_dist)
  } else {
    return(list(null_comm = null_comm, null_dist = null_dist))
  }
}

summStoc <- function(ob_dis, per_dis_l, stats = c("RC", "SES", "ST", "MST", "NST")) {
  ob_sim <- 1 - ob_dis
  # per_ma_l <- lapply(1:length(per_madis), function(k) per_madis[[k]][[1]])
  # per_dis_l <- lapply(1:length(per_madis), function(k) per_madis[[k]][[2]])
  per_dis <- Reduce("+", per_dis_l) / length(per_dis_l)
  per_sim <- 1 - per_dis

  RC_Radio <- SES_Radio <- ST_Radio <- MST_Radio <- NST_Radio <- NA
  xx <- NULL

  if ("RC" %in% stats) {
    r1 <- sapply(1:length(per_dis_l), function(k) ifelse(ob_dis > per_dis_l[[k]], 1, 0)) %>% rowSums()
    r2 <- sapply(1:length(per_dis_l), function(k) ifelse(ob_dis == per_dis_l[[k]], 1, 0)) %>% rowSums()

    RC_Radio <- ((r1 + r2 * 0.5) / length(per_dis_l) - 0.5) * 2
  }

  if ("SES" %in% stats) {
    # bdisp_ob <- vegan::betadisper(ob_dis, group = rep(1, nrow(as.matrix(ob_sim))), type = "centroid")$distances
    # bdisp_per <- lapply(per_dis_l, vegan::betadisper, group = rep(1, times = nrow(as.matrix(ob_sim))), type = "centroid") %>% sapply("[[", "distances")
    SES_Radio <- (c(ob_dis) - rowMeans(sapply(per_dis_l, c))) / apply(sapply(per_dis_l, c), 1, stats::sd)
    # SES_Radio <- (bdisp_ob - rowMeans(bdisp_per)) / apply(bdisp_per, 1, stats::sd)
    # SES_Radio <- c(SES_Radio, rep(NA, length(ob_sim) - length(SES_Radio)))
    # SES_Radio <- as.vector((ob_dis - per_dis) / apply(Reduce("cbind", per_dis_l), 1, sd))
  }

  if ("ST" %in% stats) {
    ST_Radio1 <- as.vector(per_sim / ob_sim)
    ST_Radio <- ST_Radio1
    ST_Radio[which(ST_Radio1 > 1)] <- ((1 - per_sim) / (1 - ob_sim))[which(ST_Radio1 > 1)]
    ST_Radio[which(is.nan(ST_Radio1))] <- ((1 - per_sim) / (1 - ob_sim))[which(is.nan(ST_Radio1))]
  }

  if ("MST" %in% stats) {
    MST_Radio <- as.vector((per_sim / ob_sim) * ((1 - ob_sim) / (1 - per_sim)))
    MST_Radio[which(MST_Radio > 1)] <- ((ob_sim / per_sim) * ((1 - per_sim) / (1 - ob_sim)))[which(MST_Radio > 1)]
    MST_Radio[which(is.nan(MST_Radio))] <- 1
  }

  if ("NST" %in% stats) {
    NST_Radio <- NST::cNST(ob_dis, per_dis_l, group = data.frame(rep(1, nrow(as.matrix(ob_sim)))))$index.pair.grp$NST.ij.dis
  }
  xx <- data.frame(RC_Radio, SES_Radio, ST_Radio, MST_Radio, NST_Radio)[, match(stats, c("RC", "SES", "ST", "MST", "NST")), drop = F]
  names(xx) <- stats
  as.data.table(xx)
}


nullStoc <- function(comm, null_model, sp_freq = "fix", samp_rich = "fix", dist_method = "bray", reps = 1000, nworker = 1, stats = c("RC", "SES", "ST", "MST", "NST")) {
  comm <- as.matrix(comm)
  comm <- comm[rowSums(comm) > 0, ]
  comm <- comm[, colSums(comm) > 0]

  # make_null:
  if (nworker > 1) {
    # parallel computing
    cat(paste("Now parallel computation randomization. Begin at ", date(), ". Please wait ...", sep = ""), "\n")
    cl <- makeCluster(nworker, type = "PSOCK", setup_strategy = "sequential") # PSOCK (windows); FORK (linux);
    registerDoParallel(cl)
    per_dis_l <- foreach(1:reps, .packages = "MicroEcoTk") %dopar% nullModel(comm = comm, null_model = null_model, sp_freq = sp_freq, samp_rich = samp_rich, output_nullcomm = F, dist_method = dist_method)
    stopCluster(cl)
  } else {
    per_dis_l <- lapply(1:reps, function(k) nullModel(comm = comm, null_model = null_model, sp_freq = sp_freq, samp_rich = samp_rich, output_nullcomm = F, dist_method = dist_method))
  }

  ob_dis <- vegan::vegdist(comm, dist_method)
  res <- summStoc(ob_dis, per_dis_l, stats)
  cat(paste0("Ending at ", date(), ". "), "\n")
  return(res)
}
