# @Project: DOM-in-MicroEco
# @File: VPA_dbmem.R
# @Author: Kai Ma
# Date: 16/10/2020


VPA_dbmem <- function(comm, env, geo, climate = NULL, for_sel = T, alpha = 0.05, Lon_Lat = T, permutations = 9999) {
  
  if (!require(SoDA)) {devtools::install_github('https://github.com/cran/SoDA')}; require(SoDA)
  if (!require(vegan)) {install.packages("vegan")}; require(vegan)
  if (!require(tidyfst)) {install.packages("tidyfst")}; require(tidyfst)
  if (!require(adespatial)) {install.packages("adespatial")}; require(adespatial)
  if (!require(geosphere)) {install.packages("geosphere")}; require(geosphere)
  
  env <- as.data.frame(env)
  geo <- as.data.frame(geo)
  comm <- as.data.frame(comm)
  
  dist_comm <- vegdist(comm, method = "bray") %>% as.matrix()
  if (for_sel) {
    (env_fwd <- forward.sel(dist_comm, env, alpha = alpha, nperm = permutations))
    Env_se <- env[ ,c(sort(env_fwd[ , 2])), drop = F]
  } else {
    Env_se <- env
  }

  PCsign_cap <- capscale(dist_comm~., data = Env_se, add = T)
  Env_P <- anova(PCsign_cap, permutations = permutations)[[4]][1]

  if (!is.null(climate)) {
    climate <- as.data.frame(climate)
    if (for_sel) {
      (cli_fwd <- forward.sel(dist_comm, climate, alpha = alpha, nperm = permutations))
      Cli_se <- climate[ ,c(sort(cli_fwd[ , 2])), drop = F]
    } else {
      Cli_se <- climate
    }
    
    PCsign_cap <- capscale(dist_comm~., data = Cli_se, add = T)
    Cli_P <- anova(PCsign_cap, permutations = permutations)[[4]][1]
  }
  
  if (Lon_Lat) {
    dbmem <- distm(geo, fun = distHaversine) %>% as.dist() %>% dbmem() 
  } else {
    dbmem <- dist(geo) %>% dbmem()
  }
  
  if (for_sel) {
    (pc_fwd <- forward.sel(dist_comm, dbmem, alpha = alpha, nperm = permutations))
    Geo_se <- dbmem[ ,c(sort(pc_fwd[ ,2])), drop = F]
  } else {
    Geo_se <- dbmem
  }

  PCsign_cap <- capscale(dist_comm~., data = Geo_se, add = T)
  Geo_P <- anova(PCsign_cap, permutations = permutations)[[4]][1]
  
  
  if (is.null(climate)) {
    Var_mod <- varpart(dist_comm, Env_se, Geo_se)
    Env_R2 <- Var_mod[1]$part[[2]][1, 3]
    Geo_R2 <- Var_mod[1]$part[[2]][2, 3]
    
    Spe_R2 <- Var_mod[1]$part[[3]][1, 3]
    Dis_R2 <- Var_mod[1]$part[[3]][3, 3]
    Com_R2 <- Var_mod[1]$part[[3]][2, 3]
    Res <- Var_mod[1]$`part`[[3]][4, 3]
    
    Spe_P <- anova(capscale(dist_comm~.+Condition(as.matrix(Geo_se)), data = Env_se), permutations = permutations)[[4]][1]
    Dis_P <- anova(capscale(dist_comm~.+Condition(as.matrix(Env_se)), data = Geo_se), permutations = permutations)[[4]][1]
    
    result <- mat.or.vec(6, 2)
    result[1, 1] <- Env_R2
    result[2, 1] <- Geo_R2
    result[3, 1] <- Spe_R2
    result[4, 1] <- Dis_R2
    result[5, 1] <- Com_R2
    result[6, 1] <- Res
    
    result[1, 2] <- Env_P
    result[2, 2] <- Geo_P
    result[3, 2] <- Spe_P
    result[4, 2] <- Dis_P
    result[5, 2] <- NA
    result[6, 2] <- NA
    colnames(result) <- c("Adj.R.squared","p-value")
    xx <- data.table(x = c("Enviromental","Geographic","Species sorting","Dispersal limitation","Combinated fraction", "Residuals "), result)
    xx$Adj.R.squared <- round(xx$Adj.R.squared*100, 3)
    xx$`p-value` <- round(xx$`p-value`, 5)
  } else {
    Var_mod <- varpart(dist_comm, Env_se, Geo_se, Cli_se)
    # Var_mod$part$indfract
    Spe_P <- anova(capscale(dist_comm~.+Condition(as.matrix(Geo_se), as.matrix(Cli_se)), data = Env_se), permutations = permutations)[[4]][1]
    Dis_P <- anova(capscale(dist_comm~.+Condition(as.matrix(Env_se), as.matrix(Cli_se)), data = Geo_se), permutations = permutations)[[4]][1]
    cli_P <- anova(capscale(dist_comm~.+Condition(as.matrix(Env_se), as.matrix(Geo_se)), data = Cli_se), permutations = permutations)[[4]][1]
    xx <- cbind(Var_mod$part$indfract, `p-value` = c(Spe_P, Dis_P, cli_P, rep('', 5)))
    xx$Adj.R.square <- round(xx$Adj.R.square*100, 3)
  }
  return(xx)
}
