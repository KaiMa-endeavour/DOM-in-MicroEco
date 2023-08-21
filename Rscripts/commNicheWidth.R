# @Project: DOM-in-MicroEco
# @File: commNicheWidth.R
# @Author: Kai Ma
# Date: 26/06/2023


commNicheWidth <- function(comm, nworker = 4) {
  require(parallel)
  
  comm <- as.data.frame(comm)
  taxa_nicheW <- spaa::niche.width(comm, "levins")
  comm[comm == 0] <- NA
  cl <- makeCluster(nworker, type = "PSOCK")
  xx <- parallel::parSapply(cl, 1:nrow(comm), function(r) mean(taxa_nicheW[!is.na(comm[r, ])]))
  stopCluster(cl)
  return(xx)
}
