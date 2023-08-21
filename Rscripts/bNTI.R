# @Project: DOM-in-MicroEco
# @File: bNTI.R
# @Author: Kai Ma
# Date: 18/7/2021


bNTI <- function(comm, tree, nworker = 4, reps = 1000, outgroup = NULL, output_path = NULL, abundance.weighted = T) {
  require(tidyfst)
  require(picante)
  require(ape)
  
  match_phylo <- match.phylo.comm(tree, comm)
  bMNTD <- as.vector(comdistnt(match_phylo$comm, stats::cophenetic(match_phylo$phy), abundance.weighted = abundance.weighted))
  if (!is.null(outgroup)) match_phylo$phy <- root(match_phylo$phy, outgroup = outgroup, resolve.root = T)
  
  if (nworker > 1) {
    require(parallel)
    require(foreach)
    require(doParallel)
    
    cl <- makeCluster(nworker, type = "PSOCK")
    registerDoParallel(cl)
    bMNTD_rand <- foreach(1:reps, .packages = "picante", .combine = cbind) %dopar% picante::comdistnt(match_phylo$comm, picante::taxaShuffle(stats::cophenetic(match_phylo$phy)), abundance.weighted = abundance.weighted, exclude.conspecifics = F)
    stopCluster(cl)
  } else {
    bMNTD_rand <- sapply(1:reps, function(k) {
      picante::comdistnt(match_phylo$comm, picante::taxaShuffle(stats::cophenetic(match_phylo$phy)), abundance.weighted = abundance.weighted, exclude.conspecifics = F)
    })
  }
  res <- (bMNTD - apply(bMNTD_rand, 1, mean)) / apply(bMNTD_rand, 1, stats::sd)
  if (!is.null(output_path)) {
    fwrite(data.frame(bMNTD = bMNTD, bNTI = res), output_path)
  } else {
    data.frame(bMNTD = bMNTD, bNTI = res)
  }
}
