# @Project: DOM-in-MicroEco
# @File: commDispersal.R
# @Author: Kai Ma
# Date: 23/05/2023


commDispersal <- function(comm, d, nworker = 4) {
  # A community of metacommunities: exploring patterns in species distributions across large geographical areas
  require(parallel)
  require(foreach)
  require(doParallel)
  
  comm[comm > 0] <- 1
  d <- d[match(rownames(comm), rownames(d)), match(rownames(comm), colnames(d))]
  xx <- matrix(NA, nrow(comm), nrow(comm))
  colnames(xx) <- rownames(comm)
  rownames(xx) <- rownames(comm)

  destination <- function(j) {
    two_sample_spec <- comm[i, ] + comm[j, ]
    two_sample_spec[two_sample_spec > 0] <- 1
    m <- sum(two_sample_spec)
    first <- sum(sapply(1:m, function(k) comm[j, k] * 1 / exp(-log10(d[i, j]))))
    second <- first / m / (nrow(comm) - 1)
    second
  }
  for (i in 1:(nrow(comm) - 1)) {
    # i=1
    cat(paste0("Start calculating dispersal (", i, "/", nrow(comm) - 1, ")    ", date(), " ..."), "\n")
    # second <- rep(0, nrow(comm))
    cl <- makeCluster(nworker, type = "PSOCK")
    registerDoParallel(cl)
    second <- foreach(j = (i + 1):nrow(comm), .combine = "c") %dopar% destination(j)
    stopCluster(cl)
    xx[(i + 1):nrow(comm), i] <- second
    cat(c(i, "Done.   ", date()), "\n")
  }
  stats::as.dist(xx)
}
