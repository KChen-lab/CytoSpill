#' @export
GetSpillMat <- function(data = NULL, cols, n, file = NULL, threshold, flexrep, neighbor) {
  if (is.null(data)){
    data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))
  }
  cutoffs <- .DeriveCutoffs(data, cols, n, flexrep)
  model <- .EstimateSpill(data, cutoffs, cols, upperbound = threshold, neighbor = neighbor)
  estimates <- model[[1]]
  xcols <- model[[2]]
  spillmat <- diag(length(xcols))
  for (i in 1:length(xcols)) {
    if (!is.na(xcols[[i]][1])) {
      for (j in 1:length(xcols[[i]])) {
        if (!is.na(estimates[[i]][j])){
          spillmat[xcols[[i]][j],i] <- min(estimates[[i]][j],threshold)
        }
      }
    }
  }
  return(list(spillmat,cutoffs))
}
