#' @export
SpillComp <- function(file = NULL, data = NULL, cols, n, output = NULL, threshold = 0.1, flexrep = 10, neighbor = 1) {
  if (is.null(data)){
    data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))
  }
  spillmat_results <- GetSpillMat(data, cols, n, threshold = threshold, flexrep = flexrep, neighbor = neighbor)
  spillmat <- spillmat_results[[1]]
  cutoffs <- spillmat_results[[2]]
  data_compensated <- t(apply(data[,cols], 1, function(row) nnls::nnls(t(spillmat), row)$x))
  data_colnames <- colnames(data)
  data[,cols] <- data_compensated
  colnames(data) <- data_colnames
  if (!is.null(output)) {flowCore::write.FCS(flowCore::flowFrame(data), filename = output)}
  return(list(flowCore::flowFrame(data), spillmat, cutoffs))
}
