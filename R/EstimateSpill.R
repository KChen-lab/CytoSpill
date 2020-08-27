.EstimateSpill <- function(data, cutoffs, cols, upperbound, neighbor) {
  results <- list()
  data <- data[,cols]
  spill_cols <- .SpillColsData(data, l= CytoSpill::isotope_list, nneighbor=neighbor)
  xcols <- .GetFmla(data, spill_cols = spill_cols)

  for (i in 1:ncol(data)) {
    if (!is.na(xcols[[i]][1])) {

      A = as.matrix(data[which(data[,i] < cutoffs[i]), c(i,xcols[[i]])])
      for (j in seq_along(xcols[[i]])){
        A = A[which(A[,j+1]>=cutoffs[xcols[[i]][j]]),]
        if (is.null(dim(A))) {A = matrix(A, nrow = 1)}
      }
      ##
      print(dim(A))
      ##
      if (nrow(A) == 0) { #where there is only 1 row remains
        results[[i]] <- NA
        print("no datapoint available for this channel")
      } else {
        if (nrow(A) == 1) {
          b = A[,1]
          A = matrix(A[,-1], nrow = 1)
        } else {
          b = A[,1]
          A = as.matrix(A[,-1])
        }
        x0 = runif(ncol(A), min = 0, max = upperbound)
        fn = function(x) {
          vec = A%*%x-b
          norm(vec, type = "2")
        }
        result = try(nloptr::slsqp(x0, fn, lower=rep(0, length(x0)), upper = rep(upperbound, length(x0))))
        if (isTRUE(class(result) == "try-error")) {
          result <- NA
          xcols[[i]] <- NA
        } else{
          result <- result$par
        }
        results[[i]] <- result
        print("sqp modelling done")
      }
    } else {
      results[[i]] <- NA
      print("no spillover cols for this channel")
    }
  }
  return(list(results, xcols))
}

.GetFmla <- function(data, spill_cols) {
  fmlacols <- list()
  cs <- c(1:ncol(data))
  for (i in 1:ncol(data)) {
    cols <- NULL
    for (j in seq_along(spill_cols)) {
      if (i %in% spill_cols[[j]]) {
        cols <- cbind(cols, j)
      }
    }
    if(is.null(cols)) {
      fmlacols[[i]] <- NA
    } else{
      fmlacols[[i]] <- cols
    }
  }
  return(fmlacols)
}

.SpillColsData <- function(data, l, nneighbor) {
  # get ms and mets
  chs <- colnames(data)
  # metal mass number like 167，***
  ms <- as.numeric(regmatches(chs, gregexpr("[0-9]+", chs)))
  # metal name
  mets <- gsub("[[:digit:]]+Di", "", chs)
  # get spillover cols
  spill_cols <- vector("list", length(ms))
  if (nneighbor == 1){
    for (i in seq_along(ms)) {
      p1  <- m1  <- ox <- iso <- NULL
      if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
      if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1))
      if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
      iso <- l[[mets[i]]]
      iso <- which(ms %in% iso[iso != ms[i]])
      spill_cols[[i]] <- unique(c(m1, p1, iso, ox))
    }
  }

  if (nneighbor == 2){
    for (i in seq_along(ms)) {
      p1 <- p2 <- m1 <- m2 <- ox <- iso <- NULL
      if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
      if ((ms[i] + 2)  %in% ms) p2 <- which(ms == (ms[i] + 2))
      if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1))
      if ((ms[i] - 2)  %in% ms) m2 <- which(ms == (ms[i] - 2))
      if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
      iso <- l[[mets[i]]]
      iso <- which(ms %in% iso[iso != ms[i]])
      spill_cols[[i]] <- unique(c(m1, m2, p1, p2, iso, ox))
    }
  }

  return(spill_cols)
}
