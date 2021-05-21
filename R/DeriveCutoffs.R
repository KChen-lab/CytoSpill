.DeriveCutoffs <- function(data, cols, n, flexrep){
  data_sample <- data[sample(1:nrow(data), n, replace = FALSE), ]
  data_work <- data_sample[,cols]
  cutoffs <- .DeriveCutoffsHelper(x = data_work, quantile = 0.1, flexrep = flexrep)
  return(cutoffs)
}

.DeriveCutoffsHelper <- function(x, quantile, flexrep){
  cutoffs <- rep(NA, dim(x)[2])

  for (i in 1:dim(x)[2]) {
    s <- x[, i][which(x[, i] > 0)]

    model_fit <- flexmix::initFlexmix(s~1, k = c(1,2,3,4,5), nrep = flexrep)
    model_chosen <- flexmix::getModel(model_fit, which = "ICL")

    if (length(model_chosen@prior)>1){
      param_list <- list()
      param <- flexmix::parameters(model_chosen)
      func_list <- list()
      #get parameters for each component, prior probability, mean, sd
      k <- length(model_chosen@prior)
      if (k == 2){
        d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
        d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
        func_list = list(d1,d2)
      }
      if (k == 3){
        d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
        d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
        d3 <- function(x) model_chosen@prior[3]*dnorm(x, param[1,3], param[2,3])
        func_list = list(d1,d2,d3)
      }
      if (k == 4){
        d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
        d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
        d3 <- function(x) model_chosen@prior[3]*dnorm(x, param[1,3], param[2,3])
        d4 <- function(x) model_chosen@prior[4]*dnorm(x, param[1,4], param[2,4])
        func_list = list(d1,d2,d3,d4)
      }
      if (k == 5){
        d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
        d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
        d3 <- function(x) model_chosen@prior[3]*dnorm(x, param[1,3], param[2,3])
        d4 <- function(x) model_chosen@prior[4]*dnorm(x, param[1,4], param[2,4])
        d5 <- function(x) model_chosen@prior[5]*dnorm(x, param[1,5], param[2,5])
        func_list = list(d1,d2,d3,d4,d5)
      }

      i_min <- which.min(param[1,])
      i_max <- which.max(param[1,])
      min_modal <- eval(parse(text=paste("d",i_min,sep = "")))
      max_modal <- eval(parse(text=paste("d",i_max,sep = "")))
      func_list <- list(min_modal, max_modal)
      sum_func <- function (x){sum(sapply(func_list, function(f) f(x)))}
      type1 <- 0.05
      p <- function(x) {
        min_modal(x)/(sum_func(x))-1+type1
      }
      interval <- c(0, max(param[1,]))

      if (p(0)>0){
        root <- uniroot(p,interval)$root[1]
        cutoffs[i] <- root
      } else{
        cutoffs[i] <- quantile(s, probs = quantile)
        print("using quantile")
      }
    } else{
      cutoffs[i] <- quantile(s, probs = quantile)
      print("using quantile")
    }
    print(cutoffs[i])
  }
  return(cutoffs)
}
