.DeriveCutoffs <- function(data, cols, n, flexrep){
  data_sample <- data[sample(1:nrow(data), n, replace = FALSE), ]
  data_work <- data_sample[,cols]
  # cutoffs <- .DeriveCutoffsHelper(x = data_work, emt = 0.001, quantile = 0.1)
  cutoffs <- .DeriveCutoffsHelper(x = data_work, quantile = 0.1, flexrep = flexrep)
  return(cutoffs)
}

.DeriveCutoffsHelper <- function(x, quantile, flexrep){
  cutoffs <- rep(NA, dim(x)[2])

  for (i in 1:dim(x)[2]) {
    s <- x[, i][which(x[, i] > 0)]

    # s <- asinh(s/5)

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


      # i_min <- which.min(param[1,])
      # sum_func <- function (x){sum(sapply(func_list, function(f) f(x)))}
      # type1 <- 0.05
      # min_modal <- eval(parse(text=paste("d",i_min,sep = "")))
      # p <- function(x) {
      #   min_modal(x)/(sum_func(x))-1+type1
      # }
      # interval <- c(0, min(param[1,][param[1,]!=min(param[1,])]))

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
      # if (p(0) < (1-type1)){
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
    # cutoffs[i] <- sinh(cutoffs[i])*5
    print(cutoffs[i])
  }
  return(cutoffs)
}
#
# .DeriveCutoffsHelper <- function(x, quantile){
#   cutoffs <- rep(NA, dim(x)[2])
#
#   for (i in 1:dim(x)[2]) {
#     s <- x[, i][which(x[, i] > 0)]
#     s_asinh <- asinh(s / 5)
#     q = quantile(s_asinh, 0.999)
#     s_asinh[s_asinh > q] = q
#     model_fit <- flexmix::initFlexmix(s_asinh~1, k = c(1,2,3,4,5))
#     model_chosen <- flexmix::getModel(model_fit, which = "ICL")
#
#     if (length(model_chosen@prior)>1){
#       param_list <- list()
#       param <- flexmix::parameters(model_chosen)
#       func_list <- list()
#       #get parameters for each component, prior probability, mean, sd
#       k <- length(model_chosen@prior)
#       if (k == 2){
#         d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
#         d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
#         func_list = list(d1,d2)
#       }
#       if (k == 3){
#         d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
#         d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
#         d3 <- function(x) model_chosen@prior[3]*dnorm(x, param[1,3], param[2,3])
#         func_list = list(d1,d2,d3)
#       }
#       if (k == 4){
#         d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
#         d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
#         d3 <- function(x) model_chosen@prior[3]*dnorm(x, param[1,3], param[2,3])
#         d4 <- function(x) model_chosen@prior[4]*dnorm(x, param[1,4], param[2,4])
#         func_list = list(d1,d2,d3,d4)
#       }
#       if (k == 5){
#         d1 <- function(x) model_chosen@prior[1]*dnorm(x, param[1,1], param[2,1])
#         d2 <- function(x) model_chosen@prior[2]*dnorm(x, param[1,2], param[2,2])
#         d3 <- function(x) model_chosen@prior[3]*dnorm(x, param[1,3], param[2,3])
#         d4 <- function(x) model_chosen@prior[4]*dnorm(x, param[1,4], param[2,4])
#         d5 <- function(x) model_chosen@prior[5]*dnorm(x, param[1,5], param[2,5])
#         func_list = list(d1,d2,d3,d4,d5)
#       }
#       # for (i in seq_along(model_chosen@prior)){
#       #   #prior, mean, sd
#       #   param_list[[i]] <- c(model_chosen@prior[i], param[1,i], param[2,i])
#       #   var <- paste("d", i, sep = "")
#       #   assign(var, (function(x) model_chosen@prior[i]*dnorm(x, param[1,i], param[2,i])))
#       #   # assign(var, (function(x) substitute(model_chosen@prior[i]*dnorm(x, param[1,i], param[2,i]), list(i=i))))
#       #   func_list[[i]] <- eval(parse(text=var))
#       # }
#       i_min <- which.min(param[1,])
#       sum_func <- function (x){sum(sapply(func_list, function(f) f(x)))}
#       type1 <- 0.05
#       min_modal <- eval(parse(text=paste("d",i_min,sep = "")))
#       p <- function(x) {
#         min_modal(x)/(sum_func(x))-1+type1
#       }
#       interval <- c(0, min(param[1,][param[1,]!=min(param[1,])]))
#       if (p(0)>0){
#         root <- uniroot(p,interval)$root
#         cutoffs[i] <- (sinh(root)) * 5
#       } else{
#         cutoffs[i] <- quantile(s, probs = quantile)
#       }
#     } else{
#       cutoffs[i] <- quantile(s, probs = quantile)
#     }
#   }
#   return(cutoffs)
# }

# .DeriveCutoffsHelper <- function(x, emt, quantile) {
#   # arguments: x: data matrix, emt: EM criterion
#   # em, cutoff are the functions from package cutoff which are edited to suit the current R version
#   cutoffs <- rep(NA, dim(x)[2])
#   for (i in 1:dim(x)[2]) {
#     a <- x[, i][which(x[, i] > 0)]
#     a_asinh <- asinh(a / 5)
#     if (diptest::dip.test(a_asinh)[2] > 0.05) {
#       cutoff_channeli <- quantile(a_asinh, probs = quantile)
#     } else {
#       fit_channeli <- try(em(a_asinh, "normal", "log-normal", t = 0.001), silent=FALSE)
#       if (isTRUE(class(fit_channeli) == "try-error")) {
#         cutoff_channeli <- quantile(a_asinh, probs = quantile)
#       } else {
#         cutoff_channeli <- try(cutoff(fit_channeli, t = 0.001, nb = 10, distr = 1, type1 = 0.05, level = 0.95))
#         if (isTRUE(class(cutoff_channeli) == "try-error")) {
#           cutoff_channeli <- fit_channeli$param[1] + fit_channeli$param[2]
#         } else{
#           cutoff_channeli <- cutoff_channeli[1]
#         }
#       }
#     }
#     cutoff_channeli <- (sinh(cutoff_channeli)) * 5
#     if (cutoff_channeli <= min(a)){
#       cutoff_channeli <- quantile(a[which(a > 0)], probs = quantile)
#     }
#     cutoffs[i] <- cutoff_channeli
#   }
#   return(cutoffs)
# }
