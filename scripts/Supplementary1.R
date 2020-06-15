#generate supplementary 1 results (based on reviewer 2 comment 3)
#simulation
r2c3_simulation_data <- function(){
  Nd142Di <- c(rnorm(10000, 200, 50), rep(0, 10000))
  Nd142Di[which(Nd142Di < 0)] <- 0

  Nd144Di <- c(rnorm(10000, 300, 75), rep(0, 10000))
  Nd144Di[which(Nd144Di < 0)] <- 0

  Nd145Di <- c(rnorm(10000, 400, 100), rep(0, 10000))
  Nd145Di[which(Nd145Di < 0)] <- 0

  Nd143Di <- c(rep(0, 10000), rnorm(10000, 300, 75))
  Nd143Di[which(Nd143Di < 0)] <- 0

  data_1 <- cbind(Nd142Di, Nd143Di, Nd144Di, Nd145Di)
  return(data_1)
}

set.seed(1)
data_channelwithvalue <- r2c3_simulation_data()
#different value
data_channelwithvalue_level2 <- data_channelwithvalue
data_channelwithvalue_level2[,2] <- data_channelwithvalue[,2]*0.5

data_channelwithvalue_level3 <- data_channelwithvalue
data_channelwithvalue_level3[,2] <- data_channelwithvalue[,2]*2

data_channelwithvalue_level4 <- data_channelwithvalue
data_channelwithvalue_level4[,2] <- data_channelwithvalue[,2]*0.2

#empty channel no value present
data_emptychannel <- data_channelwithvalue
data_emptychannel[,2] <- rep(0, 20000)

#get_spill_cols from CATALYST package
get_spill_cols <- function(ms, mets, l=CATALYST::isotope_list) {
  ms <- as.numeric(ms)
  spill_cols <- vector("list", length(ms))
  for (i in seq_along(ms)) {
    p1 <- m1 <- ox <- iso <- NULL
    if ((ms[i] + 1)  %in% ms) p1 <- which(ms == (ms[i] + 1))
    if ((ms[i] - 1)  %in% ms) m1 <- which(ms == (ms[i] - 1))
    if ((ms[i] + 16) %in% ms) ox <- which(ms == (ms[i] + 16))
    iso <- l[[mets[i]]]
    iso <- which(ms %in% iso[iso != ms[i]])
    spill_cols[[i]] <- unique(c(m1, p1, iso, ox))
  }
  spill_cols
}

spill_cols_data <- function(data){
  ##get ms and mets
  chs <- colnames(data)
  ####metal mass number like 167，***
  ms <- as.numeric(regmatches(chs, gregexpr("[0-9]+", chs)))
  ####metal name
  mets <- gsub("[[:digit:]]+Di", "", chs)
  ##get spillover cols
  spill_cols <- get_spill_cols(ms, mets)
  return(spill_cols)
}

#function generate a simulated spillover matrix
sm_simulation <- function(data, n, ...){
  ##get ms and mets
  chs <- colnames(data)
  ####metal mass number like 167，***
  ms <- as.numeric(regmatches(chs, gregexpr("[0-9]+", chs)))
  ####metal name
  mets <- gsub("[[:digit:]]+Di", "", chs)
  ##get spillover cols
  spill_cols <- get_spill_cols(ms, mets)

  sm <- vector("list",n)
  for (i in 1:n){
    m = diag(length(spill_cols))
    for (j in seq_along(spill_cols)){
      for (k in spill_cols[[j]]){
        m[j,k] <- runif(1, min=0, max=0.1)
      }
    }
    sm[[i]] <- m
  }
  return(sm)
}

#function generate simulated data
data_simulation <- function(signal, sm, ...){
  data <- vector("list",length(sm))
  for (i in seq_along(sm)){
    data[[i]] = signal%*%sm[[i]]
    colnames(data[[i]]) = colnames(signal)
  }
  return(data)
}

sm_simulated <- sm_simulation(data = data_channelwithvalue, n = 10)
##simulated data with different settings
data_simulated_1 <- data_simulation(signal = data_channelwithvalue, sm = sm_simulated)
data_simulated_2 <- data_simulation(signal = data_emptychannel, sm = sm_simulated)
data_simulated_level2 <- data_simulation(signal = data_channelwithvalue_level2, sm = sm_simulated)
data_simulated_level3 <- data_simulation(signal = data_channelwithvalue_level3, sm = sm_simulated)
data_simulated_level4 <- data_simulation(signal = data_channelwithvalue_level4, sm = sm_simulated)

estimate_sm <- function(data, spill_cols){
  my_sm = vector("list",length(data))
  my_sm_nmf = vector("list",length(data))
  for (i in seq_along(data)){
    cutoffs_sim <- .DeriveCutoffs(data=data[[i]], cols=c(1:ncol(data[[i]])), n = 20000, flexrep = 10)
    threshold <- 0.1
    model <- .EstimateSpill(data=data[[i]], cutoffs=cutoffs_sim, cols=c(1:ncol(data[[i]])), upperbound = threshold, neighbor = 1)
    estimates <- model[[1]]
    xcols <- model[[2]]
    spillmat <- diag(length(xcols))
    for (m in 1:length(xcols)) {
      if (!is.na(xcols[[m]][1])) {
        for (j in 1:length(xcols[[m]])) {
          if (!is.na(estimates[[m]][j])){
            spillmat[xcols[[m]][j],m] <- min(estimates[[m]][j],threshold)
          }
        }
      }
    }
    sm_sim <- spillmat
    my_sm[[i]] <- sm_sim
  }
  return(my_sm)
}

spill_cols <- spill_cols_data(data_channelwithvalue)
estimate_sm_withvalue <- estimate_sm(data_simulated_1, spill_cols)
estimate_sm_withoutvalue <- estimate_sm(data_simulated_2, spill_cols)
estimate_sm_withvalue_level2 <- estimate_sm(data_simulated_level2, spill_cols)
estimate_sm_withvalue_level3 <- estimate_sm(data_simulated_level3, spill_cols)
estimate_sm_withvalue_level4 <- estimate_sm(data_simulated_level4, spill_cols)

channel2_comparison <- function(x, y, spill_cols){
  comparison = NULL
  for (i in seq_along(x)){
    single <- NULL
    for (j in c(1,3,4)){
      for (k in c(2)){
        single <- rbind(single,c(x[[i]][j,k],y[[i]][j,k]))
      }
    }
    comparison <- rbind(comparison,single)
  }
  colnames(comparison) <- c('simulation','my_method')
  return(comparison)
}

withvalue_comparison <- channel2_comparison(sm_simulated, estimate_sm_withvalue, spill_cols)
withoutvalue_comparison <- channel2_comparison(sm_simulated, estimate_sm_withoutvalue, spill_cols)
withvaluevswithoutvalue <- channel2_comparison(estimate_sm_withvalue,estimate_sm_withoutvalue, spill_cols)

withvalue_level2_comparison <- channel2_comparison(sm_simulated, estimate_sm_withvalue_level2, spill_cols)
withvalue_level3_comparison <- channel2_comparison(sm_simulated, estimate_sm_withvalue_level3, spill_cols)
withvalue_level4_comparison <- channel2_comparison(sm_simulated, estimate_sm_withvalue_level4, spill_cols)


ggplot(as.data.frame(withvalue_comparison), aes(x=my_method, y=simulation)) +
  geom_point(size=0.1) + theme_classic() + xlab('Channel with signal estimated spillover effects') +
  ylab('Simulated spillover effects') + coord_fixed() + geom_smooth()

ggplot(as.data.frame(withoutvalue_comparison), aes(x=my_method, y=simulation)) +
  geom_point(size=0.1) + theme_classic() + xlab('Channel without signal estimated spillover effects') +
  ylab('Simulated spillover effects') + coord_fixed() + geom_smooth()

# ggplot(as.data.frame(withvaluevswithoutvalue), aes(x=my_method, y=simulation)) +
#   geom_point(size=0.1) + theme_classic() + xlab('Channel with signal estimated spillover effects') +
#   ylab('Channel without signal estimated spillover effects') + coord_fixed()

ggplot(as.data.frame(withvalue_level2_comparison), aes(x=my_method, y=simulation)) +
  geom_point(size=0.1) + theme_classic() + xlab('Channel with signal level2 estimated spillover effects') +
  ylab('Simulated spillover effects') + coord_fixed() + geom_smooth()

ggplot(as.data.frame(withvalue_level3_comparison), aes(x=my_method, y=simulation)) +
  geom_point(size=0.1) + theme_classic() + xlab('Channel with signal level3 estimated spillover effects') +
  ylab('Simulated spillover effects') + coord_fixed() + geom_smooth()

ggplot(as.data.frame(withvalue_level4_comparison), aes(x=my_method, y=simulation)) +
  geom_point(size=0.1) + theme_classic() + xlab('Channel with signal level4 estimated spillover effects') +
  ylab('Simulated spillover effects') + coord_fixed() + geom_smooth()


cor(withvalue_comparison[,1],withvalue_comparison[,2])^2
#0.9940178

cor(withoutvalue_comparison[,1],withoutvalue_comparison[,2])^2
#0.9940971

cor(withvalue_level2_comparison[,1],withvalue_comparison[,2])^2
#0.9940178

cor(withvalue_level3_comparison[,1],withvalue_comparison[,2])^2
#0.9940178

cor(withvalue_level4_comparison[,1],withvalue_comparison[,2])^2
#0.9940178
