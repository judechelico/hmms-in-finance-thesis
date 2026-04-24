library(depmixS4)
library(ggplot2)

# GPT corrected trading algorithm - review and understand
# only let pvals = TRUE if you have one tolerance and two states
# set as true by default to avoid time-expensive mistakes
tradingAlg = function(returns,
                      tols = c(0.005, 0.01, 0.015),
                      days = 30, nStates = 3,
                      pvals = TRUE) {
  
  # a smaller tolerance makes it less likely for 
  # the algorithm to buy 
  # days is the number of days back the HMM goes 
  # per iteration
  
  len = length(returns)
  lentols = length(tols)
  nTrades = len - days
  
  ## storage
  price = matrix(NA, nrow = lentols, ncol = nTrades)
  decision = matrix(NA, nrow = lentols, ncol = nTrades)
  buyAndHold = numeric(nTrades)
  expVolNext = numeric(nTrades)
  if (nStates == 2){
    runMeans = numeric(nTrades)
    revMeans = numeric(nTrades)
    runSds = numeric(nTrades)
    revSds = numeric(nTrades)
  }
  
  if (pvals == TRUE){
    pvalList = numeric(nTrades)
  }
  
  for (i in (days + 1):len) {
    
    t = i - days
    
    # progress checker 
    prog = round(i / len, 3) * 100 
    if ((prog %% 5) == 0){ 
      cat(paste(prog, "% complete\n")) 
    }
    
    ## rolling window (past only)
    sub = returns[t:(i - 1)]
    df = data.frame(sub)
    
    ## fit HMM
    mod = depmix(sub ~ 1, data = df,
                 nstates = nStates, family = gaussian())
    
    for (j in 1:100) {
      ok = TRUE
      tryCatch(
        fit <- fit(mod, verbose = FALSE),
        error = function(e) ok <<- FALSE
      )
      if (ok) break
    }
    
    ## extract state SDs safely
    sds = matrix(sapply(1:nStates, function(s) {
      fit@response[[s]][[1]]@parameters$sd
    }), nrow = nStates)
    
    if (nStates == 2){
      ## extract state means safely
      means = matrix(sapply(1:nStates, function(s) {
        fit@response[[s]][[1]]@parameters$coefficients
      }), nrow = nStates)
      runMeans[t] = max(means)
      revMeans[t] = min(means)
      runSds[t]   = sds[which.max(means)]
      revSds[t]   = sds[which.min(means)]
    }
    
    ## posterior probs for last observation
    post  <- posterior(fit, type = "filtering")
    last_probs  <- matrix(as.numeric(post[nrow(post), -1]),
                          nrow = nStates)
    
    # transition matrix
    tmat = matrix(
      as.numeric(getpars(fit)[
        (nStates + 1):(nStates + nStates^2)
        ]), nrow = nStates, byrow = TRUE
    )
    
    ## expected volatility next period
    evNext = sum(t(last_probs) %*% tmat %*% sds)
    expVolNext[t] = evNext
    
    ## buy-and-hold recursion
    buyAndHold[t] = ifelse(
      t == 1,
      100 * (1 + returns[i]),
      buyAndHold[t - 1] * (1 + returns[i])
    )
    
    if (pvals == TRUE){
      if (lentols > 1 | nStates != 2){
        stop(
          "Only set pvals = TRUE for one tolerance and two states!"
          )
      }
      state = posterior(fit)$state
      in1 = sub[state == 1]
      in2 = sub[state == 2]
      tryCatch(
        {pval = t.test(in1, in2)$p.val},
        error = function(e){pval = NA},
        finally = {pvalList[t] = pval}
      )
    }
    
    ## grid over tolerances
    for (k in 1:lentols) {
      
      position = as.integer(evNext < tols[k])
      decision[k, t] = position
      
      price[k, t] = ifelse(
        t == 1,
        100 * (1 + position * returns[i]),
        price[k, t - 1] * (1 + position * returns[i])
      )
    }
  }
  
  if (pvals == TRUE){
    returnMe = list(
      algoPaths   = price,
      bnhPath     = buyAndHold,
      finalAlgo   = price[, nTrades],
      finalBnH    = buyAndHold[nTrades],
      meanEV      = mean(expVolNext),
      sdEV        = sd(expVolNext),
      meanPos     = rowMeans(decision),
      pvals       = pvalList
    )
  } else {
    returnMe = list(
      algoPaths   = price,
      bnhPath     = buyAndHold,
      finalAlgo   = price[, nTrades],
      finalBnH    = buyAndHold[nTrades],
      meanEV      = mean(expVolNext),
      sdEV        = sd(expVolNext),
      meanPos     = rowMeans(decision)
    )
  }
  if (nStates == 2){
    returnMe[['runMeans']] = runMeans
    returnMe[['revMeans']] = revMeans
    returnMe[['runSds']]   = runSds
    returnMe[['revSds']]   = revSds
  }
  
  return(returnMe)
}

setwd("~/Courses/Senior/Math Thesis/Data")

data = read.csv("SPY6MoDaily.csv")

# data comes in reverse index order, so we have to reverse it
data$Exchange.Date = rev(data$Exchange.Date)
data$Close = rev(data$Close)
data$Change = rev(data$Change)

ret = data$Change
traded = tradingAlg(ret, tol = c(0.01),
                    days = 30, nStates = 2)

traded2 = tradingAlg(ret, tol = seq(0.005, 0.05, by = 0.001),
                    days = 30, nStates = 2, pvals = FALSE)

# the BIG test (daily)
spy = read.csv("spy.csv")
spy$X.Chg = as.numeric(gsub("%", "", spy$X.Chg)) / 100
spy$X.Chg[is.na(spy$X.Chg)] = 0
spyRet = rev(spy$X.Chg)

tol = seq(0.001, 0.075, by = 0.001)
traded3S = tradingAlg(spyRet, 
                    tol = tol, pvals = FALSE,
                    days = 100, nStates = 3)

# for pvals (2 states), best tol is 0.028
traded2 = tradingAlg(spyRet, tol = 0.028, days = 100,
                     pvals = TRUE, nStates = 2)

tol = seq(0.001, 0.1, by = 0.001)
traded2S = tradingAlg(spyRet[1501:length(spyRet)], 
                      tol = tol, pvals = FALSE,
                      days = 100, nStates = 2)

traded2STrain = tradingAlg(spyRet[1:1500],
                           tol = tol, pvals = FALSE,
                           days = 100, nStates = 2)
best = tol[which.max(traded2STrain$finalAlgo)]

traded2STest = tradingAlg(spyRet[1501:length(spyRet)],
                          tol = best, pvals = FALSE,
                          days = 100, nStates = 2)

# export data
names = paste("tol", tol)
export = as.data.frame(
  matrix(rep(0, ncol(traded2S$algoPaths) * (length(tol) + 8)), 
         nrow = ncol(traded2S$algoPaths))
)
names(export) = c("index", "runmean", "revmean",
                  "runsd", "revsd", "meanEV", "sdEV", 
                  "bnh", names)
export$index = c(1:ncol(traded2S$algoPaths))
export$bnh = traded2S$bnhPath
export$runmean = traded2S$runMeans
export$revmean = traded2S$revMeans
export$runsd = traded2S$runSds
export$revsd = traded2S$revSds
export$meanEV = traded2S$meanEV
export$sdEV = traded2S$meanEV
for (i in 1:length(tol)){
  export[i + 8] = traded2S$algoPaths[i, ]
}

write.csv(export, "bigSPY2sWithParams1501MoreTols.csv", 
          row.names = FALSE)

# the less big test (weekly)
weekSPY = read.csv("SPY5YearWeekly.csv")
weekSPYRet = rev(weekSPY$X.Chg)

traded = tradingAlg(weekSPYRet,
                    # expect higher variances week by
                    # week instead of day by day
                    tol = seq(0.005, 0.05, by = 0.001),
                    # here, 'days' is actually weeks
                    days = 52, nStates = 3)

p = 100
for (i in 53:100){
  p = p * (1 + weekSPYRet[i])
}

# minutely testing
minuteSPY = read.csv("SPY1mo5minly.csv")
minuteSPY$X.Chg = as.numeric(gsub("%", "", minuteSPY$X.Chg)) / 100
minuteSPY$X.Chg[is.na(minuteSPY$X.Chg)] = 0
minSPYRet = rev(minuteSPY$X.Chg)

tradedMin = tradingAlg(minSPYRet, 
                       tol = seq(0.0001, 0.0015, by = 0.00005),
                       # 'days' is actually 5min blocks
                       days = 60, nStates = 2, pvals = FALSE)
