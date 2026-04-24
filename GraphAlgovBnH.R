source("Trading.R")

# run this after running 'trading' file to get the function ----
# (and spyRet)
tol = seq(0.001, 0.1, by = 0.001)
traded2STrain = tradingAlg(spyRet[1:1500],
                           tol = tol, pvals = FALSE,
                           days = 100, nStates = 2)
tol = tol[which.max(traded2STrain$finalAlgo)]

traded2STest = tradingAlg(spyRet[1501:length(spyRet)],
                          tol = tol, pvals = FALSE,
                          days = 100, nStates = 2)

# export data
names = paste("tol", tol)
export = as.data.frame(
  matrix(rep(0, ncol(traded2STest$algoPaths) * (length(tol) + 8)), 
         nrow = ncol(traded2STest$algoPaths))
)
names(export) = c("index", "runmean", "revmean",
                  "runsd", "revsd", "meanEV", "sdEV", 
                  "bnh", names)
export$index = c(1:ncol(traded2STest$algoPaths))
export$bnh = traded2STest$bnhPath
export$runmean = traded2STest$runMeans
export$revmean = traded2STest$revMeans
export$runsd = traded2STest$runSds
export$revsd = traded2STest$revSds
export$meanEV = traded2STest$meanEV
export$sdEV = traded2STest$meanEV
for (i in 1:length(tol)){
  export[i + 8] = traded2STest$algoPaths[i, ]
}

write.csv(export, "bigSPY2sTest.csv", 
          row.names = FALSE)
# --------------------------------------------------------------

setwd("~/Courses/Senior/Math Thesis")
# or whatever directory you keep the
# algorithm in

# get the algorithm
source("Cappe2011Online.R")

setwd("~/Courses/Senior/Math Thesis/Data")
# or whatever directory you keep the downloaded
# csv in

# the BIG test (daily for like a billion years)
spy = read.csv("spy.csv")
spy$X.Chg = as.numeric(gsub("%", "", spy$X.Chg)) / 100
spy$X.Chg[is.na(spy$X.Chg)] = 0
spyRet = rev(spy$X.Chg)
dates = as.Date(rev(spy$Exchange.Date[1601:8299]), "%m/%d/%Y")

# this one is not necessary to run for you
tol = seq(0.001, 0.075, by = 0.0001)
traded = cappe2011Algo2(spyRet, tol, lookBack = 100, nStates = 2,
                        bn = c(0.02))

# THIS BLOCK TAKES A LONG TIME TO RUN
# it's just here for my testing but you shouldn't
# need to run it
bns = seq(0.1, 0.49, by = 0.01)
win.graph()
par(mfrow = c(5, 8))
tol = seq(0.001, 0.075, by = 0.0001)
maxes = numeric(length(bn))
tolMax = numeric(length(bn))
for (i in 1:length(bn)){
  temp = cappe2011Algo2(spyRet[1:1500], tol, lookBack = 100, nStates = 2,
                        bn = bns[i])
  plot(temp$finalAlgo ~ tol, type = "l")
  abline(h = temp$finalBnH, col = "blue", lty = 2)
  maxes[i] = max(temp$finalAlgo)
  tolMax[i] = tol[which.max(temp$finalAlgo)]
}
bestTol = tolMax[which.max(maxes)]
bestBn = bns[which.max(maxes)]
# train the model on the first 1500 observations
# to get an idea of the best bn and tol
# we end up with tol = 0.0303 and bn ~ 0.14

# run this one to see the efficiency of the
# algorithm with the best-performing params
tradedBest = cappe2011Algo2(spyRet[1501:length(spyRet)], 
                            tols = c(0.0303),
                            lookBack = 100, bn = c(0.14),
                            nStates = 2)

twoSTest = read.csv("bigSPY2sTest.csv")
traded2STest = data.frame(
  algoPaths = matrix(0, nrow = 6699, ncol = 1))
traded2STest$algoPaths = twoSTest$tol.0.047

win.graph()
plot(tradedBest$algoPaths[1,] ~ dates, type = "l", col = "red",
     main = "Returns over Time", xlab = "Date",
     ylab = "Cumulative Returns")
lines(tradedBest$bnhPath ~ dates, col = "blue")
lines(traded2STest$algoPaths ~ dates, col = "black")
legend("topleft", legend = c("Buy and Hold", "Online",
                             "Primitive"),
       col = c("blue", "red", "black"), lty = c(1, 1))

hist(tradedBest$lowVolMeans - tradedBest$highVolMeans,
     main = "Difference in State Means (l - h)", 
     xlab = "Difference",
     breaks = seq(
       min(tradedBest$lowVolMeans - tradedBest$highVolMeans) - 0.0001,
       max(tradedBest$lowVolMeans - tradedBest$highVolMeans) + 0.0001,
       by = 0.00025
     ))

hist(tradedBest$lowVolSds - tradedBest$highVolSds,
     main = "Difference in State Standard Deviations (l - h)", 
     xlab = "Difference",
     breaks = seq(
       min(tradedBest$lowVolSds - tradedBest$highVolSds) - 0.01, 
       0.001,
       by = 0.0001
     ), xlim = c(-0.003, 0.0001))

plot(tradedBest$highVolMeans,
     main = "High Volatility State Means over Time",
     xlab = "Time", ylab = "Mean", type = "l")

plot(tradedBest$highVolSds,
     main = "High Volatility State Standard Deviations over Time",
     xlab = "Time", ylab = "Mean", type = "l")

primPath = traded2STest$algoPaths[1, ]
onlinePath = tradedBest$algoPaths[1, ]
bnhPath = traded2STest$bnhPath
primCount = 0
primSum = 0
onlineCount = 0
onlineSum = 0
for (i in 2:length(path)){
  primDiff = primPath[i] - primPath[i - 1]
  onlineDiff = onlinePath[i] - onlinePath[i - 1]
  if (primDiff == 0){
    primSum = primSum + (
      (bnhPath[i] - bnhPath[i - 1]) / bnhPath[i - 1])
    primCount = primCount + 1
  }
  if (onlineDiff == 0){
    onlineSum = onlineSum + (
      (bnhPath[i] - bnhPath[i - 1]) / bnhPath[i - 1])
    onlineCount = onlineCount + 1
  }
}
primMean = primSum / primCount
onlineMean = onlineSum / onlineCount
primMean
onlineMean