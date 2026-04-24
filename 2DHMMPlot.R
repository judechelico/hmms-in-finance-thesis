source("GPTCorrectedGMMReplication.R")
par(mfrow = c(1, 1))

set.seed(1)

sig1 = crossprod(
  matrix(rnorm(4, mean = 0, sd = 2), 2)
)
sig2 = crossprod(
  matrix(rnorm(4, mean = 0, sd = 2), 2)
)
sig3 = crossprod(
  matrix(rnorm(4, mean = 0, sd = 2), 2)
)

mu1 = rnorm(2, 0, 7)
mu2 = rnorm(2, 0, 7)
mu3 = rnorm(2, 0, 7)

s1 = mvrnorm(10, mu1, sig1)
s2 = mvrnorm(15, mu2, sig2)
s3 = mvrnorm(12, mu3, sig3)

X = c(s1[1:5, 1], s2[1:10, 1], s1[6:10, 1], s3[1:7, 1], 
      s2[11:15, 1], s3[8:12, 1])
Y = c(s1[1:5, 2], s2[1:10, 2], s1[6:10, 2], s3[1:7, 2], 
      s2[11:15, 2], s3[8:12, 2])
col = c(rep("red", 5), rep("blue", 10), rep("red", 5),
        rep("green", 7), rep("blue", 5), rep("green", 5))
pch = c(rep(15, 5), rep(16, 10), rep(15, 5),
        rep(17, 7), rep(16, 5), rep(17, 5))

len = length(X)
plot(Y[2:(len - 1)] ~ X[2:(len - 1)], 
     col = col[2:(len - 1)], pch = pch[2:(len - 1)],
     xlim = c(min(X) - 1, max(X) + 1.75),
     ylim = c(min(Y) - 1, max(Y) + 0.05),
     xlab = "Simulated x-values",
     ylab = "Simulated y-values",
     main = "2D HMGMM Visualization")

points(ellipse(sig1, center = mu1), 
       type = 'l', col = 'red', lty = 2, lwd = 2)

points(ellipse(sig2, center = mu2), 
       type = 'l', col = 'blue', lty = 2, lwd = 2)

points(ellipse(sig3, center = mu3), 
       type = 'l', col = 'green', lty = 2, lwd = 2)

for (n in 2:length(X)){
  lines(c(Y[n - 1], Y[n]) ~ c(X[n - 1], X[n]))
}

points(Y[2:(len - 1)] ~ X[2:(len - 1)], 
       col = col[2:(len - 1)], pch = pch[2:(len - 1)])
points(c(Y[1], Y[length(Y)]) ~ c(X[1], X[length(X)]),
       col = "purple", bg = "purple", pch = 23)