library(depmixS4)

# no pvals as of now as idk how to tell which points are from
# which cluster without the fitted model object
# GPT suggests using a fixed bn to avoid static regimes
cappe2011Algo2 = function(returns, tols = c(0.005, 0.01, 0.015),
                          lookBack = 30, nStates = 3, bn = 0.15){
  base = returns[1:lookBack]
  nTrades = length(returns) - length(base)
  df = data.frame(base)
  lentols = length(tols)
  
  ## storage
  price = matrix(NA, nrow = lentols, ncol = nTrades)
  decision = matrix(NA, nrow = lentols, ncol = nTrades)
  bnh = numeric(nTrades)
  expVolNext = numeric(nTrades)
  if (nStates == 2){
    lowVolMeans = numeric(nTrades)
    highVolMeans = numeric(nTrades)
    lowVolSds = numeric(nTrades)
    highVolSds = numeric(nTrades)
  } else {
    stateMeans = matrix(NA, nrow = nStates, ncol = nTrades)
    stateSDs = matrix(NA, nrow = nStates, ncol = nTrades)
  }
  
  # begin initialization
  
  # fit initial HMM
  mod = depmix(base ~ 1, data = df,
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
  
  ## extract state means safely
  means = matrix(sapply(1:nStates, function(s) {
    fit@response[[s]][[1]]@parameters$coefficients
  }), nrow = nStates)
  
  # transition matrix
  tmat = matrix(
    as.numeric(getpars(fit)[
      (nStates + 1):(nStates + nStates^2)
    ]), nrow = nStates, byrow = TRUE
  )
  
  pi0 = matrix(
    as.numeric(fit@prior@parameters$coefficients),
    nrow = 1)
  
  # pi0 and matrix both 1xnStates
  # numerator = pi0 * matrix(pnorm(x[1], mean = means, sd = sds), nrow = 1)
  # denominator = sum(numerator)
  gamma0 = matrix( # numerator / denominator
    as.numeric(fit@posterior[nrow(fit@posterior), -1]),
    nrow = 1
  )
  
  rhoq0 = list()
  for (i in 1:nStates){
    rhoq0[[i]] = matrix(0, nrow = nStates, ncol = nStates)
  }
  
  rhog0 = list()
  for (d in 0:2){
    rhog0[[d + 1]] = matrix(0, nrow = nStates, ncol = nStates)
    for (i in 1:nStates){
      for (k in 1:nStates){
        delta = ifelse(i - k == 0, 1, 0)
        rhog0[[d + 1]][i, k] = delta*(base[length(base)]^d)
      } 
    }
  }
  
  gammanprev = gamma0
  rhoqnprev = rhoq0
  rhognprev = rhog0
  Sqn = matrix(0, nrow = nStates, ncol = nStates)
  # what if we didn't update qn?
  # qn = matrix(0, nrow = nStates, ncol = nStates)
  qn = tmat
  Sgn = matrix(0, nrow = 3, ncol = nStates)
  mun = means
  # no idea what v0 is supposed to be so we'll try this
  vn = sds^2
  
  expVoln = sum(gammanprev %*% tmat %*% sqrt(vn))
  
  # end initialization
  
  # GPT fix - define a new list to index properly
  online = returns[(length(base)+1):length(returns)]
  
  # the core updating and trading loop ====================================
  for (n in 1:nTrades){
    # progress checker
    prog = round(n / nTrades, 3) * 100 
    if ((prog %% 5) == 0){ 
      cat(paste(prog, "% complete\n")) 
    }
    
    # trading decision
    first = (n == 1)
    for (i in 1:lentols){
      decision[i, n] = tols[i] >= expVoln
      price[i, n] = ifelse(first, 100 * (1 + decision[i, n]*online[n]),
                         price[i, n - 1] * (1 + decision[i, n]*online[n]))
    }
    bnh[n] = ifelse(first, 100 * (1 + online[n]),
                    bnh[n - 1] * (1 + online[n]))
    
    # update model
    # GPT fix - define gn and calculate gamma properly
    # also use logs and exp for numerical stability
    gn = dnorm(online[n], mean = mun, sd = sqrt(vn))
    pred = gammanprev %*% tmat
    logjoint = log(pred + 1e-12) + t(log(gn + 1e-12))
    joint = exp(logjoint - max(logjoint))
    gamman = joint / sum(joint)
    
    rn = matrix(0, nrow = nStates, ncol = nStates)
    
    for (i in 1:nStates){
      for (j in 1:nStates){
        
        # GPT fix - correct denominator
        denom = sum(gammanprev * tmat[,j])
        
        if (denom < 1e-12) denom = 1e-12
        
        rn[i,j] = gammanprev[i] * tmat[i,j] / denom
      }
    }
    # we needed the old gamman until now, but now we can replace it
    gammanprev = gamman
    
    rhoqn = list()
    for (k in 1:nStates){
      rhoqn[[k]] = matrix(0, nrow = nStates, ncol = nStates)
      for (i in 1:nStates){
        for (j in 1:nStates){
          term1 = bn*ifelse(j - k == 0, 1, 0)*rn[i, j]
          # GPT fix - I was computing term2 wrong due to indexing mistakes
          term2 = (1 - bn) * sum(
            rhoqnprev[[k]][, j] * rn[i, ]
          )
          rhoqn[[k]][i, j] = term1 + term2
        }
      }
    }
    rhoqnprev = rhoqn
    
    rhogn = list()
    for (d in 0:2){
      rhogn[[d + 1]] = matrix(0, nrow = nStates, ncol = nStates)
      for (k in 1:nStates){
        for (i in 1:nStates){
          # GPT fix - logic error in term1
          term1 = bn * rn[i,k] * (online[n]^d)
          # GPT fix - also computing term2 wrong here
          term2 = (1 - bn) * sum(
            rhognprev[[d + 1]][, k] * rn[i, ]
          )
          rhogn[[d + 1]][i, k] = term1 + term2
        }
      }
    }
    rhognprev = rhogn
    
    for (i in 1:nStates){
      for (j in 1:nStates){
        Sqn[i, j] = sum(rhoqn[[i]][j, ] %*% gamman)
      }
    }
    
    for (i in 1:nStates){
      for (j in 1:nStates){
        qn[i, j] = Sqn[i, j]/sum(Sqn[i, ])
      }
    }
    # what if we didn't update the tmat?
    # tmat = qn # GPT fix - I forgot this line originally
    
    for (d in 0:2){
      for (i in 1:nStates){
        Sgn[d + 1, i] = sum(rhogn[[d + 1]][i, ] %*% gamman)
      }
    }
    # GPT fix - add a floor to avoid division by 0
    eps = 1e-8
    for (i in 1:nStates){
      Sgn[1,i] = max(Sgn[1,i], eps)
    }
    
    for (i in 1:nStates){ # the problem is that we end up dividing by 0 here
      mun[i] = Sgn[2, i]/Sgn[1, i]
    }
    
    
    # GPT fix - calculate vn as a vector properly and avoid underflow
    for (i in 1:nStates){
      vn[i, 1] = (Sgn[3,i] - 2*mun[i]*Sgn[2,i] + mun[i]^2*Sgn[1,i]) /
        Sgn[1,i]
      vn[i] = max(vn[i], 1e-8)
    }
    
    if (nStates == 2){
      lowVolMeans[n]  = mun[which.min(vn)]
      highVolMeans[n] = mun[which.max(vn)]
      lowVolSds[n] = sqrt(min(vn))
      highVolSds[n] = sqrt(max(vn))
    } else {
      for (i in 1:nStates){
        stateMeans[i, n] = mun[i]
        stateSDs[i, n] = sqrt(vn[i])
      }
    }
    
    # calculate expected next volatility
    predNext = gamman %*% tmat
    expVoln  = sum(predNext * t(sqrt(vn)))
    expVolNext[n] = expVoln
  }
  
  # the return block
  returnMe = list(
    algoPaths   = price,
    bnhPath     = bnh,
    finalAlgo   = price[, nTrades],
    finalBnH    = bnh[nTrades],
    evNext      = expVolNext,
    meanPos     = rowMeans(decision),
    finalTmat   = tmat
  )
  if (nStates == 2){
    returnMe[['lowVolMeans']] = lowVolMeans
    returnMe[['highVolMeans']] = highVolMeans
    returnMe[['lowVolSds']]   = lowVolSds
    returnMe[['highVolSds']]   = highVolSds
  } else {
    returnMe[['means']] = stateMeans
    returnMe[['sds']] = stateSDs
  }
  
  return(returnMe)
}