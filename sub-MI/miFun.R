# require(ChaoEntropy)


#######################################################################
#                                                                     #
#                        Preparation Function                         #
#                                                                     #
#######################################################################
entropyFun <- function(p) {
  p <- p[p > 0]
  out <- -sum(p * log(p))
  return(out)
}

f0Fun <- function(x) {
  n <- sum(x)
  x <- x[x > 0]
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (f2 > 0) 
    f0 <- (n-1) / n * f1^2 / (2 * f2)
  if (f2 == 0)
    f0 <- (n-1) / n * f1 * (f1 - 1) / 2
  f0 <- ceiling(f0)
  return(f0)
}

Chao1Fun <- function(x) {
  n <- sum(x)
  x <- x[x > 0]
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (f2 > 0) 
    f0 <- (n-1) / n * f1^2 / (2 * f2)
  if (f2 == 0)
    f0 <- (n-1) / n * f1 * (f1 - 1) / 2
  Shat <- sum(x > 0) + f0
  return(Shat)
}

Candf0Fun <- function(f1, f2, n) {
  if (f2 > 0) {
    C <- 1 - f1 / n * ((n - 1) * f1 / ((n - 1) * f1 + 2 * f2))
    f0 <- (n - 1) / n * f1^2 / (2 * f2)
  } else if (f2 == 0 & f1 != 0) {
    C <- 1 - f1 / n * ((n - 1) * (f1 - 1) / ((n - 1) * (f1 - 1) + 2))
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
  } else {
    C <- 1
    f0 <- (n - 1) / n * f1 * (f1 - 1) / 2
  }
  f0 <- ceiling(f0)
  return(c(C, f0))
}

MIBootstrapFun_MLE <- function(mat, B, FunName) {
  n <- sum(mat)
  prob.hat <- mat / n
  W <- rmultinom(B, n, as.numeric(prob.hat))
  
  se <- sd(apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(mat))
    FunName(w1)
  }))
  return(se)
}

MIBootstrapProposedFun <- function(mat, B, FunName) {
  x <- as.numeric(mat)
  n <- sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  tmp <- Candf0Fun(f1, f2, n)
  Chat <- tmp[1]
  f0 <- tmp[2]
  
  f0x <- f0Fun(apply(mat, 1, sum))
  f0y <- f0Fun(apply(mat, 2, sum))
  blank <- f0x * ncol(mat) + f0y * nrow(mat) + f0x * f0y
  if (f0 > blank) f0 <- blank
  if (f0 != 0) {
    lambda <- (1 - Chat) / sum(x / n * (1 - x / n)^n)
    pi <- x / n * (1 - lambda * (1 - x /n)^n)
    
    extend <- c(pi, rep(0, f0x * ncol(mat) + f0y * nrow(mat) + f0x * f0y))
    zeroPos <- which(extend == 0)
    pos <- sample(zeroPos, f0)
    extend[pos] <- (1 - Chat) / f0
    
    set.seed(123)
    W1 <- rmultinom(B, n, extend)
    
    se <- sd(apply(W1, 2, function(w) {
      tab1 <- matrix(w[1:length(x)], ncol=ncol(mat))
      remain <- w[-(1:length(x))]
      
      if (f0y == 0) {
        part1 <- NA
      } else {
        part1 <- remain[1:(f0y * nrow(mat))]
      }
      
      add1 <- matrix(part1, ncol=f0y, nrow=nrow(mat))
      tab2 <- cbind(tab1, add1)
      
      if (f0y == 0) {
        part2 <- remain
      } else {
        part2 <- remain[-(1:(f0y * nrow(mat)))]
      }
      
      
      add2 <- matrix(part2, ncol=(ncol(mat) + f0y))
      AddTable <- rbind(tab2, add2)
      FunName(AddTable)
    }))
  } else {
    se <- MIBootstrapFun_MLE(mat, B, FunName)
  }
  
  return(se)
}


#######################################################################
#                                                                     #
#                         Estimation Function                         #
#                                                                     #
#######################################################################
EstMLEFun <- function(mat) {  # MLE
  n <- sum(mat)
  prob.hat <- mat / n
  px.hat <- apply(prob.hat, 1, sum)
  py.hat <- apply(prob.hat, 2, sum)
  I.hat <- entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
  # MLE of Mutual Information!
  return(I.hat)
}

EstMLEbc2Fun <- function(mat) {
  n <- sum(mat)
  prob.hat <- mat / n
  px.hat <- apply(prob.hat, 1, sum)
  py.hat <- apply(prob.hat, 2, sum)
  I.hat <- entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
  # MLE of Mutual Information!
  
  x <- Chao1Fun(apply(mat, 1, sum))
  y <- Chao1Fun(apply(mat, 2, sum))
  xy <- Chao1Fun(mat)
  bc <- (xy - x - y + 1) / (2 * n)
  
  I.hat.bc <- I.hat - bc
  return(I.hat.bc)
}

EstMLEbc4Fun <- function(mat) {
  n <- sum(mat)
  prob.hat <- mat / n
  px.hat <- apply(prob.hat, 1, sum)
  py.hat <- apply(prob.hat, 2, sum)
  I.hat <- entropyFun(px.hat) + entropyFun(py.hat) - entropyFun(prob.hat)
  # MLE of Mutual Information!
  
  x <- Chao1Fun(apply(mat, 1, sum))
  y <- Chao1Fun(apply(mat, 2, sum))
  xy <- Chao1Fun(mat)
  bc <- (xy - x - y + 1) / (2 * n)
  
  temp1 <- prob.hat[prob.hat > 0]
  temp2 <- px.hat[px.hat > 0]
  temp3 <- py.hat[py.hat > 0]
  bc2 <- (sum(1/temp2) + sum(1/temp3) - sum(1/temp1) - 1) / (12 * n^2)
  I.hat.bc <- I.hat - bc + bc2
  return(I.hat.bc)
}

EstJKFun <- function(mat) {
  n <- sum(mat)
  x <- apply(mat, 1, sum)
  y <- apply(mat, 2, sum)
  xy <- as.numeric(mat)
  Hx <-ChaoEntropy(x, method="Jackknife", se=F)
  Hy <- ChaoEntropy(y, method="Jackknife", se=F)
  Hxy <- ChaoEntropy(xy, method="Jackknife", se=F)
  I.hat <- as.numeric(Hx + Hy - Hxy)
  return(I.hat)
  
}

EstHTFun <- function(mat) {
  n <- sum(mat)
  x <- apply(mat, 1, sum)
  y <- apply(mat, 2, sum)
  xy <- as.numeric(mat)
  Hx <-ChaoEntropy(x, method="ChaoShen", se=F)
  Hy <- ChaoEntropy(y, method="ChaoShen", se=F)
  Hxy <- ChaoEntropy(xy, method="ChaoShen", se=F)
  I.hat <- as.numeric(Hx + Hy - Hxy)
  return(I.hat)
}

EstMEEFun <- function(mat) {
  n <- sum(mat)
  x <- apply(mat, 1, sum)
  y <- apply(mat, 2, sum)
  xy <- as.numeric(mat)
  Hx <-ChaoEntropy(x, method="Chao", se=F)
  Hy <- ChaoEntropy(y, method="Chao", se=F)
  Hxy <- ChaoEntropy(xy, method="Chao", se=F)
  I.hat <- as.numeric(Hx + Hy - Hxy)
  return(I.hat)
}

XY2PmatFun <- function(X, Y) {
  X <- round(X)
  Y <- round(Y)
  
  i <- rep(unique(X), length(unique(Y)))
  j <- rep(unique(Y), each=length(unique(X)))
  
  pxy <- mapply(function(i, j) sum(X == i & Y == j) / length(X), i, j)
  pmat <- matrix(pxy, ncol=length(unique(Y)))
  return(pmat)
}

#######################################################################
#                                                                     #
#                        Simulation Function                          #
#                                                                     #
#######################################################################
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


MI.Data.Pic.Fun <- function(X, Y) {
  df1 <- data.frame(X, Y)
  p1 <- ggplot(df1, aes(x=X, y=Y)) + geom_point()  
  df2 <- round(df1)
  p2 <- ggplot(df2, aes(x=X, y=Y)) + geom_point()
  multiplot(p1, p2, cols=2)
}

MI.Pmat.I.Fun <- function(pmat, rd=5) {
  pxy <- as.numeric(pmat)
  px <- apply(pmat, 1, sum)
  py <- apply(pmat, 2, sum)
  
  Dx <- sum(px > 0)
  Dy <- sum(py > 0)
  D <- sum(as.numeric(pmat) > 0)
  
  Hx <- entropyFun(px)
  Hy <- entropyFun(py)
  Hxy <- entropyFun(as.numeric(pmat))
  I <- Hx + Hy - Hxy
  return(round(I, rd))
}

Plot.4.Sample.Size.Fun <- function(temp, N) {
  bias <- as.numeric(temp)
  index <- rep(N, each=nrow(temp))
  method <- rep(c("MLE", "MLE-bc", "MLE-bc2", "Jackknife", "Chao_Shen", "Chao (2013)"), ncol(temp))
  df <- data.frame(method, index, bias)
  
  p <- ggplot(df, aes(x=index, y=bias, group=method))
  p + geom_hline(yintercept=0, lty=2) +   
    geom_line(aes(colour=method),size=1) + 
    geom_point(aes(colour=method), size=3)
}

MI.Pmat.Fun <- function(pmat, R, n, B, rd=3) {
  pxy <- as.numeric(pmat)
  px <- apply(pmat, 1, sum)
  py <- apply(pmat, 2, sum)
  
  Dx <- sum(px > 0)
  Dy <- sum(py > 0)
  D <- sum(as.numeric(pmat) > 0)
  
  Hx <- entropyFun(px)
  Hy <- entropyFun(py)
  Hxy <- entropyFun(as.numeric(pmat))
  I <- Hx + Hy - Hxy
  
  temp1 <- pxy[pxy > 0]
  temp2 <- px[px > 0]
  temp3 <- py[py > 0]
  bc1 <- (D - Dx - Dy + 1) / (2 * n)
  bc2 <- bc1 - (sum(1/temp2) + sum(1/temp3) - sum(1/temp1) - 1) / (12 * n^2)
  thmbc <- round(c(bc1, bc2), 5)
  
  #   set.seed(12345)
  W <- rmultinom(R, n, as.numeric(pmat))
  
  MLEtemp <- apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(pmat))
    est <- EstMLEFun(w1)
    se <- MIBootstrapProposedFun(w1, B, FunName=EstMLEFun)
    return(c(est, se))
  })
  MLE <- mean(MLEtemp[1, ])
  MLEse <- sd(MLEtemp[1, ])
  MLEese <- mean(MLEtemp[2, ])
  MLErmse <- sqrt(mean((MLEtemp[1, ] - I)^2))
  
  
  MLEbc2temp <- apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(pmat))
    est <- EstMLEbc2Fun(w1)
    se <- MIBootstrapProposedFun(w1, B, FunName=EstMLEbc2Fun)
    return(c(est, se))
  })
  MLEbc2 <- mean(MLEbc2temp[1, ])
  MLEbc2se <- sd(MLEbc2temp[1, ])
  MLEbc2ese <- mean(MLEbc2temp[2, ])
  MLEbc2rmse <- sqrt(mean((MLEbc2temp[1, ] - I)^2))
  
  MLEbc4temp <- apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(pmat))
    est <- EstMLEbc4Fun(w1)
    se <- MIBootstrapProposedFun(w1, B, FunName=EstMLEbc4Fun)
    return(c(est, se))  
  })
  MLEbc4 <- mean(MLEbc4temp[1, ])
  MLEbc4se <- sd(MLEbc4temp[1, ])
  MLEbc4ese <- mean(MLEbc4temp[2, ])
  MLEbc4rmse <- sqrt(mean((MLEbc4temp[1, ] - I)^2))
  
  JKtemp <- apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(pmat))
    est <- EstJKFun(w1)
    se <- MIBootstrapProposedFun(w1, B, FunName=EstJKFun)
    return(c(est, se))  
  })
  JK <- mean(JKtemp[1, ])
  JKse <- sd(JKtemp[1, ])
  JKese <- mean(JKtemp[2, ])
  JKrmse <- sqrt(mean((JKtemp[1, ] - I)^2))
  
  HTtemp <- apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(pmat))
    est <- EstHTFun(w1)
    se <- MIBootstrapProposedFun(w1, B, FunName=EstHTFun)
    return(c(est, se))  
  })
  HT <- mean(HTtemp[1, ])
  HTse <- sd(HTtemp[1, ])
  HTese <- mean(HTtemp[2, ])
  HTrmse <- sqrt(mean((HTtemp[1, ] - I)^2))
  
  MEEtemp <- apply(W, 2, function(w) {
    w1 <- matrix(w, ncol=ncol(pmat))
    est <- EstMEEFun(w1)
    se <- MIBootstrapProposedFun(w1, B, FunName=EstMEEFun)
    return(c(est, se))
  })
  MEE <- mean(MEEtemp[1, ])
  MEEse <- sd(MEEtemp[1, ])
  MEEese <- mean(MEEtemp[2, ])
  MEErmse <- sqrt(mean((MEEtemp[1, ] - I)^2))
  
  est <- c(MLE, MLEbc2, MLEbc4, JK, HT, MEE)
  bias <- est - I
  se <- c(MLEse, MLEbc2se, MLEbc4se, JKse, HTse, MEEse)
  ese <- c(MLEese, MLEbc2ese, MLEbc4ese, JKese, HTese, MEEese)
  rmse <- c(MLErmse, MLEbc2rmse, MLEbc4rmse, JKrmse, HTrmse, MEErmse)
  output <- round(data.frame(est, bias, se, ese, rmse), rd)
  rownames(output) <- c("MLE", "MLE-bc", "MLE-bc2", "Jackknife", "Chao_Shen", "Chao (2013)")
  colnames(output) <- c("Avg Estimator", "Avg Bias", "Sample sd", "Avg est sd", "RMSE")
  return(list(Mutual_Information=I, output=output, thmBC=thmbc))
}


PlotFun <- function(dataName, cutby, col="deepskyblue1", howbig=7) {
  subdata <- new.fs1[which(new.fs1$sp == dataName), ]
  
  gx <- subdata$gx
  gy <- subdata$gy
  
  k <- seq(0, 500-cutby, by=cutby)
  h <- seq(cutby, 500, by=cutby)
  mat <- mapply(function(k, h) {
    i <- seq(0, 500-cutby, by=cutby)
    j <- seq(cutby, 500, by=cutby)
    mapply(function(i, j) sum(gx >= k & gx < h & gy >= i & gy < j), i, j)
  }, k, h)
  
  df <- data.frame(gx, gy)
  
  x <- rep(seq(cutby/2, 500-cutby/2, by=cutby), each=500/cutby)
  y <- rep(seq(cutby/2, 500-cutby/2, by=cutby), 500/cutby)
  ds <- data.frame(x, y, label=as.numeric(mat))
  
  dl <- data.frame(h=seq(cutby, 500-cutby, by=cutby), 
                   v=seq(cutby, 500-cutby, by=cutby))
  
  p <- ggplot(df, aes(x=gx, y=gy))
  pic <- p + geom_point(colour=col, alpha = I(0.5)) + 
    geom_text(data=ds, aes(x=x, y=y, label=label), fontface=2, size=howbig) + 
    labs(title = paste(dataName, "Scatter Plot")) + 
    geom_hline(data=dl, aes(yintercept = h), lty=2) + 
    geom_vline(data=dl, aes(xintercept = v), lty=2)
  pmat <- mat / sum(mat)
  px <- apply(pmat, 1, sum)
  py <- apply(pmat, 2, sum)
  
  Dx <- sum(px > 0)
  Dy <- sum(py > 0)
  D <- sum(as.numeric(pmat) > 0)
  
  Hx <- entropyFun(px)
  Hy <- entropyFun(py)
  Hxy <- entropyFun(as.numeric(pmat))
  I <- Hx + Hy - Hxy
  
  return(list(pic, I, pmat))
}
