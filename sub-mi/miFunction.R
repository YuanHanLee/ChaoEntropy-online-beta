ChaoEntropy.Ind <-
  function(data, method = c("all", "Chao", "ChaoShen", "Grassberger", 
                            "Jackknife", "Zhang", "Observed"), 
           B = 200, conf = 0.95, se = TRUE) {
    method <- match.arg(method)
    if (method == "all") {
      a <- entropy_Chao(data, B, conf)
      b <- entropy_ChaoShen(data, B, conf)
      c <- entropy_Grassberger(data, B, conf)
      d <- entropy_Jackknife(data, B, conf)
      e <- entropy_Zhang(data, B, conf)
      f <- entropy_Observed(data, B, conf)
      out <- rbind(a, b, c, d, e, f)
    }
    if (method == "Chao")
      out <- entropy_Chao(data, B, conf)
    if (method == "ChaoShen")
      out <- entropy_ChaoShen(data, B, conf)
    if (method == "Grassberger")
      out <- entropy_Grassberger(data, B, conf)
    if (method == "Jackknife") 
      out <- entropy_Jackknife(data, B, conf)
    if (method == "Zhang")
      out <- entropy_Zhang(data, B, conf)
    if (method == "Observed")
      out <- entropy_Observed(data, B, conf)
    
    if (se == FALSE) 
      out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
    
    return(out)
  }

ChaoEntropy.Sam <-
  function(data, method = c("all", "Chao", "Observed"),
           B = 200, conf = 0.95, se = TRUE) {
    method <- match.arg(method)
    if (method == "all") {
      a <- entropy_ChaoIn(data, B, conf)
      f <- entropy_ObservedIn(data, B, conf)
      out <- rbind(a, f)
    }
    if (method == "Chao")
      out <- entropy_ChaoIn(data, B, conf)
    if (method == "Observed")
      out <- entropy_ObservedIn(data, B, conf)
    if (se == FALSE) 
      out <- data.frame(Estimator = out[, 1], row.names = rownames(out))
    
    return(out)
  }

ChaoEntropyMI <-
  function(data, datatype = c("abundance", "incidence"), 
           method = c("all", "Chao", "ChaoShen", "Grassberger", 
                      "Jackknife", "Zhang", "Observed"), 
           se = TRUE, nboot = 200, conf = 0.95) {
    if (is.matrix(data) == TRUE || is.data.frame(data) == TRUE) {
      data <- as.matrix(data)
      if (ncol(data) != 1 & nrow(data) != 1)
        stop("Error: The data format is wrong.")
      if (ncol(data) == 1) {
        data <- data[, 1]
      } else {
        data <- data[1, ]
      } 
    } 
    
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
      #       cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
      #           "\n")
      #       cat("          We use \"conf\" = 0.95 to calculate!", 
      #           "\n\n")
      conf <- 0.95
    }
    
    if (se == TRUE) {
      B <- nboot
      if (nboot < 1)
        nboot <- 1
      #       if (nboot == 1)
      #         cat("Warning: When \"nboot\" =" ,B, ", the bootstrap s.e. and confidence interval can't be calculated.", 
      #             "\n\n")  
    }
    if (se == FALSE)
      nboot <- 1
    method <- match.arg(method)
    datatype <- match.arg(datatype)
    if (datatype == "abundance") {
      #       if (sum(data > 0) == 1) {
      #         cat("Warning: When the individual-based (abundance) data only have \"ONE\" species.", 
      #             "\n")
      #         cat("         ALL estimator are equal to 0,and the standard error will meaningless.", 
      #             "\n\n")
      #       }
      if (sum(data) == 0)
        stop("Error: The data didn't have enough information.")
      out <- ChaoEntropy.Ind(data, method, nboot, conf, se)
    }
    if (datatype == "incidence") {
      if (sum(data[1] < data[-1]) != 0)
        stop("Error: Total number of sampling units should be greater than the species incidence frequency.")
      if (length(data) == 1) 
        stop("Error: The input format of first entry should be total number of sampling units, and followed by species incidence frequency. You only type the total number of sampling units.")
      if (sum(data[-1]) == 0)
        stop("Error: The data didn't have enough information.")
      if (method[1] == "ChaoShen" || method == "Grassberger" || 
            method == "Jackknife" || method == "Zhang") {
        #         cat("Warning: Sample-based incidence data doesn't have this estimator.",
        #             "\n")
        #         cat("         We use Chao and Observed estimator below.", "\n\n")
        method <- "all"
      }
      out <- ChaoEntropy.Sam(data, method, nboot, conf, se)
    }
    return(out)
  }

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
  Hx <-ChaoEntropyMI(x, method="Jackknife", se=F)
  Hy <- ChaoEntropyMI(y, method="Jackknife", se=F)
  Hxy <- ChaoEntropyMI(xy, method="Jackknife", se=F)
  I.hat <- as.numeric(Hx + Hy - Hxy)
  return(I.hat)
  
}

EstHTFun <- function(mat) {
  n <- sum(mat)
  x <- apply(mat, 1, sum)
  y <- apply(mat, 2, sum)
  xy <- as.numeric(mat)
  Hx <-ChaoEntropyMI(x, method="ChaoShen", se=F)
  Hy <- ChaoEntropyMI(y, method="ChaoShen", se=F)
  Hxy <- ChaoEntropyMI(xy, method="ChaoShen", se=F)
  I.hat <- as.numeric(Hx + Hy - Hxy)
  return(I.hat)
}

EstMEEFun <- function(mat) {
  n <- sum(mat)
  x <- apply(mat, 1, sum)
  y <- apply(mat, 2, sum)
  xy <- as.numeric(mat)
  Hx <-ChaoEntropyMI(x, method="Chao", se=F)
  Hy <- ChaoEntropyMI(y, method="Chao", se=F)
  Hxy <- ChaoEntropyMI(xy, method="Chao", se=F)
  I.hat <- as.numeric(Hx + Hy - Hxy)
  return(I.hat)
}

MI_Chao <- function(mat, B=200, conf=0.95) {
  mydata <- as.matrix(mat)
  est <- EstMEEFun(mydata)
  se <- MIBootstrapProposedFun(mydata, B, EstMEEFun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Chao_MI (2013)")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}

MI_MLE <- function(mat, B=200, conf=0.95) {
  mydata <- as.matrix(mat)
  est <- EstMLEFun(mydata)
  se <- MIBootstrapProposedFun(mydata, B, EstMLEFun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Observed_MI")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}

MI_MLEbc1 <- function(mat, B=200, conf=0.95) {
  mydata <- as.matrix(mat)
  est <- EstMLEbc2Fun(mydata)
  se <- MIBootstrapProposedFun(mydata, B, EstMLEbc2Fun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Bias Correct 1")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}

MI_MLEbc2 <- function(mat, B=200, conf=0.95) {
  mydata <- as.matrix(mat)
  est <- EstMLEbc4Fun(mydata)
  se <- MIBootstrapProposedFun(mydata, B, EstMLEbc4Fun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Bias Correct 2")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}

MI_JK <- function(mat, B=200, conf=0.95) {
  mydata <- as.matrix(mat)
  est <- EstJKFun(mydata)
  se <- MIBootstrapProposedFun(mydata, B, EstJKFun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Zahl (1977) Jackknife")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}

MI_HT <- function(mat, B=200, conf=0.95) {
  mydata <- as.matrix(mat)
  est <- EstHTFun(mydata)
  se <- MIBootstrapProposedFun(mydata, B, EstHTFun)
  z <- qnorm(1 - (1 - conf)/2)
  CI <- c(max(est - z * se, 0), est + z * se)
  out <- matrix(c(est, se, CI), nrow = 1)
  rownames(out) <- c("Chao_Shen (2003)")
  colnames(out) <- c("Estimator", "Bootstrap s.e.",
                     paste(conf*100, "% Lower"), paste(conf*100, "% Upper"))
  return(out)
}



ChaoMI <- function(data, method = c("Chao", "ChaoShen", "Jackknife",
                                    "Bias Correct 1", "Bias Correct 2", 
                                    "Observed"), 
                   se=TRUE, nboot=200, conf=0.95) {
  data <- as.matrix(data)
  
  if (is.numeric(conf) == FALSE || conf > 1 || conf < 0) {
    cat("Warning: \"conf\"(confidence level) must be a numerical value between 0 and 1, e.g. 0.95.",
        "\n")
    cat("          We use \"conf\" = 0.95 to calculate!", 
        "\n\n")
    conf <- 0.95
  }
  if (se == TRUE) {
    B <- nboot
    if (nboot < 1)
      nboot <- 1
    if (nboot == 1)
      cat("Warning: When \"nboot\" =" ,B, ", the bootstrap s.e. and confidence interval can't be calculated.", 
          "\n\n")  
  }
  if (se == FALSE)
    nboot <- 1
  
  if ("Chao" %in% method) {
    a <- MI_Chao(data, B, conf)
  } else {
    a <- NULL
  }
  if ("ChaoShen" %in% method) {
    b <- MI_HT(data, B, conf)
  } else {
    b <- NULL
  }
  if ("Jackknife" %in% method) {
    c <- MI_JK(data, B, conf)
  } else {
    c <- NULL
  }
  if ("Bias Correct 1" %in% method) {
    d <- MI_MLEbc1(data, B, conf)
  } else {
    d <- NULL
  }
  if ("Bias Correct 2" %in% method) {
    e <- MI_MLEbc2(data, B, conf)
  } else {
    e <- NULL
  }
  if ("Observed" %in% method) {
    f <- MI_MLE(data, B, conf)
  } else {
    f <- NULL
  }
  out <- rbind(a, b, c, d, e, f)
  return(out)
}
