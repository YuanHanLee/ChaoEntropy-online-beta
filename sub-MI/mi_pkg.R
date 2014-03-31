# require(ChaoEntropy)

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



ChaoMI <- function(data, method = c("all", "Chao", "ChaoShen", "Jackknife",
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
  
  method <- match.arg(method)
  if (method == "all") {
    a <- MI_Chao(data, B, conf)
    b <- MI_HT(data, B, conf)
    c <- MI_JK(data, B, conf)
    d <- MI_MLEbc1(data, B, conf)
    e <- MI_MLEbc2(data, B, conf)
    f <- MI_MLE(data, B, conf)
    out <- rbind(a, b, c, d, e, f)
  }
  if (method == "Chao")
    out <- MI_Chao(data, B, conf)
  if (method == "ChaoShen")
    out <- MI_HT(data, B, conf)
  if (method == "Jackknife")
    out <- MI_JK(data, B, conf)
  if (method == "Bias Correct 1") 
    out <- MI_MLEbc1(data, B, conf)
  if (method == "Bias Correct 2")
    out <- MI_MLEbc2(data, B, conf)
  if (method == "Observed")
    out <- MI_MLE(data, B, conf)
  
  return(out)
}
