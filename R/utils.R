# Compute a nonparametric reference interval 
RefIntervalNP <- function(x){
  ri <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  return(ri)
}

# Compute a parametric reference interval
RefIntervalP <- function(x){
  m <- mean(x, na.rm = TRUE)
  sd <- sd(x, na.rm = TRUE)
  ri <- m + c(-1, 1) * qnorm(p = 0.975) * sd
  return(ri)
}

# Compute a bootstrap reference interval
RefIntervalBoot <- function(x, B = 1000, ci = FALSE, conf = 0.90){
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  f.quantile <- function(x, ind, ...){
    quantile(x[ind], ...)
  }
  b.quant <- boot::boot(x, f.quantile, R = B, probs = c(0.025, 0.975), na.rm = TRUE)
  if (ci){
    ri <- apply(b.quant$t, MARGIN = 2, FUN = quantile, probs = c(1/2, probs))
  }else{
    ri <- apply(b.quant$t, MARGIN = 2, FUN = quantile, probs = 1/2)
  }
  return(ri)
}

# Compute a robust reference interval
RefIntervalRob <- function(x, c1 = 3.7, c2 = 205.6, eps = 1e-05){
  median <- median(x, na.rm = TRUE)
  mad <- mad(x, center = median, constant = 0.6745, na.rm = TRUE)
  n <- sum(!is.na(x))
  T_b <- median
  repeat{
    T_b_old <- T_b
    u <- (x - T_b) / (c1 * mad)
    w <- ifelse(abs(u) < 1 , (1 - u^2)^2, 0)
    T_b <- sum(w * x) / sum(w)
    change <- (T_b - T_b_old)
    if (change < eps) break
  }
  u <- (x - median) / (c2 * mad)
  sb_c2 <- c2 * mad * sqrt(
    (n * sum( (u^2 * (1 - u^2)^4)[abs(u) < 1] ) / 
       ( sum( ((1 - u^2) * (1 - 5 * u^2))[abs(u) < 1] ) * 
           max( c(1, - 1 + sum( ((1 - u^2) * (1 - 5 * u^2))[abs(u) < 1] ) ) ) ) ) )
  u <- (x - median) / (c1 * mad)
  sb_c1 <- c1 * mad * sqrt(
    (n * sum( (u^2 * (1 - u^2)^4)[abs(u) < 1] ) / 
       ( sum( ((1 - u^2) * (1 - 5 * u^2))[abs(u) < 1] ) * 
           max( c(1, - 1 + sum( ((1 - u^2) * (1 - 5 * u^2))[abs(u) < 1] ) ) ) ) ) )
  u <- (x - T_b) / (c1 * sb_c1)
  ST_c1 <- c1 * sb_c1 * sqrt(
    sum( (u^2 * (1 - u^2)^4)[abs(u) < 1] / 
           ( sum( ((1 - u^2) * (1 - 5 * u^2))[abs(u) < 1] ) * 
               max( c(1, - 1 + sum( ((1 - u^2) * (1 - 5 * u^2))[abs(u) < 1] ) ) ) ) ) )
  ri <- T_b + c(-1, 1) * qt(p = 0.975, df = (n - 1)) * sqrt(sb_c2^2 + ST_c1^2)
  return(ri)
}

# Compute bootstrap confidence intervals
# for the limits of the reference interval
RI_ConfidenceInterval <- function(x, B = 1000, conf = 0.90, method = "non-parametric"){
  # method = c("non-parametric", "parametric", "robust")
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  RefIntervalNP.boot <- function(x, ind, ...){
    RefIntervalNP(x[ind], ...)
  }
  RefIntervalP.boot <- function(x, ind, ...){
    RefIntervalP(x[ind], ...)
  }
  RefIntervalBoot.boot <- function(x, ind, ...){
    RefIntervalBoot(x[ind], ...)
  }
  RefIntervalRob.boot <- function(x, ind, ...){
    RefIntervalRob(x[ind], ...)
  }
  if (method == "non-parametric"){
    b.quant <- boot::boot(x, RefIntervalNP.boot, R = B)
    ci <- apply(b.quant$t, MARGIN = 2, FUN = quantile, probs = probs)
  }else{
    if (method == "parametric"){
      b.quant <- boot::boot(x, RefIntervalP.boot, R = B)
      ci <- apply(b.quant$t, MARGIN = 2, FUN = quantile, probs = probs)
    }else{
      b.quant <- boot::boot(x, RefIntervalRob.boot, R = B)
      ci <- apply(b.quant$t, MARGIN = 2, FUN = quantile, probs = probs)
    }
  }
  return(ci)
}

# Combines reference interval with their respective CIs
ComputeRefInterval <- function(x, type.ri = "all", conf = 0.90, B = 1000){
  # type.ri = c("all", "non-parametric", "parametric", "bootstrap", "robust")
  if (type.ri == "all"){
    ri <- list()
    # Non-parametric
    ri.df <- as.data.frame(
      matrix(0, nrow = 2, ncol = 4,
             dimnames = list(NULL, c("Fraction", "Reference Limits",
                                     paste("Lower Limit", round(conf, 2), "CI"),
                                     paste("Upper Limit", round(conf, 2), "CI")))))
    ri.df[,1] <- c("2.5%", "97.5%")
    ri.df[,2] <- RefIntervalNP(x)
    ri.df[,3:4] <- t(RI_ConfidenceInterval(x = x, method = "non-parametric", B = B, conf = conf))
    ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
    ri$nonparametric <- ri.df
    # Parametric
    ri.df <- as.data.frame(
      matrix(0, nrow = 2, ncol = 4,
             dimnames = list(NULL, c("Fraction", "Reference Limits",
                                     paste("Lower Limit", round(conf, 2), "CI"),
                                     paste("Upper Limit", round(conf, 2), "CI")))))
    ri.df[,1] <- c("2.5%", "97.5%")
    ri.df[,2] <- RefIntervalP(x)
    ri.df[,3:4] <- t(RI_ConfidenceInterval(x = x, method = "parametric", B = B, conf = conf))
    ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
    ri$parametric <- ri.df
    # Bootstrap
    ri.df <- as.data.frame(
      matrix(0, nrow = 2, ncol = 4,
             dimnames = list(NULL, c("Fraction", "Reference Limits",
                                     paste("Lower Limit", round(conf, 2), "CI"),
                                     paste("Upper Limit", round(conf, 2), "CI")))))
    ri.df[,1] <- c("2.5%", "97.5%")
    aux <- RefIntervalBoot(x = x, ci = TRUE, conf = conf, B = B)
    ri.df[,2] <- aux[1,]
    ri.df[,3:4] <- t(aux[2:3,])
    rm(aux)
    ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
    ri$bootstrap <- ri.df
    # Robust
    ri.df <- as.data.frame(
      matrix(0, nrow = 2, ncol = 4,
             dimnames = list(NULL, c("Fraction", "Reference Limits",
                                     paste("Lower Limit", round(conf, 2), "CI"),
                                     paste("Upper Limit", round(conf, 2), "CI")))))
    ri.df[,1] <- c("2.5%", "97.5%")
    ri.df[,2] <- RefIntervalRob(x)
    ri.df[,3:4] <- t(RI_ConfidenceInterval(x = x, method = "robust", B = B, conf = conf))
    ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
    ri$robust <- ri.df
  }else{
    if (type.ri == "non-parametric"){
      ri <- list()
      # Non-parametric
      ri.df <- as.data.frame(
        matrix(0, nrow = 2, ncol = 4,
               dimnames = list(NULL, c("Fraction", "Reference Limits",
                                       paste("Lower Limit", round(conf, 2), "CI"),
                                       paste("Upper Limit", round(conf, 2), "CI")))))
      ri.df[,1] <- c("2.5%", "97.5%")
      ri.df[,2] <- RefIntervalNP(x)
      ri.df[,3:4] <- t(RI_ConfidenceInterval(x = x, method = "non-parametric", B = B, conf = conf))
      ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
      ri$nonparametric <- ri.df
    }else{
      if (type.ri == "parametric"){
        ri <- list()
        # Parametric
        ri.df <- as.data.frame(
          matrix(0, nrow = 2, ncol = 4,
                 dimnames = list(NULL, c("Fraction", "Reference Limits",
                                         paste("Lower Limit", round(conf, 2), "CI"),
                                         paste("Upper Limit", round(conf, 2), "CI")))))
        ri.df[,1] <- c("2.5%", "97.5%")
        ri.df[,2] <- RefIntervalP(x)
        ri.df[,3:4] <- t(RI_ConfidenceInterval(x = x, method = "parametric", B = B, conf = conf))
        ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
        ri$parametric <- ri.df
      }else{
        if (type.ri == "bootstrap"){
          ri <- list()
          # Bootstrap
          ri.df <- as.data.frame(
            matrix(0, nrow = 2, ncol = 4,
                   dimnames = list(NULL, c("Fraction", "Reference Limits",
                                           paste("Lower Limit", round(conf, 2), "CI"),
                                           paste("Upper Limit", round(conf, 2), "CI")))))
          ri.df[,1] <- c("2.5%", "97.5%")
          aux <- RefIntervalBoot(x = x, ci = TRUE, conf = conf, B = B)
          ri.df[,2] <- aux[1,]
          ri.df[,3:4] <- t(aux[2:3,])
          rm(aux)
          ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
          ri$bootstrap <- ri.df
        }else{
          ri <- list()
          # Robust
          ri.df <- as.data.frame(
            matrix(0, nrow = 2, ncol = 4,
                   dimnames = list(NULL, c("Fraction", "Reference Limits",
                                           paste("Lower Limit", round(conf, 2), "CI"),
                                           paste("Upper Limit", round(conf, 2), "CI")))))
          ri.df[,1] <- c("2.5%", "97.5%")
          ri.df[,2] <- RefIntervalRob(x)
          ri.df[,3:4] <- t(RI_ConfidenceInterval(x = x, method = "robust", B = B, conf = conf))
          ri.df[,2:4] <- format(round(ri.df[,2:4], 3), nsmall = 3)
          ri$robust <- ri.df
        }
      }
    }
    
  }
  return(ri)
}
