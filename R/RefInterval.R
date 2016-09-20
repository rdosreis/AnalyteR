#' Analyte distribution plot
#'
#' Produces a histogram plot of the analyte.
#'
#' @param x a numeric vector with the analyte values.
#' @param analyte.name an optional character string giving the analyte name.
#' @param main an overall title for the plot.
#' @param logical; probability if TRUE, probability densities are plotted
#' and a normal density curve is added computed with sample
#' \code{\link[base]{mean}} and standard deviation (\code{\link[base]{sd}});
#' if FALSE, the histogram graphic is a representation of frequencies and
#' no curve is plotted.
#' @return A histogram plot of the analyte.
#' @export AnalyteDistPlot
#' @examples
#' data(ca)
#' AnalyteDistPlot(ca$value)
#' AnalyteDistPlot(ca$value, analyte.name = "Calcium (mg/dL)")
#'
#' data(alt)
#' AnalyteDistPlot(alt$value[alt$sex == "m"],
#' analyte.name = "Alanine Aminotrasferase (U/L)", main = "Male", probability = FALSE)
AnalyteDistPlot <- function(x, analyte.name = "Analyte", main = "", probability = TRUE){
  if (!probability){
    hist(x, col = "grey", border = "white",
         xlab = analyte.name, main = "")
  }else{
    m <- mean(x, na.rm = TRUE)
    sd <- sd(x, na.rm = TRUE)
    hist(x, col = "grey", border = "white", probability = TRUE,
         xlab = analyte.name, main = main)
    lines(x, dnorm(x, mean = m, sd = sd), col = "steelblue", lwd = 2, lty = 2)
  }
}

#' Analyte summary statistics
#'
#' Describes analyte distribution with summary statistics.
#'
#' @param x a numeric vector with the analyte values.
#' @return A vector with the analyte summary statistics: valid sample size, mean, standard deviation
#' range, median, first and third quartiles, and number of missing observations.
#' @export AnalyteSummary
#' @examples
#' data(ca)
#' AnalyteSummary(ca$value)
#'
#' data(alt)
#' AnalyteSummary(alt$value[alt$sex == "f"])
AnalyteSummary <- function(x){
  out  <- c(length(x[!is.na(x)]),
            mean(x, na.rm = TRUE),
            sd(x, na.rm = TRUE),
            range(x, na.rm = TRUE),
            quantile(x, probs = I(1:3)/4, na.rm = TRUE),
            sum(is.na(x)))
  name <- c("N","Mean","SD","Min","Max","Q1","Median","Q2","NAs")
  names(out) <- name
  return(out)
}

#' Analyte reference intervals
#'
#' Produces reference intervals for the analyte.
#'
#' @param data a dataframe.
#' @param x a numeric vector with the analyte values.
#' @param type.ri a character string giving a type for computing the reference interval.
#' This must be one of the strings "all" (default), "non-parametric", "parametric", "bootstrap",
#' "robust".
#' @param conf.ci numeric; the confidence level for computing confidence intervals for the
#' limits of the reference interval. This must be between zero and one.
#' @param B the number of bootstrap replicates for computing confidence intervals for the limits
#' of reference interval, and reference interval type "bootstrap".
#' @param analyte.name an optional character string giving the analyte name.
#' @param partition an optional character string giving the partition variable name. Default is
#' NULL. If partition is not NULL, then reference intervals for each level of the partition
#' are computed.
#' @return Summary stistics, distribution plots, and reference intervals with confidence intervals.
#' @export RefInterval
#' @examples
#' data(ca)
#' RefInterval(data = ca, x = "value", analyte.name = "Calcium (mg/dL)")
#' RefInterval(data = ca, x = "value",
#' type.ri = "non-parametric", analyte.name = "Calcium (mg/dL)",
#' partition = "sex")
#'
#' # Robust ref. intervals for small data sets
#' data(casub)
#' RefInterval(x = "casub",
#' type.ri = "robust", analyte.name = "Calcium (mg/dL)")
#' @references CLSI (...)
RefInterval <- function(data = NULL, x, type.ri = "all", conf.ci = 0.90, B = 1000, analyte.name = "Analyte", partition = NULL){
  if (is.null(data)){
    x <- get(x)
  }else{
    x <- data[, x]
  }
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (length(x[!is.na(x)]) <= 1)
    stop("'x' must have more than one non-NA value")
  na.type <- pmatch(type.ri,
                    c("all", "non-parametric", "parametric",
                      "bootstrap", "robust"))
  if (is.na(na.type))
    stop("invalid 'type.ri' argument")
  if (conf.ci < 0 | conf.ci > 1)
    stop("'conf.ci' must be between zero and one")
  # if (!is.integer(B))
  #   stop("'B' must be integer")
  # B <- as.integer(B)
  if (B <= 1)
    stop("'B' must be greater than one")
  conf <- conf.ci
  out <- list()
  out$combined <- list()
  aux <- suppressWarnings(AnalyteSummary(x))
  aux[1] <- format(round(aux[1], 0), nsmall = 0)
  aux[2:8] <- format(round(as.numeric(aux[2:8]), 2), nsmall = 2)
  out$combined$summaries <- data.frame(value = t(t(aux)))
  out$combined$RefInterval <-  ComputeRefInterval(x, type.ri = type.ri, conf = conf, B = B)
  AnalyteDistPlot(x, analyte.name = analyte.name, main = "Combined")
  # ---------------------------------------------
  # Partition
  # ---------------------------------------------
  if (!is.null(partition)){
    label.partition <- partition
    partition <- factor(data[, partition])
    partition.levels <- unique(partition)
    plot(x ~ partition, xlab = label.partition, ylab = analyte.name, col = "grey", pch = 16)
    for (i in 1:length(partition.levels)){
      pli <- as.character(partition.levels[i])
      x.part <- x[partition == partition.levels[i]]
      out[[pli]] <- list()
      aux <- suppressWarnings(AnalyteSummary(x.part))
      aux[1] <- format(round(aux[1], 0), nsmall = 0)
      aux[2:8] <- format(round(as.numeric(aux[2:8]), 2), nsmall = 2)
      out[[pli]]$summaries <- data.frame(value = t(t(aux)))
      out[[pli]]$RefInterval <-  ComputeRefInterval(x.part, type.ri = type.ri, conf = conf, B = B)
      AnalyteDistPlot(x.part, analyte.name = analyte.name, main = pli)
    }
  }
  return(out)
}
