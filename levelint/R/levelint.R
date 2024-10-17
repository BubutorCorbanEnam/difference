#' Check the Level of intergration and Check for Stationarity
#'
#' This function makes a time series stationary by using the Augmented Dickey-Fuller (ADF),
#' Phillips-Perron (PP), and KPSS tests.
#'
#' @param timeseries A numeric vector or time series object.
#' @param alpha Significance level for the tests (e.g., 0.05).
#'
#' @return A differenced time series that has been made stationary.
#' @examples
#' library(tseries)
#' library(urca)
#' library(levelint)
#' data <- cumsum(rnorm(100)) # Random walk
#' levelint(data, alpha = 0.05)
#'
#' @export
levelint <- function(timeseries, alpha = 0.05) {
  if (!requireNamespace("tseries", quietly = TRUE)) {
    stop("Package 'tseries' is required but not installed.")
  }
  if (!requireNamespace("urca", quietly = TRUE)) {
    stop("Package 'urca' is required but not installed.")
  }

  df <- timeseries

  # ADF Test
  differences <- 0
  result <- tseries::adf.test(df)
  p_value <- result$p.value
  while (p_value >= alpha) {
    differences <- differences + 1
    df <- diff(df, differences = differences)
    result <- tseries::adf.test(df)
    p_value <- result$p.value
    cat(sprintf("ADF test: P-value %f at difference %d.\n", p_value, differences))
  }

  # PP Test
  df <- timeseries  # Reset df for PP test
  differences <- 0
  result <- urca::ur.pp(df, type = "Z-tau", model = "constant", lags = "short")
  p_value <- result@teststat[1]
  critical_value <- result@cval[1]

  while (p_value >= critical_value) {
    differences <- differences + 1
    df <- diff(df, differences = differences)
    result <- urca::ur.pp(df, type = "Z-tau", model = "constant", lags = "short")
    p_value <- result@teststat[1]
    critical_value <- result@cval[1]
    cat(sprintf("PP test: Test statistic %f at difference %d.\n", p_value, differences))
  }

  # KPSS Test
  df <- timeseries  # Reset df for KPSS test
  differences <- 0
  result <- tseries::kpss.test(df, null = "Trend")
  p_value <- result$p.value

  while (p_value <= alpha) {
    differences <- differences + 1
    df <- diff(df, differences = differences)
    result <- tseries::kpss.test(df, null = "Trend")
    p_value <- result$p.value
    cat(sprintf("KPSS test: P-value %f at difference %d.\n", p_value, differences))
  }

  return(head(df, 5))  # Return the first five rows of the stationary time series
}
