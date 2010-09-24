set.breaks <- function (mat, interval = .1) {
  if (max(abs(mat)) > 1) {
    m <- floor(max(abs(mat)))
  } else {
    m <- round(max(abs(mat)), nchar(1/interval) - 1)
  }

  mm <- m + interval/2
  vals <- seq(interval/2,mm,interval)
  # Note: the first and last values mimic infinity
  mybreaks  <- c(-(m + 1e6), c(-rev( vals ), vals), m + 1e6)
  mybreaks
}

