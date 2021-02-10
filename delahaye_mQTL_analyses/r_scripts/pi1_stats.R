require(splines)
# assess replication of p-values in study 1 from study 2
# i.e. the proportion of p
pi1 <- function(p, m) {
  p <- sort(p)
  vals <- data.frame(lambda = seq(0, 0.95, 0.01))
  resp <- unlist(lapply(vals$lambda, function(x) sum(p > x) / (m * (1 - x))))
  vals$pi_all <- resp
  fit <- lm(pi_all ~ ns(lambda, df = 3), data = vals)
  plot(vals$lambda, vals$pi_all)
  lines(vals$lambda, p1 <- predict(fit, vals))
  return(1 - predict(fit, data.frame(lambda = 1)))
}
