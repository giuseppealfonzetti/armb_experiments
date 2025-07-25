lab_experiment <- 'bin_id'
path <- paste0('output/', lab_experiment, '/')
# dir.create(path, recursive = TRUE)

# data generation function for the experiment
genData_experiment <- function(PAR, N, P, M, SZ = 20, SEED = 123) {
  beta0 <- PAR[1]
  beta_ext <- PAR[2:(P + 1)] / sqrt(P)
  beta_g <- c(0, PAR[(P + 2):(P + M)])
  logsd <- PAR[length(PAR)]
  set.seed(SEED)
  X <- as.data.frame(matrix(rnorm(N * (P)), N, P))
  raneff <- rnorm(N, 0, exp(logsd))
  mat <- cbind(X, g = 1:N)
  df <- cbind(mat[sort(rep(1:N, M)), ], period = rep(1:M, N), size = SZ)

  linpred <- sapply(1:nrow(df), function(r) {
    fixed <- beta0
    if (P > 0) {
      fixed <- fixed + t(as.numeric((X[df[r, "g"], ]))) %*% beta_ext
    }
    fixed + raneff[df[r, "g"]] + beta_g[df[r, "period"]]
  })

  MU <- sapply(linpred, plogis)
  Y <- sapply(MU, function(mui) rbinom(1, SZ, prob = mui))

  df$period <- factor(df$period)

  return(list('X' = df, 'Y' = Y))
}

save(genData_experiment, file = paste0(path, "setup.rda"))

# settings investigated
sims_setup <- qs::qread("output/sims_ri_setup.qs")
sims <- sims_setup |>
  dplyr::mutate(
    theta = purrr::map2(
      p,
      m,
      ~ c(rep(1, .x + 1), seq(-1, 0, length.out = (.y - 1)), 0)
    )
  )

resp <- "LOGRI"
save(sims, resp, genData_experiment, file = paste0(path, "setup.rda"))
