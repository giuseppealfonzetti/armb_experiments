lab_experiment <- 'ber_eq_fxd'
path <- paste0('output/', lab_experiment, '/')
# dir.create(path, recursive = TRUE)

resp <- binomial(link="logit")
# data generation function for the experiment
genData_experiment <- function(PAR, N, P, SEED = 123){
  set.seed(SEED)
  S <- matrix(.2, P,P)
  diag(S) <- 1
  X <- mvtnorm::rmvnorm(N, mean = rep(0,P), sigma = S)
  MU <- resp$linkinv(X%*%PAR/sqrt(P))
  Y <- sapply(MU, function(mui) rbinom(1, 1, prob = mui))

  return(list('X' = X, 'Y' = Y))
}

save(genData_experiment, file = paste0(path, "setup.rda"))

# settings investigated
sims_setup <- qs::qread("output/sims_setup.qs")
sims <- sims_setup |>
  dplyr::mutate(
    theta = purrr::map(p, ~rep(1, .x))
  )

save(sims, resp, genData_experiment, file = paste0(path, "setup.rda"))
