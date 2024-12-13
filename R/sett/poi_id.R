lab_experiment <- 'poi_id'
path <- paste0('output/', lab_experiment, '/')
# dir.create(path, recursive = TRUE)

resp <- poisson(link="log")
# data generation function for the experiment
genData_experiment <- function(PAR, N, P, SEED = 123){
  set.seed(SEED)
  X <- matrix(rnorm(N*P), N, P)
  MU <- resp$linkinv(X%*%PAR/sqrt(P))
  Y <- sapply(MU, function(mui) rpois(1, lambda = mui))

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
