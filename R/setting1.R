lab_experiment <- 'setting1'
path <- paste0('output/', lab_experiment, '/')
# dir.create(path, recursive = TRUE)

resp <- binomial(link="logit")
# data generation function for the experiment
genData_experiment <- function(PAR, N, P, SEED = 123){
  set.seed(SEED)
  X <- matrix(rt(N*P, df = 3), N, P)
  Y <- armb::simy(FAMILY = resp, ETA = X%*%PAR/sqrt(P))

  return(list('X' = X, 'Y' = Y))
}

save(genData_experiment, file = paste0(path, "setup.rda"))

# settings investigated
sims <- tidyr::expand_grid(n = c(1e3, 2e3), p = c(10, 15), seed = 1:2) |>
  dplyr::mutate(
    theta = purrr::map(p, ~rep(1, .x))
  )

save(sims, resp, genData_experiment, file = paste0(path, "setup.rda"))

