if(commandArgs(trailingOnly = TRUE)[2]>1)pbo = pbapply::pboptions(type="txt")
# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))

# Resample and fit
blb <- sims |> tidyr::expand_grid(
  nreps = c(100, 200),
  gammas = c(.7, .8),
  s = c(1, 2, 3))

blb_fit <- pbapply::pblapply(purrr::transpose(blb), FUN = function(setting){
  try({
    cat(paste0("setting n:",setting$n, ", p:", setting$p, ", R:",setting$nreps, ", gamma:", setting$gammas, ", subsets:", setting$s, ", seed:",setting$seed, " | started at ", format(Sys.time(), format = "%F %R")))
  dt <- genData_experiment(PAR = setting$theta, N = setting$n, P = setting$p, SEED = setting$seed)
  start <- Sys.time()
  blbfit <- rSW2utils::blb(
    data = as.data.frame(cbind(y=dt$Y, dt$X)),
    fun_estimator = function(x, weights) {
      coef(glm(y~.-1, data = x, weights = weights, family = resp))
    },
    fun_metric = function(x) {
      out <- c('sd' = sd(x), diff(quantile(x, probs = c(0.025, 0.975))))
      names(out) <- c('sd', 'ci')
      out
    },
    n_resamples = setting$nreps,
    subset_size_b = round(setting$n^setting$gammas,0),
    window_subsets = setting$s-1,
    n_subsets = setting$s
  )
  out <- list(time = difftime(Sys.time(), start, units = 'secs'), out = blbfit)
  cat(paste0(" | Done in ", round(out$time,2), " secs.\n"))

  return(out)})
}, cl = commandArgs(trailingOnly = TRUE)[2])

blb$fit <- blb_fit
blb <- blb |>
  dplyr::mutate(
    method = 'blb',
    lab = paste0('blb_R', nreps, 'g', gammas*10, 's', s))
qs::qsave(blb, file = paste0(path,'blb.qs'))
