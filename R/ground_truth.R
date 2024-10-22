pbo = pbapply::pboptions(type="txt")

# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))


# Compute oracle estimates
nrep <- 10
gt <- sims |>
  dplyr::select(-seed) |>
  dplyr::distinct() |>
  tidyr::expand_grid(id = 1:nrep)

gt$est <- pbapply::pblapply(purrr::transpose(gt), FUN = function(x){
  dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$id)
  out <- biglm::bigglm(formula = terms(y~.-1, data = data.frame(y=dt$Y, dt$X)),
                       data = data.frame(y=dt$Y, dt$X),
                       family = resp)
  return(coef(out))
}, cl = commandArgs(trailingOnly = TRUE)[2])
gt <- gt |>
  dplyr::group_by(n, p) |>
  dplyr::summarise(gt = list(Reduce(rbind, est)))
qs::qsave(gt, file = paste0(path,'gt.qs'))
