pbo = pbapply::pboptions(type="txt")

# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]#"experiment1"#
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))

# Compute the mle via glm
mle <- sims
mle$mle <- pbapply::pblapply(purrr::transpose(sims), FUN = function(x){
  dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$seed)
  out <- glm(y~.-1, family = resp, data = data.frame(y=dt$Y, dt$X))
  return(coef(out))
}, cl = commandArgs(trailingOnly = TRUE)[2])
qs::qsave(mle, file = paste0(path,'mle.qs'))
