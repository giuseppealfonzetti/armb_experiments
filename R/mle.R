# This script computes the MLE in all simulation experiments.
pbo = pbapply::pboptions(type = "txt")

# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1] #"experiment1"#
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))

# type of model
type <- "glm"
if (is.character(resp)) {
  type <- "glmm"
}
# Compute the mle via glm
mle <- sims
mle$mle <- pbapply::pblapply(
  purrr::transpose(sims),
  FUN = function(x) {
    if (type == "glmm") {
      dt <- genData_experiment(
        PAR = x$theta,
        N = x$n,
        P = x$p,
        M = x$m,
        SZ = x$sz,
        SEED = x$seed
      )
      fml <- "cbind(y, size - y) ~"
      if (x$p > 0) {
        for (var_idx in 1:x$p) {
          fml <- paste0(fml, " V", var_idx, " + ")
        }
      }
      fml <- paste0(fml, "period + (1 | g)")
      df <- data.frame(y = dt$Y, dt$X)
      fit <- lme4::glmer(
        as.formula(fml),
        data = df,
        family = binomial(),
        verbose = 0,
        nAGQ = 9
      )
      return(c(
        summary(fit)$coef[, 1],
        logsd = log(attributes(lme4::VarCorr(fit)$g)$stddev)
      ))
    } else {
      dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$seed)
      out <- glm(y ~ . - 1, family = resp, data = data.frame(y = dt$Y, dt$X))
      return(coef(out))
    }
  },
  cl = commandArgs(trailingOnly = TRUE)[2]
)
qs::qsave(mle, file = paste0(path, 'mle.qs'))
