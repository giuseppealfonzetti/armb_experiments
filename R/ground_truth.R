# This script computes the parametric bootstrap used as oracle in the simulation experiments
pbo = pbapply::pboptions(type = "txt")

# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))

# type of model
type <- "glm"
if (is.character(resp)) {
  type <- "glmm"
}
# Compute oracle estimates
nrep <- 500
gt <- sims |>
  dplyr::select(-seed) |>
  dplyr::distinct() |>
  tidyr::expand_grid(id = 1:nrep)

gt$est <- pbapply::pblapply(
  purrr::transpose(gt),
  FUN = function(x) {
    if (type == "glmm") {
      dt <- genData_experiment(
        PAR = x$theta,
        N = x$n,
        P = x$p,
        M = x$m,
        SZ = x$sz,
        SEED = x$id
      )
      fml <- "cbind(y, size - y) ~"
      if (x$p > 0) {
        for (var_idx in 1:x$p) {
          fml <- paste0(fml, " V", var_idx, " + ")
        }
      }
      fml <- paste0(fml, "period + (1 | g)")
      # cat(fml, "\n")
      df <- data.frame(y = dt$Y, dt$X)
      # cat(colnames(df),"\n")
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
      # return(fit)
    } else {
      dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$id)
      out <- glm(y ~ . - 1, family = resp, data = data.frame(y = dt$Y, dt$X))
      return(coef(out))
    }
  },
  cl = commandArgs(trailingOnly = TRUE)[2]
)
gt <- gt |>
  dplyr::group_by(n, p) |>
  dplyr::summarise(gt = list(Reduce(rbind, est)))
qs::qsave(gt, file = paste0(path, 'gt.qs'))
