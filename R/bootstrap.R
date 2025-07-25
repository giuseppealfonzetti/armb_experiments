# This script compute nonparametric bootstrap in all simulation experiments
if (commandArgs(trailingOnly = TRUE)[2] > 1) {
  pbo = pbapply::pboptions(type = "txt")
}
# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))
mle <- qs::qread(file = paste0(path, "mle.qs"))

# Fitting function to pass to boot
glmFun <- function(D, IDX, START) {
  mod <- glm(y ~ . - 1, family = resp, data = D[IDX, ], start = START)
  out <- coef(mod)
  return(out)
}

glmerFun <- function(D, IDX, START, FML) {
  dmat <- purrr:::reduce(D[IDX], dplyr::bind_rows)
  mod <- lme4::glmer(as.formula(FML), family = binomial(), data = dmat)
  out <- return(c(
    summary(mod)$coef[, 1],
    logsd = log(attributes(lme4::VarCorr(mod)$g)$stddev)
  ))
  return(out)
}
# type of model
type <- "glm"
if (is.character(resp)) {
  type <- "glmm"
}
# Resample and fit
btstrp <- mle |>
  tidyr::expand_grid(nreps = c(100, 200))
fit <- pbapply::pblapply(
  purrr::transpose(btstrp),
  FUN = function(x) {
    cat(paste0(
      "setting n:",
      x$n,
      ", p:",
      x$p,
      ", R:",
      x$nreps,
      ", seed:",
      x$seed,
      " | started at ",
      format(Sys.time(), format = "%F %R")
    ))
    if (type == "glmm") {
      dt <- genData_experiment(
        PAR = x$theta,
        N = x$n,
        P = x$p,
        M = x$m,
        SZ = x$sz,
        SEED = x$seed
      )
      df <- data.frame(y = dt$Y, dt$X)
      list_df <- lapply(by(df, df$g, identity), as.data.frame)

      fml <- "cbind(y, size - y) ~"
      if (x$p > 0) {
        for (var_idx in 1:x$p) {
          fml <- paste0(fml, " V", var_idx, " + ")
        }
      }
      fml <- paste0(fml, "period + (1 | g)")

      start <- Sys.time()
      bt <- boot::boot(
        data = list_df,
        statistic = glmerFun,
        R = x$nreps,
        parallel = 'multicore',
        ncpus = 1,
        START = x$mle,
        FML = fml
      )
      out <- list(
        time = difftime(Sys.time(), start, units = 'secs'),
        pars = bt[['t']]
      )
      cat(paste0(" | Completed after ", round(out$time, 2), " secs.\n"))

      return(out)
    } else {
      dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$seed)
      start <- Sys.time()
      bt <- boot::boot(
        data = data.frame(y = dt$Y, dt$X),
        statistic = glmFun,
        R = x$nreps,
        parallel = 'multicore',
        ncpus = 1,
        START = x$mle
      )
      out <- list(
        time = difftime(Sys.time(), start, units = 'secs'),
        pars = bt[['t']]
      )
      cat(paste0(" | Completed after ", round(out$time, 2), " secs.\n"))

      return(out)
    }
  },
  cl = commandArgs(trailingOnly = TRUE)[2]
)


btstrp$fit <- fit
btstrp <- btstrp |>
  dplyr::mutate(
    method = "bootstrap",
    lab = paste0(method, '_R', nreps)
  )

qs::qsave(btstrp, file = paste0(path, 'bt.qs'))
