# This script computes ARM-B for all simulation experiments
if (commandArgs(trailingOnly = TRUE)[2] > 1) {
  pbo = pbapply::pboptions(type = "txt")
}
# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]
# lab_experiment <- "setting1"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))
mle <- qs::qread(file = paste0(path, "mle.qs"))

# Resample and fit
type <- "glm"
if (is.character(resp)) {
  type <- "glmm"
}
armb <- mle |>
  tidyr::expand_grid(
    nreps = c(100, 200),
    step0 = c(
      NA,
      1 /
        2^{
          0:5
        }
    )
  )
armb_fit <- pbapply::pblapply(
  purrr::transpose(armb),
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
      fml <- "cbind(Y, size - Y) ~"
      if (x$p > 0) {
        for (var_idx in 1:x$p) {
          fml <- paste0(fml, " V", var_idx, " + ")
        }
      }
      fml <- paste0(fml, "period")

      step0 <- dplyr::if_else(is.na(x$step0), 1, x$step0)
      arm_ctrl <- list(
        BURN = 1 * x$n,
        STEPSIZE0 = step0,
        TRIM = .25 * x$n,
        VERBOSE = FALSE,
        SEED = 123,
        CONV_CHECK = TRUE,
        CONV_WINDOW = 1,
        TOL = 1 / x$nreps
      )
      tune_ctrl <- list(
        LENGTH = 0,
        BURN = 1 * x$n,
        MAXA = 10,
        SCALE = .5,
        AUTO = FALSE,
        VERBOSE = FALSE
      )
      armbt <- armb::armb4(
        Y = as.numeric(dt$Y),
        X = dt$X,
        FAMILY = binomial(),
        NSIM = x$nreps,
        MLE = x$mle,
        TUNE_STEP = is.na(x$step0),
        NCORES = 1,
        SEED = 123,
        ARM_CONTROL = arm_ctrl,
        TUNE_CONTROL = tune_ctrl,
        VERBOSE = TRUE,
        TYPE = "GLMM",
        FRML = fml
      )
      cat(paste0(" | Done in ", round(armbt$time, 2), " secs.\n"))

      armbt$batch <- arm_ctrl$BATCH
      return(armbt)
    } else {
      dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$seed)

      step0 <- dplyr::if_else(is.na(x$step0), 1, x$step0)
      arm_ctrl <- list(
        BURN = 1 * x$n,
        STEPSIZE0 = step0,
        TRIM = .1 * x$n,
        VERBOSE = FALSE,
        SEED = 123,
        CONV_CHECK = TRUE,
        CONV_WINDOW = 3,
        TOL = 1 / x$nreps
      )
      tune_ctrl <- list(
        LENGTH = 0, #,
        BURN = 1 * x$n,
        MAXA = 10,
        SCALE = .5,
        AUTO = FALSE,
        VERBOSE = FALSE
      )
      armbt <- armb::armb4(
        Y = as.numeric(dt$Y),
        X = dt$X,
        FAMILY = resp,
        NSIM = x$nreps,
        MLE = x$mle,
        TUNE_STEP = is.na(x$step0),
        NCORES = 1,
        SEED = 123,
        ARM_CONTROL = arm_ctrl,
        TUNE_CONTROL = tune_ctrl,
        VERBOSE = TRUE
      )

      cat(paste0(" | Done in ", round(armbt$time, 2), " secs.\n"))

      armbt$batch <- arm_ctrl$BATCH
      return(armbt)
    }
  },
  cl = commandArgs(trailingOnly = TRUE)[2]
)
armb$fit <- armb_fit
armb <- armb |>
  dplyr::mutate(
    method = "armb",
    lab = paste0(method, '_R', nreps)
  )

qs::qsave(armb, file = paste0(path, 'armb.qs'))
