if(commandArgs(trailingOnly = TRUE)[2]>1)pbo = pbapply::pboptions(type="txt")
# Identify the experiment
lab_experiment <- commandArgs(trailingOnly = TRUE)[1]
# lab_experiment <- "setting1"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))
mle <- qs::qread(file = paste0(path,"mle.qs"))

# Resample and fit

armb <- mle |>
  tidyr::expand_grid(nreps = c(100, 200))
armb_fit <- pbapply::pblapply(purrr::transpose(armb), FUN = function(x){

  cat(paste0("setting n:",x$n, ", p:", x$p, ", R:",x$nreps, ", seed:",x$seed, " | started at ", format(Sys.time(), format = "%F %R")))
  dt <- genData_experiment(PAR = x$theta, N = x$n, P = x$p, SEED = x$seed)
  brn <- .25*x$n
  arm_ctrl <- list(
    MAXT = brn + x$n,
    BURN = brn,
    BATCH = 1,
    STEPSIZE0 = .01,
    PAR1 = 1,
    PAR2 = 1e-4,
    PAR3 = .75,
    PATH_WINDOW = x$n,
    VERBOSE_WINDOW = 10,
    VERBOSE = F,
    SEED = 123
  )

  tune_ctrl <- list(
    MAXT    = .5*x$n,
    MAXA    = 50,
    SCALE   = .8,
    AUTO    = TRUE,
    VERBOSE = FALSE,
    CONV_WINDOW = 1000,
    CONV_TOL = 1e-9
  )

  start <- Sys.time()
  armbt <- armb::armb(
    Y = dt$Y,
    X = dt$X,
    FAMILY = resp,
    NSIM = x$nreps,
    MLE = x$mle,
    TUNE_STEP = TRUE,
    TUNE_GAMMA = FALSE,
    NCORES = 1,
    SEED = 123,
    ARM_CONTROL = arm_ctrl,
    TUNE_CONTROL = tune_ctrl)

  cat(paste0(" | Done in ", round(armbt$time,2), " secs.\n"))

  armbt$batch <- arm_ctrl$BATCH
  return(armbt)
}, cl = commandArgs(trailingOnly = TRUE)[2])
armb$fit <- armb_fit
armb <- armb |>
  dplyr::mutate(
    method = "armb",
    lab = paste0(method, '_R', nreps)
  )

qs::qsave(armb, file = paste0(path,'armb.qs'))
