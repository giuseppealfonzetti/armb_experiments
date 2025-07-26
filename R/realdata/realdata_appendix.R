# This script benchmarks ARM-B on three additional datasets
library(tidyverse) |> suppressPackageStartupMessages()
library(evtree) |> suppressPackageStartupMessages()


# region LIKE
mblrdata <- as_tibble(base::readRDS(url(
  "https://slcladal.github.io/data/mbd.rda",
  "rb"
)))
data.like <- mblrdata |>
  dplyr::mutate_if(is.character, factor) |>
  dplyr::mutate(Age = relevel(Age, "Young")) |>
  dplyr::arrange(ID) |>
  dplyr::mutate(across(where(is.factor), ~ as.numeric(.x) - 1)) |>
  dplyr::mutate(ID = ID + 1, size = 1) |>
  dplyr::rename(g = ID)

n <- length(unique(data.like$g))
fml.base <- "(Gender + Age + ConversationType + Priming)"
fml.glm <- paste0("SUFlike ~ ", fml.base)
fml.glmm <- paste0(fml.glm, "+ (1 | g)")
mle.fit <- lme4::glmer(
  fml.glmm,
  data = data.like,
  family = binomial
)
mle <- c(
  summary(mle.fit)$coef[, 1],
  logsd = log(attributes(lme4::VarCorr(mle.fit)$g)$stddev)
)
glmerFun <- function(D, IDX, START, FML) {
  dmat <- purrr:::reduce(D[IDX], dplyr::bind_rows)
  mod <- lme4::glmer(as.formula(FML), family = binomial(), data = dmat)
  out <- return(
    c(
      summary(mod)$coef[, 1],
      logsd = try(log(attributes(lme4::VarCorr(mod)$g)$stddev))
    ) *
      sqrt(nrow(dmat) / nrow(data.like))
  )
  return(out)
}

list_df <- lapply(by(data.like, data.like$g, identity), as.data.frame)
tictoc::tic()
bt <- boot::boot(
  data = list_df,
  statistic = glmerFun,
  R = 100,
  parallel = 'multicore',
  ncpus = 1,
  START = mle,
  FML = fml.glmm
)
btsd <- apply(bt[["t"]], 2, sd)
bt.like.time <- tictoc::toc()

step0 <- 1
arm_ctrl <- list(
  BURN = 1 * n,
  STEPSIZE0 = step0,
  TRIM = .25 * n,
  VERBOSE = FALSE,
  SEED = 123,
  CONV_CHECK = TRUE,
  CONV_WINDOW = 1,
  TOL = 1e-2
)
tune_ctrl <- list(
  LENGTH = 1 * n,
  BURN = 0,
  MAXA = 10,
  SCALE = .5,
  AUTO = TRUE,
  VERBOSE = FALSE
)
nreps <- 100
tictoc::tic()
armbt <- armb::armb4(
  Y = as.numeric(data.like$SUFlike),
  X = data.like[, -6],
  FAMILY = binomial(link = "logit"),
  NSIM = nreps,
  MLE = mle,
  TUNE_STEP = TRUE,
  NCORES = 1,
  SEED = 123,
  ARM_CONTROL = arm_ctrl,
  TUNE_CONTROL = tune_ctrl,
  VERBOSE = FALSE,
  TYPE = "GLMM",
  FRML = paste0("cbind(Y, size - Y) ~ ", fml.base)
)
armbsd <- apply(armbt$pars, 2, sd, na.rm = TRUE)
armb.like.time <- tictoc::toc()

gg_track <- map_dfr(
  1:length(armbt$chains),
  ~ {
    tibble(
      chain_id = .x,
      iter = armbt$chains[[.x]]$path_iters,
      deltaL2 = armbt$chains[[.x]]$path_norm,
      avdeltaL2 = sapply(armbt$chains[[.x]]$path_delta, function(x) {
        sqrt(sum(x^2))
      }),
      diffLInf = armbt$chains[[.x]]$path_diff,
      grLInf = armbt$chains[[.x]]$path_ngr,
      nll = armbt$chains[[.x]]$path_nll,
      pnll = armbt$chains[[.x]]$path_pf,
    )
  }
) |>
  pivot_longer(
    c(deltaL2, avdeltaL2, grLInf, diffLInf, nll, pnll),
    names_to = "metric",
    values_to = "val"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("deltaL2", "avdeltaL2", "diffLInf", "grLInf", "nll", "pnll"),
      labels = c(
        latex2exp::TeX(
          "$||\\Delta_t||_2$"
        ),
        latex2exp::TeX("$||\\bar{\\Delta}_t||_2$"),
        latex2exp::TeX("$||\\Delta_t-\\Delta_{t-1}||_{\\infty}$"),
        latex2exp::TeX("$||(\\Delta_t-\\Delta_{t-1})/\\eta_t||_{\\infty}$"),
        latex2exp::TeX("$f_t$"),
        latex2exp::TeX("$|f_t-f_{t-1}|/|f_{t-1}|$")
      )
    )
  ) |>
  # filter(iter > .1 * n) |>
  ggplot(aes(x = iter, y = val, col = as.factor(chain_id))) +
  facet_wrap(~metric, labeller = "label_parsed", scales = "free") +
  geom_line(alpha = .5) +
  # geom_hline(
  #   data = tibble(levs = c(1e-1, 1e-2, 1e-3, 1e-4)),
  #   aes(yintercept = levs),
  #   linetype = "dashed",
  #   color = "lightgrey"
  # ) +
  # geom_vline(aes(xintercept = fit$burn), color = "skyblue") +
  labs(x = "Iteration", y = "") +
  theme_bw() +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_viridis_d()
# gg_track

res.like <- tibble(
  var_lab = names(mle),
  mle = mle,
  sd_btstrp = btsd,
  sd_armb = armbsd,
  sd_blb = NULL
)
logsd <- res.like[res.like$var_lab == "logsd.(Intercept)", ]
newsd <- logsd |>
  mutate(
    var_lab = "sd",
    mle = exp(mle),
    sd_btstrp = sd_btstrp * mle,
    sd_armb = sd_armb * mle,
  )

gg.like <- res.like |>
  bind_rows(newsd) |>
  filter(var_lab != "logsd.(Intercept)") |>
  select(var_lab, mle, sd_btstrp, sd_armb) |>
  pivot_longer(
    cols = starts_with("sd"),
    names_to = "method",
    values_to = "sd"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("sd_armb", "sd_btstrp", "sd_blb"),
      labels = c("ARM-B", "Oracle", "BLB")
    )
  ) |>
  ggplot(aes(x = mle, y = var_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "lightgrey") +
  geom_segment(
    aes(x = mle - 1.96 * sd, xend = mle + 1.96 * sd, col = method),
    position = position_dodge(width = 0.2)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), #remove y axis labels
    axis.ticks.y = element_blank(), #remove y axis ticks
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  scale_color_grey(start = .2, end = .8) +
  labs(
    x = "MLE",
    y = "Parameters",
    col = "",
    subtitle = "95% Confidence Intervals"
  )
#gg.like
ggsave(gg.like, filename = "output/like.pdf", width = 10, height = 6)

# endregion
# region MAGIC
data("MAGICGammaTelescope", package = "evtree")
data.magic <- as_tibble(MAGICGammaTelescope) |>
  mutate_at("class", ~ as.numeric(factor(.x)) - 1) |>
  mutate(across(-class, ~ scale(.)[, 1]))

n <- nrow(data.magic)
p <- ncol(data.magic)
fam = binomial(link = "logit")
frml <- "class~(fLength + fWidth + fSize +  fConc + fConc1 +  fAsym + fM3Long + fM3Trans + fAlpha + fDist)"
mle <- glm(formula = frml, data = data.magic, family = fam)

glmFun <- function(D, IDX, START) {
  mod <- glm(frml, family = fam, data = D[IDX, ], start = START)
  out <- coef(mod)
  return(out)
}

tictoc::tic()
bt <- boot::boot(
  data = data.magic,
  statistic = glmFun,
  R = 100,
  parallel = 'multicore',
  ncpus = 1,
  START = coef(mle)
)
btsd <- apply(bt[["t"]], 2, sd)
bt.magic.time <- tictoc::toc()

tictoc::tic()
blbfit <- rSW2utils::blb(
  data = as.data.frame(cbind(class = data.magic$class, model.matrix(mle))),
  fun_estimator = function(x, weights, start) {
    coef(glm(
      as.formula(frml),
      data = x,
      weights = weights,
      family = fam,
      start = coef(mle)
    ))
  },
  fun_metric = function(x) {
    out <- c('sd' = sd(x))
    names(out) <- c('sd')
    out
  },
  n_resamples = 100,
  subset_size_b = round(n^.7, 0),
  window_subsets = 2,
  n_subsets = 3
)
blb.magic.time <- tictoc::toc()


step0 <- 1
arm_ctrl <- list(
  BURN = 1 * n,
  STEPSIZE0 = step0,
  TRIM = .25 * n,
  VERBOSE = FALSE,
  SEED = 123,
  CONV_CHECK = TRUE,
  CONV_WINDOW = 1,
  TOL = 1e-2
)
tune_ctrl <- list(
  LENGTH = 1 * n,
  BURN = 0,
  MAXA = 10,
  SCALE = .5,
  AUTO = TRUE,
  VERBOSE = FALSE
)
nreps <- 100
tictoc::tic()
armbt <- armb::armb4(
  Y = data.magic$class,
  X = model.matrix(mle),
  FAMILY = fam,
  NSIM = nreps,
  MLE = coef(mle),
  TUNE_STEP = TRUE,
  NCORES = 1,
  SEED = 123,
  ARM_CONTROL = arm_ctrl,
  TUNE_CONTROL = tune_ctrl,
  VERBOSE = FALSE
)
armbsd <- apply(armbt$pars, 2, sd)
armb.magic.time <- tictoc::toc()

gg_track <- map_dfr(
  1:length(armbt$chains),
  ~ {
    tibble(
      chain_id = .x,
      iter = armbt$chains[[.x]]$path_iters,
      deltaL2 = armbt$chains[[.x]]$path_norm,
      avdeltaL2 = sapply(armbt$chains[[.x]]$path_delta, function(x) {
        sqrt(sum(x^2))
      }),
      diffLInf = armbt$chains[[.x]]$path_diff,
      grLInf = armbt$chains[[.x]]$path_ngr,
      nll = armbt$chains[[.x]]$path_nll,
      pnll = armbt$chains[[.x]]$path_pf,
    )
  }
) |>
  pivot_longer(
    c(deltaL2, avdeltaL2, grLInf, diffLInf, nll, pnll),
    names_to = "metric",
    values_to = "val"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("deltaL2", "avdeltaL2", "diffLInf", "grLInf", "nll", "pnll"),
      labels = c(
        latex2exp::TeX(
          "$||\\Delta_t||_2$"
        ),
        latex2exp::TeX("$||\\bar{\\Delta}_t||_2$"),
        latex2exp::TeX("$||\\Delta_t-\\Delta_{t-1}||_{\\infty}$"),
        latex2exp::TeX("$||(\\Delta_t-\\Delta_{t-1})/\\eta_t||_{\\infty}$"),
        latex2exp::TeX("$f_t$"),
        latex2exp::TeX("$|f_t-f_{t-1}|/|f_{t-1}|$")
      )
    )
  ) |>
  filter(iter > .1 * n) |>
  ggplot(aes(x = iter, y = val, col = as.factor(chain_id))) +
  facet_wrap(~metric, labeller = "label_parsed", scales = "free") +
  geom_line(alpha = .5) +
  labs(x = "Iteration", y = "") +
  theme_bw() +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_viridis_d()
# gg_track
res.magic <- tibble(
  var_lab = names(coef(mle)),
  mle = coef(mle),
  sd_btstrp = btsd,
  sd_armb = armbsd,
  sd_blb = as.numeric(blbfit)
)

gg.magic <- res.magic |>
  select(var_lab, mle, sd_btstrp, sd_armb, sd_blb) |>
  pivot_longer(
    cols = starts_with("sd"),
    names_to = "method",
    values_to = "sd"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("sd_armb", "sd_btstrp", "sd_blb"),
      labels = c("ARM-B", "Oracle", "BLB")
    )
  ) |>
  ggplot(aes(x = mle, y = var_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "lightgrey") +
  geom_segment(
    aes(x = mle - 1.96 * sd, xend = mle + 1.96 * sd, col = method),
    position = position_dodge(width = 0.2)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), #remove y axis labels
    axis.ticks.y = element_blank(), #remove y axis ticks
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  scale_color_grey(start = .8, end = .2) +
  labs(
    x = "MLE",
    y = "Regression coefficients",
    col = "",
    subtitle = "95% Confidence Intervals"
  )
# gg.magic
ggsave(gg.magic, filename = "output/magic.pdf", width = 10, height = 6)

# endregion
# region WILT

wilt.train <- read.csv("data/wilt/training.csv")
wilt.test <- read.csv("data/wilt/testing.csv")
data.wilt <- as_tibble(wilt.train) |>
  bind_rows(wilt.test) |>
  mutate_at("class", ~ as.numeric(factor(.x)) - 1) |>
  mutate(across(-class, ~ scale(.)[, 1]))

n <- nrow(data.wilt)
p <- ncol(data.wilt)
fam = binomial(link = "logit")
frml <- "class~(GLCM_pan + Mean_Green + Mean_Red + Mean_NIR + SD_pan)-1"
mle <- glm(formula = frml, data = data.wilt, family = fam)

glmFun <- function(D, IDX, START) {
  mod <- glm(frml, family = fam, data = D[IDX, ], start = START)
  out <- coef(mod)
  return(out)
}

tictoc::tic()
bt <- boot::boot(
  data = data.wilt,
  statistic = glmFun,
  R = 100,
  parallel = 'multicore',
  ncpus = 1,
  START = coef(mle)
)
btsd <- apply(bt[["t"]], 2, sd)
bt.wilt.time <- tictoc::toc()

tictoc::tic()
blbfit <- rSW2utils::blb(
  data = data.wilt,
  fun_estimator = function(x, weights, start) {
    coef(glm(
      class ~ . - 1,
      data = x,
      weights = weights,
      family = fam,
      start = coef(mle)
    ))
  },
  fun_metric = function(x) {
    out <- c('sd' = sd(x))
    names(out) <- c('sd')
    out
  },
  n_resamples = 100,
  subset_size_b = round(n^.8, 0),
  window_subsets = 2,
  n_subsets = 3
)
blb.wilt.time <- tictoc::toc()


step0 <- 5
arm_ctrl <- list(
  BURN = 1 * n,
  STEPSIZE0 = step0,
  TRIM = .1 * n,
  VERBOSE = FALSE,
  SEED = 123,
  CONV_CHECK = TRUE,
  CONV_WINDOW = 3,
  TOL = 1e-2
)
tune_ctrl <- list(
  LENGTH = 1 * n,
  BURN = 0,
  MAXA = 10,
  SCALE = .5,
  AUTO = TRUE,
  VERBOSE = FALSE
)
nreps <- 100
tictoc::tic()
armbt <- armb::armb4(
  Y = data.wilt$class,
  X = model.matrix(mle),
  FAMILY = fam,
  NSIM = nreps,
  MLE = coef(mle),
  TUNE_STEP = TRUE,
  NCORES = 1,
  SEED = 123,
  ARM_CONTROL = arm_ctrl,
  TUNE_CONTROL = tune_ctrl,
  VERBOSE = FALSE
)
armbsd <- apply(armbt$pars, 2, sd)
armb.wilt.time <- tictoc::toc()


gg_track <- map_dfr(
  1:length(armbt$chains),
  ~ {
    tibble(
      chain_id = .x,
      iter = armbt$chains[[.x]]$path_iters,
      deltaL2 = armbt$chains[[.x]]$path_norm,
      avdeltaL2 = sapply(armbt$chains[[.x]]$path_delta, function(x) {
        sqrt(sum(x^2))
      }),
      diffLInf = armbt$chains[[.x]]$path_diff,
      grLInf = armbt$chains[[.x]]$path_ngr,
      nll = armbt$chains[[.x]]$path_nll,
      pnll = armbt$chains[[.x]]$path_pf,
    )
  }
) |>
  pivot_longer(
    c(deltaL2, avdeltaL2, grLInf, diffLInf, nll, pnll),
    names_to = "metric",
    values_to = "val"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c("deltaL2", "avdeltaL2", "diffLInf", "grLInf", "nll", "pnll"),
      labels = c(
        latex2exp::TeX(
          "$||\\Delta_t||_2$"
        ),
        latex2exp::TeX("$||\\bar{\\Delta}_t||_2$"),
        latex2exp::TeX("$||\\Delta_t-\\Delta_{t-1}||_{\\infty}$"),
        latex2exp::TeX("$||(\\Delta_t-\\Delta_{t-1})/\\eta_t||_{\\infty}$"),
        latex2exp::TeX("$f_t$"),
        latex2exp::TeX("$|f_t-f_{t-1}|/|f_{t-1}|$")
      )
    )
  ) |>
  filter(iter > .1 * n) |>
  ggplot(aes(x = iter, y = val, col = as.factor(chain_id))) +
  facet_wrap(~metric, labeller = "label_parsed", scales = "free") +
  geom_line(alpha = .5) +
  labs(x = "Iteration", y = "") +
  theme_bw() +
  scale_y_log10() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  ) +
  scale_color_viridis_d()
gg_track

res.wilt <- tibble(
  var_lab = names(coef(mle)),
  mle = coef(mle),
  sd_btstrp = btsd,
  sd_armb = armbsd,
  sd_blb = as.numeric(blbfit)
)

gg.wilt <- res.wilt |>
  select(var_lab, mle, sd_btstrp, sd_armb, sd_blb) |>
  pivot_longer(
    cols = starts_with("sd"),
    names_to = "method",
    values_to = "sd"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("sd_armb", "sd_btstrp", "sd_blb"),
      labels = c("ARM-B", "Oracle", "BLB")
    )
  ) |>
  ggplot(aes(x = mle, y = var_lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "lightgrey") +
  geom_segment(
    aes(x = mle - 1.96 * sd, xend = mle + 1.96 * sd, col = method),
    position = position_dodge(width = 0.2)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_blank(), #remove y axis labels
    axis.ticks.y = element_blank(), #remove y axis ticks
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  scale_color_grey(start = .8, end = .2) +
  labs(
    x = "MLE",
    y = "Regression coefficients",
    col = "",
    subtitle = "95% Confidence Intervals"
  )
# gg.wilt
ggsave(gg.wilt, filename = "output/wilt.pdf", width = 10, height = 6)

# endregion
