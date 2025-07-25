# This script generates all the plots related to the simulations experiments
# for logistic regression, poisson regression and binomial regression with random intercepts

library(tidyverse) |> suppressPackageStartupMessages()
library(latex2exp) |> suppressPackageStartupMessages()
library(ggh4x) |> suppressPackageStartupMessages()
library(ggpubr) |> suppressPackageStartupMessages()

# region: LOGISTIC | INDEPENDENT #####
# Identify the experiment
lab_experiment <- "ber_id"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path, "gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS) {
  apply(PARS, 2, function(x) {
    c('sd' = sd(x))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path, "armb.qs"))
armb <- armb |>
  select(-theta, -mle) |>
  filter(!is.null(fit)) |>
  mutate(
    perf = pmap(list(fit, n), function(FIT, N) {
      apply(na.omit(FIT$pars), 2, function(x) {
        c(
          'sd' = sd(x, na.rm = TRUE)
        )
      })
    })
  ) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path, "bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(
    perf = map(
      fit,
      ~ apply(.x$pars, 2, function(x) {
        c(
          'sd' = sd(x)
        )
      })
    )
  ) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path, "blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~ as.numeric(.x$out[1, ]))) |>
  select(-fit)


res <- armb |>
  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(mse = map2_dbl(gt, perf, ~ mean((.x - .y)^2))) |>
  group_by(n, p, nreps, step0, method, gammas, s, lab) |>
  summarise(mse = mean(mse, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path, "armb.qs")) |>
  mutate(times = map_dbl(fit, ~ if_else(is.null(.x$time), NA, .x$time))) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, step0, method, lab) |>
  summarise(times = mean(times, na.rm = TRUE))


bt_times <- qs::qread(file = paste0(path, "bt.qs")) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path, "blb.qs")) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times)

resID <- res |>
  left_join(times) |>
  mutate(design = "Identity")

# endregion
# region: LOGISTIC | EQUICORRELATION #####

# Identify the experiment
lab_experiment <- "ber_eq"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path, "gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS) {
  apply(PARS, 2, function(x) {
    c('sd' = sd(x))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path, "armb.qs"))
armb <- armb |>
  select(-theta, -mle) |>
  mutate(null = map_lgl(fit, ~ is.null(.x))) |>
  filter(!null) |>
  mutate(
    perf = pmap(list(fit, n), function(FIT, N) {
      apply(na.omit(FIT$pars), 2, function(x) {
        c(
          'sd' = sd(x, na.rm = TRUE)
        )
      })
    })
  ) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path, "bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(
    perf = map(
      fit,
      ~ apply(.x$pars, 2, function(x) {
        c(
          'sd' = sd(x)
        )
      })
    )
  ) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path, "blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~ as.numeric(.x$out[1, ]))) |>
  select(-fit)


res <- armb |>
  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(mse = map2_dbl(gt, perf, ~ mean((.x - .y)^2))) |>
  group_by(n, p, nreps, step0, method, gammas, s, lab) |>
  summarise(mse = mean(mse, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path, "armb.qs")) |>
  mutate(null = map_lgl(fit, ~ is.null(.x))) |>
  filter(!null) |>
  mutate(times = map_dbl(fit, ~ if_else(is.null(.x$time), NA, .x$time))) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, step0, method, lab) |>
  summarise(times = mean(times, na.rm = TRUE))


bt_times <- qs::qread(file = paste0(path, "bt.qs")) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path, "blb.qs")) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times)
resEQ <- res |>
  left_join(times) |>
  mutate(design = "Equicorrelation")
# endregion
# region LOGISTIC | TOEPLITZ #####

# Identify the experiment
lab_experiment <- "ber_toe"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path, "gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS) {
  apply(PARS, 2, function(x) {
    c('sd' = sd(x))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path, "armb.qs"))
armb <- armb |>
  select(-theta, -mle) |>
  filter(!is.null(fit)) |>
  mutate(
    perf = pmap(list(fit, n), function(FIT, N) {
      apply(na.omit(FIT$pars), 2, function(x) {
        c(
          'sd' = sd(x, na.rm = TRUE)
        )
      })
    })
  ) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path, "bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(
    perf = map(
      fit,
      ~ apply(.x$pars, 2, function(x) {
        c(
          'sd' = sd(x)
        )
      })
    )
  ) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path, "blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~ as.numeric(.x$out[1, ]))) |>
  select(-fit)


res <- armb |>
  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(mse = map2_dbl(gt, perf, ~ mean((.x - .y)^2))) |>
  group_by(n, p, nreps, step0, method, gammas, s, lab) |>
  summarise(mse = mean(mse, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path, "armb.qs")) |>
  mutate(times = map_dbl(fit, ~ if_else(is.null(.x$time), NA, .x$time))) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, step0, method, lab) |>
  summarise(times = mean(times, na.rm = TRUE))


bt_times <- qs::qread(file = paste0(path, "bt.qs")) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path, "blb.qs")) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times)


resTOE <- res |>
  left_join(times) |>
  mutate(design = "Toeplitz")

# endregion
# region POISSON | INDEPENDENT ###########
# Identify the experiment
lab_experiment <- "poi_id"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path, "gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS) {
  apply(PARS, 2, function(x) {
    c('sd' = sd(x))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path, "armb.qs"))
armb <- armb |>
  select(-theta, -mle) |>
  filter(!is.null(fit)) |>
  mutate(
    perf = pmap(list(fit, n), function(FIT, N) {
      apply(na.omit(FIT$pars), 2, function(x) {
        c(
          'sd' = sd(x, na.rm = TRUE)
        )
      })
    })
  ) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path, "bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(
    perf = map(
      fit,
      ~ apply(.x$pars, 2, function(x) {
        c(
          'sd' = sd(x)
        )
      })
    )
  ) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path, "blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~ as.numeric(.x$out[1, ]))) |>
  select(-fit)


res <- armb |>
  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(mse = map2_dbl(gt, perf, ~ mean((.x - .y)^2))) |>
  group_by(n, p, nreps, step0, method, gammas, s, lab) |>
  summarise(mse = mean(mse, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path, "armb.qs")) |>
  mutate(times = map_dbl(fit, ~ if_else(is.null(.x$time), NA, .x$time))) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, step0, method, lab) |>
  summarise(times = mean(times, na.rm = TRUE))


bt_times <- qs::qread(file = paste0(path, "bt.qs")) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path, "blb.qs")) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times)

resIDp <- res |>
  left_join(times) |>
  mutate(design = "Identity")
# endregion

# region Paper Plots ########

logL2 <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  ggplot(aes(x = n, y = mse, group = interaction(method, design, p))) +
  facet_nested(design ~ p, scales = 'free') +
  coord_cartesian(ylim = c(NA, 1e-1)) +
  scale_y_log10(labels = scales::scientific) +
  geom_line(
    aes(col = method, linetype = method)
  ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = ""
  ) +
  scale_x_continuous(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

logL2
ggsave(logL2, file = "output/logistic_sims.pdf", width = 8, height = 6)

logTime <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  ggplot(aes(x = n, y = times)) +
  facet_nested(~p, scales = 'free') +
  geom_line(
    aes(linetype = design, col = method)
  ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x = "Sample size", y = "Time (s)", col = "", linetype = "") +
  scale_y_log10(labels = scales::scientific) +
  scale_x_continuous(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
logTime
ggsave(
  logTime,
  file = "output/logistic_sims_time.pdf",
  width = 8,
  height = 6
)

poiL2 <- resIDp |>
  filter(is.na(step0)) |>
  # replace_na(list("step0" = 1e-5)) |>
  # filter(step0 < 1e-2) |>
  mutate(pn = p / n) |>
  filter(
    # metric %in% c("_msqdist"),
    lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    # metric = factor(
    #   metric,
    #   levels = c("_msqdist", "_time (s)"),
    #   labels = c("L2 distance\nfrom oracle diagonal", "Time (s)")
    # ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  ggplot(aes(x = n, y = mse, group = interaction(method, design, p))) +
  facet_nested(~p, scales = 'free') +
  coord_cartesian(ylim = c(NA, 1e-1)) +
  scale_y_log10(labels = scales::scientific) +
  # geom_point(aes(col = method)) +
  geom_line(
    aes(col = method, linetype = method)
    # , col = "grey8" linetype = factor(step0),
  ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  # scale_color_viridis_d() +
  # scale_fill_grey()+
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = ""
  ) +
  scale_x_continuous(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
poiTime <- resIDp |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g8s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resIDp$p),
      labels = paste0("p=", unique(resIDp$p))
    )
  ) |>
  ggplot(aes(x = n, y = times)) +
  facet_nested(~p, scales = 'free') +
  geom_line(aes(linetype = method, col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x = "Sample size", y = "Time (s)", col = "", linetype = "") +
  scale_y_log10(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggPoi <- ggpubr::ggarrange(
  poiL2 +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  poiTime,
  nrow = 2,
  common.legend = TRUE
)
ggsave(ggPoi, file = "output/poi_sims.pdf", width = 8, height = 6)
# endregion

# region Appendix Logistic ####

ggSupp1a <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    # metric %in% c("_msqdist"),
    lab %in%
      c(
        "armb_R100",
        "blb_R100g7s1",
        "blb_R100g7s2",
        "blb_R100g7s3",
        "bootstrap_R100"
      )
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    # metric = factor(
    #   metric,
    #   levels = c("_msqdist", "_time (s)"),
    #   labels = c("L2 distance\nfrom oracle diagonal", "Time (s)")
    # ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  replace_na(list(s = 1)) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(design ~ p, scales = 'free') +
  geom_line(aes(linetype = as.factor(s), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = "s:"
  ) +
  scale_y_log10(labels = scales::scientific) +
  coord_cartesian(ylim = c(NA, 5e-2)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggSupp1b <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    # metric %in% c("_msqdist"),
    lab %in%
      c(
        "armb_R100",
        "armb_R200",
        "blb_R100g7s1",
        "blb_R200g7s1",
        "bootstrap_R100"
      )
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    # metric = factor(
    #   metric,
    #   levels = c("_msqdist", "_time (s)"),
    #   labels = c("L2 distance\nfrom oracle diagonal", "Time (s)")
    # ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  replace_na(list(gammas = .7)) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(design ~ p, scales = 'free') +
  geom_line(aes(linetype = as.factor(nreps), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  coord_cartesian(ylim = c(NA, 5e-2)) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = unname(TeX("$R$:"))
  ) +
  scale_y_log10(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggSupp1 <- ggpubr::ggarrange(
  ggSupp1a +
    theme(legend.direction = "horizontal", legend.box = "vertical"),
  ggSupp1b +
    theme(legend.direction = "horizontal", legend.box = "vertical")
)
# ggSupp1
ggsave(ggSupp1, file = "output/logistic_sims_supp1.pdf", width = 8, height = 5)


ggSupp2a <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    # metric %in% c("_msqdist"),
    lab %in% c("armb_R100", "blb_R100g7s1", "blb_R100g8s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    # metric = factor(
    #   metric,
    #   levels = c("_msqdist", "_time (s)"),
    #   labels = c("L2 distance\nfrom oracle diagonal", "Time (s)")
    # ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  replace_na(list(gammas = .7)) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(design ~ p, scales = 'free') +
  geom_line(aes(linetype = as.factor(gammas), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = unname(TeX("BLB hyperparameter $\\gamma$:"))
  ) +
  scale_y_log10(labels = scales::scientific) +
  coord_cartesian(ylim = c(NA, 5e-2)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggSupp2b <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    # metric %in% c("_time (s)"),
    lab %in% c("armb_R100", "blb_R100g7s1", "blb_R100g8s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    # metric = factor(
    #   metric,
    #   levels = c("_msqdist", "_time (s)"),
    #   labels = c("L2 distance\nfrom oracle diagonal", "Time (s)")
    # ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  replace_na(list(gammas = .7)) |>
  ggplot(aes(x = n, y = times)) +
  facet_nested(design ~ p, scales = 'free') +
  geom_line(aes(linetype = as.factor(gammas), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Time (s)",
    col = "",
    linetype = unname(TeX("BLB hyperparameter $\\gamma$:"))
  ) +
  scale_y_log10(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggSupp2 <- ggpubr::ggarrange(ggSupp2a, ggSupp2b, common.legend = TRUE)
# ggSupp2
ggsave(ggSupp2, file = "output/logistic_sims_supp2.pdf", width = 8, height = 5)

ggSupp3a <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(auto = (is.na(step0) & method == "armb")) |>
  filter(!auto) |>
  replace_na(list("step0" = min(resID$step0, na.rm = TRUE))) |>
  # filter(step0 < 1e-2) |>
  mutate(pn = p / n) |>
  filter(
    # metric %in% c("_msqdist"),
    lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    # metric = factor(
    #   metric,
    #   levels = c("_msqdist", "_time (s)"),
    #   labels = c("L2 distance\nfrom oracle diagonal", "Time (s)")
    # ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  ggplot(aes(x = n, y = mse, group = interaction(method, design, p, step0))) +
  facet_nested(design ~ p, scales = 'free') +
  coord_cartesian(ylim = c(NA, 1e-1)) +
  scale_y_log10(labels = scales::scientific) +
  # geom_point(aes(col = method)) +
  geom_line(
    aes(col = method, linetype = factor(step0))
    # , col = "grey8" linetype = factor(step0),
  ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  # scale_color_viridis_d() +
  # scale_fill_grey()+
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = unname(TeX("ARM-B, $\\eta_0$:"))
  ) +
  scale_x_continuous(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggSupp3a <- ggSupp3a +
  guides(
    color = guide_legend(
      override.aes = list(
        color = c(NA, "darkgrey", "lightgrey"), # suppress color A
        label = c(TRUE, FALSE, FALSE)
      )
    ),
    linetype = guide_legend()
  )
ggsave(ggSupp3a, file = "output/logistic_sims_supp3.pdf", width = 8, height = 5)


# endregion
# region Appendix Poisson
poiSupp1a <- resIDp |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in%
      c(
        "armb_R100",
        "blb_R100g7s1",
        "blb_R100g7s2",
        "blb_R100g7s3",
        "bootstrap_R100"
      )
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resIDp$p),
      labels = paste0("p=", unique(resIDp$p))
    )
  ) |>
  replace_na(list(s = 1)) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(~p, scales = 'free') +
  geom_line(aes(linetype = as.factor(s), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = "s:"
  ) +
  scale_y_log10(labels = scales::scientific) +
  coord_cartesian(ylim = c(NA, 5e-2)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

poiSupp1b <- resIDp |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in%
      c(
        "armb_R100",
        "armb_R200",
        "blb_R100g7s2",
        "blb_R200g7s2",
        "bootstrap_R100"
      )
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resIDp$p),
      labels = paste0("p=", unique(resIDp$p))
    )
  ) |>
  replace_na(list(gammas = .7)) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(~p, scales = 'free') +
  geom_line(aes(linetype = as.factor(nreps), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  coord_cartesian(ylim = c(NA, 5e-2)) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = unname(TeX("$R$:"))
  ) +
  scale_y_log10(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
poiSupp1 <- ggpubr::ggarrange(
  poiSupp1a +
    theme(legend.direction = "horizontal", legend.box = "vertical"),
  poiSupp1b +
    theme(legend.direction = "horizontal", legend.box = "vertical")
)
# poiSupp1
ggsave(poiSupp1, file = "output/poisson_sims_supp1.pdf", width = 8, height = 3)

poiSupp2a <- resIDp |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g7s2", "blb_R100g8s2", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  replace_na(list(gammas = .7)) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(~p, scales = 'free') +
  geom_line(aes(linetype = as.factor(gammas), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = unname(TeX("BLB hyperparameter $\\gamma$:"))
  ) +
  scale_y_log10(labels = scales::scientific) +
  coord_cartesian(ylim = c(NA, 5e-2)) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


poiSupp2b <- resIDp |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g7s2", "blb_R100g8s2", "bootstrap_R100")
  ) |>
  mutate(
    design = factor(
      design,
      levels = c("Identity", "Equicorrelation", "Toeplitz"),
      ordered = TRUE
    ),
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(resID$p),
      labels = paste0("p=", unique(resID$p))
    )
  ) |>
  replace_na(list(gammas = .7)) |>
  ggplot(aes(x = n, y = times)) +
  facet_nested(~p, scales = 'free') +
  geom_line(aes(linetype = as.factor(gammas), col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Time (s)",
    col = "",
    linetype = unname(TeX("BLB hyperparameter $\\gamma$:"))
  ) +
  scale_y_log10(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

poiSupp2 <- ggpubr::ggarrange(poiSupp2a, poiSupp2b, common.legend = TRUE)
# poiSupp2
ggsave(poiSupp2, file = "output/poisson_sims_supp2.pdf", width = 8, height = 3)

# endregion
# region Random Intercept

lab_experiment <- "bin_id"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path, "setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path, "gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS) {
  apply(PARS, 2, function(x) {
    c('sd' = sd(x))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path, "armb.qs"))


armb <- armb |>
  select(-theta, -mle) |>
  filter(!is.null(fit)) |>
  mutate(
    perf = pmap(list(fit, n), function(FIT, N) {
      apply(na.omit(FIT$pars), 2, function(x) {
        c(
          'sd' = sd(x, na.rm = TRUE)
        )
      })
    })
  ) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path, "bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(err = map_lgl(fit, ~ class(.x)[1] == "try-error")) |>
  filter(!err) |>
  mutate(
    perf = map(
      fit,
      ~ apply(.x$pars, 2, function(x) {
        c(
          'sd' = sd(x)
        )
      })
    )
  ) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path, "blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~ class(.x)[1] == "try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~ .x$out)) |>
  select(-fit)


res <- armb |>
  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(mse = map2_dbl(gt, perf, ~ mean((.x - .y)^2))) |>
  group_by(n, p, nreps, step0, method, gammas, s, lab) |>
  summarise(mse = mean(mse, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path, "armb.qs")) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab, step0) |>
  summarise(times = mean(times))


bt_times <- qs::qread(file = paste0(path, "bt.qs")) |>
  mutate(err = map_lgl(fit, ~ class(.x)[1] == "try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path, "blb.qs")) |>
  mutate(err = map_lgl(fit, ~ class(.x) == "try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~ .x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times)

gg <- res |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(res$p),
      labels = paste0("p=", unique(res$p))
    )
  ) |>
  ggplot(aes(x = n, y = mse)) +
  facet_nested(~p, scales = 'free') +
  geom_point(aes(col = method)) +
  geom_line(
    aes(linetype = method, col = method)
  ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(
    x = "Sample size",
    y = "Mean squared error from oracle diagonal",
    col = "",
    linetype = ""
  ) +
  scale_y_log10(
    labels = scales::scientific,
    breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2)
  ) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

gg_time <- times |>
  filter(is.na(step0)) |>
  mutate(pn = p / n) |>
  filter(
    lab %in% c("armb_R100", "blb_R100g8s1", "bootstrap_R100")
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("armb", "blb", "bootstrap"),
      labels = c("ARM-B", "BLB", "Bootstrap")
    ),
    p = factor(
      p,
      levels = unique(res$p),
      labels = paste0("p=", unique(res$p))
    )
  ) |>
  ggplot(aes(x = n, y = times)) +
  facet_nested(~p, scales = 'free') +
  geom_line(aes(linetype = method, col = method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x = "Sample size", y = "Time (s)", col = "", linetype = "") +
  scale_y_log10(labels = scales::scientific) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggrndint <- ggpubr::ggarrange(
  gg +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ),
  gg_time,
  nrow = 2,
  common.legend = TRUE,
  heights = c(3, 2)
)
ggrndint
ggsave(ggrndint, file = "output/rndint_sims.pdf", width = 8, height = 6)

# endregion
