library(tidyverse)|> suppressPackageStartupMessages()
library(latex2exp)|> suppressPackageStartupMessages()
library(ggh4x)|> suppressPackageStartupMessages()
library(ggpubr)|> suppressPackageStartupMessages()
###### LOGISTIC | INDEPENDENT ###########
# Identify the experiment
lab_experiment <- "ber_id"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path,"gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS){
  apply(PARS, 2, function(x) {
    c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path,"armb.qs"))
# armb |>
#   mutate(
#     gamma = map_dbl(fit, ~.x$gamma),
#     step = map_dbl(fit, ~.x$step0)
#   ) |> select(n, p, gamma, step)

armb <- armb |>
  select(-theta, -mle) |>
  mutate(perf = pmap(list(fit, n), function(FIT, N) {
    apply(na.omit(FIT$pars), 2, function(x) {
      c('sd' = sqrt(var(x, na.rm = TRUE)*(N^FIT$gamma*FIT$batch/N)), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975), na.rm = TRUE))))
    })})) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path,"bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(perf = map(fit, ~apply(.x$pars, 2,
                                function(x) {
                                  c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
                                }))) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path,"blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~.x$out)) |>
  select(-fit)











res <- armb |> mutate(gammas=NA) |>

  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(
    metric_msqdist = map2_dbl(perf, gt, ~mean((.x[1,]-.y[1,])^2, na.rm=TRUE)),
    metric_mabsdist = map2_dbl(perf, gt, ~mean(abs(.x[1,]-.y[1,]), na.rm=TRUE)),
    metric_ir = map2_dbl(perf, gt, ~mean(abs(.x[2,]-.y[2,]), na.rm=TRUE))) |>
  select(-perf, -gt) |>
  pivot_longer(cols = starts_with('metric'),
               names_to = 'metric', values_to = 'val') |>
  mutate(metric = stringi::stri_extract(metric, regex = "_(.*)$")) |>
  group_by(n, p, nreps, method, gammas, s, lab, metric) |>
  summarise(val = mean(val, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path,"armb.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))



bt_times <- qs::qread(file = paste0(path,"bt.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path,"blb.qs")) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times) |>
  mutate(metric = '_time (s)') |>
  rename('val' = times)

resID <- res |>
  bind_rows(times) |>
  mutate(design="Identity")

###### LOGISTIC | EQUICORRELATION #####

# Identify the experiment
lab_experiment <- "ber_eq"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path,"gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS){
  apply(PARS, 2, function(x) {
    c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path,"armb.qs"))
# armb |>
#   mutate(
#     gamma = map_dbl(fit, ~.x$gamma),
#     step = map_dbl(fit, ~.x$step0)
#   ) |> select(n, p, gamma, step)

armb <- armb |>
  select(-theta, -mle) |>
  mutate(perf = pmap(list(fit, n), function(FIT, N) {
    apply(na.omit(FIT$pars), 2, function(x) {
      c('sd' = sqrt(var(x, na.rm = TRUE)*(N^FIT$gamma*FIT$batch/N)), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975), na.rm = TRUE))))
    })})) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path,"bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(perf = map(fit, ~apply(.x$pars, 2,
                                function(x) {
                                  c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
                                }))) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path,"blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~.x$out)) |>
  select(-fit)











res <- armb |> mutate(gammas=NA) |>

  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(
    metric_msqdist = map2_dbl(perf, gt, ~mean((.x[1,]-.y[1,])^2, na.rm=TRUE)),
    metric_mabsdist = map2_dbl(perf, gt, ~mean(abs(.x[1,]-.y[1,]), na.rm=TRUE)),
    metric_ir = map2_dbl(perf, gt, ~mean(abs(.x[2,]-.y[2,]), na.rm=TRUE))) |>
  select(-perf, -gt) |>
  pivot_longer(cols = starts_with('metric'),
               names_to = 'metric', values_to = 'val') |>
  mutate(metric = stringi::stri_extract(metric, regex = "_(.*)$")) |>
  group_by(n, p, nreps, method, gammas, s, lab, metric) |>
  summarise(val = mean(val, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path,"armb.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))



bt_times <- qs::qread(file = paste0(path,"bt.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path,"blb.qs")) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times) |>
  mutate(metric = '_time (s)') |>
  rename('val' = times)

resEQ <- res |>
  bind_rows(times) |>
  mutate(design="Equicorrelation")

###### LOGISTIC | TOEPLITZ #####

# Identify the experiment
lab_experiment <- "ber_toe"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path,"gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS){
  apply(PARS, 2, function(x) {
    c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path,"armb.qs"))
# armb |>
#   mutate(
#     gamma = map_dbl(fit, ~.x$gamma),
#     step = map_dbl(fit, ~.x$step0)
#   ) |> select(n, p, gamma, step)

armb <- armb |>
  select(-theta, -mle) |>
  mutate(perf = pmap(list(fit, n), function(FIT, N) {
    apply(na.omit(FIT$pars), 2, function(x) {
      c('sd' = sqrt(var(x, na.rm = TRUE)*(N^FIT$gamma*FIT$batch/N)), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975), na.rm = TRUE))))
    })})) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path,"bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(perf = map(fit, ~apply(.x$pars, 2,
                                function(x) {
                                  c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
                                }))) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path,"blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~.x$out)) |>
  select(-fit)











res <- armb |> mutate(gammas=NA) |>

  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(
    metric_msqdist = map2_dbl(perf, gt, ~mean((.x[1,]-.y[1,])^2, na.rm=TRUE)),
    metric_mabsdist = map2_dbl(perf, gt, ~mean(abs(.x[1,]-.y[1,]), na.rm=TRUE)),
    metric_ir = map2_dbl(perf, gt, ~mean(abs(.x[2,]-.y[2,]), na.rm=TRUE))) |>
  select(-perf, -gt) |>
  pivot_longer(cols = starts_with('metric'),
               names_to = 'metric', values_to = 'val') |>
  mutate(metric = stringi::stri_extract(metric, regex = "_(.*)$")) |>
  group_by(n, p, nreps, method, gammas, s, lab, metric) |>
  summarise(val = mean(val, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path,"armb.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))



bt_times <- qs::qread(file = paste0(path,"bt.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path,"blb.qs")) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times) |>
  mutate(metric = '_time (s)') |>
  rename('val' = times)

resTOE <- res |>
  bind_rows(times) |>
  mutate(design="Toeplitz")


###### POISSON | INDEPENDENT ###########
# Identify the experiment
lab_experiment <- "poi_id"
path <- paste0('output/', lab_experiment, '/')
load(file = paste0(path,"setup.rda"))


# Read ground truth simulations and compute quantities of interest
gt <- qs::qread(file = paste0(path,"gt.qs"))
gt$gt <- lapply(gt$gt, function(PARS){
  apply(PARS, 2, function(x) {
    c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
  })
})

# Read armb simulations and compute quantities of interest
armb <- qs::qread(file = paste0(path,"armb.qs"))
# armb |>
#   mutate(
#     gamma = map_dbl(fit, ~.x$gamma),
#     step = map_dbl(fit, ~.x$step0)
#   ) |> select(n, p, gamma, step)

armb <- armb |>
  select(-theta, -mle) |>
  mutate(perf = pmap(list(fit, n), function(FIT, N) {
    apply(na.omit(FIT$pars), 2, function(x) {
      c('sd' = sqrt(var(x, na.rm = TRUE)*(N^FIT$gamma*FIT$batch/N)), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975), na.rm = TRUE))))
    })})) |>
  select(-fit)

# Read bootstrap simulations and compute quantities of interest
bt <- qs::qread(file = paste0(path,"bt.qs"))
bt <- bt |>
  select(-theta, -mle) |>
  mutate(perf = map(fit, ~apply(.x$pars, 2,
                                function(x) {
                                  c('sd' = sd(x), 'ir' = as.numeric(diff(quantile(x, probs = c(.025, .975)))))
                                }))) |>
  select(-fit)

# Read blb simulations and compute quantities of interest
blb <- qs::qread(file = paste0(path,"blb.qs"))
blb <- blb |>
  select(-theta) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(perf = map(fit, ~.x$out)) |>
  select(-fit)











res <- armb |> mutate(gammas=NA) |>

  bind_rows(bt) |>
  bind_rows(blb) |>
  left_join(gt, by = c('n', 'p')) |>
  mutate(
    metric_msqdist = map2_dbl(perf, gt, ~mean((.x[1,]-.y[1,])^2, na.rm=TRUE)),
    metric_mabsdist = map2_dbl(perf, gt, ~mean(abs(.x[1,]-.y[1,]), na.rm=TRUE)),
    metric_ir = map2_dbl(perf, gt, ~mean(abs(.x[2,]-.y[2,]), na.rm=TRUE))) |>
  select(-perf, -gt) |>
  pivot_longer(cols = starts_with('metric'),
               names_to = 'metric', values_to = 'val') |>
  mutate(metric = stringi::stri_extract(metric, regex = "_(.*)$")) |>
  group_by(n, p, nreps, method, gammas, s, lab, metric) |>
  summarise(val = mean(val, na.rm = TRUE))

armb_times <- qs::qread(file = paste0(path,"armb.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))



bt_times <- qs::qread(file = paste0(path,"bt.qs")) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -mle, -fit) |>
  group_by(n, p, nreps, method, lab) |>
  summarise(times = mean(times))

blb_times <- qs::qread(file = paste0(path,"blb.qs")) |>
  mutate(err = map_lgl(fit, ~class(.x)=="try-error")) |>
  filter(!err) |>
  mutate(times = map_dbl(fit, ~.x$time)) |>
  select(-theta, -fit) |>
  group_by(n, p, nreps, method, gammas, s, lab) |>
  summarise(times = mean(times))

times <- armb_times |>
  bind_rows(bt_times) |>
  bind_rows(blb_times) |>
  mutate(metric = '_time (s)') |>
  rename('val' = times)

resIDp <- res |>
  bind_rows(times) |>
  mutate(design="Identity")


#### Paper Plots ########






logL2 <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_msqdist"),
          lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(design~p, scales = 'free')+
  # geom_point(aes(col=method)) +
  geom_line(aes(linetype = method, col=method)
            # , col = "grey8"
            ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  # scale_fill_grey()+
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype="")+
  scale_y_log10(labels = scales::scientific)+
  coord_cartesian(ylim=c(NA, 1e-1))+
  # scale_x_continuous(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))
# logL2
ggsave(logL2, file="output/logistic_sims.pdf", width = 8, height = 6)

logTime <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_time (s)"),
          lab %in% c("armb_R100", "blb_R100g7s1", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  # geom_point(aes(col=method)) +
  geom_line(aes(linetype = design, col=method)
            # , col = "grey8"
  ) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  # scale_fill_grey()+
  labs(x="Sample size", y="Time (s)", col="", linetype="")+
  scale_y_log10(labels = scales::scientific)+
  # coord_cartesian(ylim=c(NA, 1e-1))+
  # scale_x_continuous(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))
# logTime
ggsave(logTime, file="output/logistic_sims_time.pdf", width = 8, height = 6)

poiL2 <- resIDp |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_msqdist"),
          lab %in% c("armb_R100", "blb_R100g8s3", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resIDp$p), labels=paste0("p=", unique(resIDp$p)))
  ) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  geom_line(aes(linetype = method, col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype="")+
  scale_y_log10(labels = scales::scientific)+
  coord_cartesian(ylim=c(NA, 1e-3))+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

poiTime <- resIDp |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_time (s)"),
          lab %in% c("armb_R100", "blb_R100g8s1", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resIDp$p), labels=paste0("p=", unique(resIDp$p)))
  ) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  geom_line(aes(linetype = method, col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="Time (s)", col="", linetype="")+
  scale_y_log10(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

ggPoi <- ggarrange(poiL2+theme(
  axis.text.x = element_blank(),
  axis.title.x = element_blank()
  ), poiTime, nrow = 2, common.legend = TRUE)
ggsave(ggPoi, file="output/poi_sims.pdf", width = 8, height = 6)

#### Appendix Logistic ####

ggSupp1a <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_msqdist"),
          lab %in% c("armb_R100", "blb_R100g7s1", "blb_R100g7s2", "blb_R100g7s3", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(s=1)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(design~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(s), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype="s:")+
  scale_y_log10(labels = scales::scientific)+
  coord_cartesian(ylim=c(NA, 5e-2))+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))
ggSupp1b <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(pn = p/n) |>
  filter(metric%in%c("_msqdist"),
         lab %in% c("armb_R100", "armb_R200", "blb_R100g7s1", "blb_R200g7s1", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(gammas=.7)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(design~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(nreps), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  coord_cartesian(ylim=c(NA, 5e-2))+
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype=unname(TeX("$R$:")))+
  scale_y_log10(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))



ggSupp1 <- ggpubr::ggarrange(
  ggSupp1a+
    theme(legend.direction = "horizontal", legend.box = "vertical"),
  ggSupp1b+
    theme(legend.direction = "horizontal", legend.box = "vertical")
)
# ggSupp1
ggsave(ggSupp1, file="output/logistic_sims_supp1.pdf", width = 8, height = 5)


ggSupp2a <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_msqdist"),
          lab %in% c("armb_R100", "blb_R100g7s1", "blb_R100g8s1", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(gammas=.7)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(design~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(gammas), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype=unname(TeX("BLB hyperparameter $\\gamma$:")))+
  scale_y_log10(labels = scales::scientific)+
  coord_cartesian(ylim=c(NA, 5e-2))+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))


ggSupp2b <- resID |>
  bind_rows(resEQ) |>
  bind_rows(resTOE) |>
  mutate(pn = p/n) |>
  filter(metric%in%c("_time (s)"),
          lab %in% c("armb_R100", "blb_R100g7s1", "blb_R100g8s1", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(gammas=.7)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(design~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(gammas), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="Time (s)", col="", linetype=unname(TeX("BLB hyperparameter $\\gamma$:")))+
  scale_y_log10(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

ggSupp2 <- ggpubr::ggarrange(ggSupp2a, ggSupp2b, common.legend = TRUE)
# ggSupp2
ggsave(ggSupp2, file="output/logistic_sims_supp2.pdf", width = 8, height = 5)






#### Appendix Poisson ####
poiSupp1a <- resIDp |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_msqdist"),
          lab %in% c("armb_R100", "blb_R100g7s1", "blb_R100g7s2", "blb_R100g7s3", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(s=1)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(s), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype="s:")+
  scale_y_log10(labels = scales::scientific)+
  coord_cartesian(ylim=c(NA, 5e-2))+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

poiSupp1b <- resIDp |>
  mutate(pn = p/n) |>
  filter(metric%in%c("_msqdist"),
         lab %in% c("armb_R100", "armb_R200", "blb_R100g7s2", "blb_R200g7s2", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(gammas=.7)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(nreps), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  coord_cartesian(ylim=c(NA, 5e-2))+
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype=unname(TeX("$R$:")))+
  scale_y_log10(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))
poiSupp1 <- ggpubr::ggarrange(
  poiSupp1a+
    theme(legend.direction = "horizontal", legend.box = "vertical"),
  poiSupp1b+
    theme(legend.direction = "horizontal", legend.box = "vertical")
)
# poiSupp1
ggsave(poiSupp1, file="output/poisson_sims_supp1.pdf", width = 8, height = 3)



poiSupp2a <- resIDp |>
  mutate(pn = p/n) |>
  filter( metric%in%c("_msqdist"),
          lab %in% c("armb_R100", "blb_R100g7s2", "blb_R100g8s2", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(gammas=.7)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(gammas), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="L2 distance from oracle diagonal", col="", linetype=unname(TeX("BLB hyperparameter $\\gamma$:")))+
  scale_y_log10(labels = scales::scientific)+
  coord_cartesian(ylim=c(NA, 5e-2))+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "none", axis.text.x = element_text(angle = 45, hjust=1))


poiSupp2b <- resIDp |>
  mutate(pn = p/n) |>
  filter(metric%in%c("_time (s)"),
         lab %in% c("armb_R100", "blb_R100g7s2", "blb_R100g8s2", "bootstrap_R100")) |>
  mutate(
    design = factor(design, levels=c("Identity", "Equicorrelation", "Toeplitz"), ordered = TRUE),
    metric = factor(metric, levels=c("_msqdist", "_time (s)"), labels=c("L2 distance\nfrom oracle diagonal", "Time (s)")),
    method = factor(method, levels=c("armb","blb", "bootstrap"), labels=c("ARM-B", "BLB", "Bootstrap")),
    p=factor(p, levels=unique(resID$p), labels=paste0("p=", unique(resID$p)))
  ) |>
  replace_na(list(gammas=.7)) |>
  ggplot(aes(x = n, y = val)) +
  facet_nested(~p, scales = 'free')+
  geom_line(aes(linetype = as.factor(gammas), col=method)) +
  theme_bw() +
  scale_color_grey(start = .1, end = .8) +
  labs(x="Sample size", y="Time (s)", col="", linetype=unname(TeX("BLB hyperparameter $\\gamma$:")))+
  scale_y_log10(labels = scales::scientific)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

poiSupp2 <- ggpubr::ggarrange(poiSupp2a, poiSupp2b, common.legend = TRUE)
# poiSupp2
ggsave(poiSupp2, file="output/poisson_sims_supp2.pdf", width = 8, height = 3)
