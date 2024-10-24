library(tidyverse)

# Identify the experiment
lab_experiment <- "poi_toe_rnd"
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
armb |>
  mutate(
    gamma = map_dbl(fit, ~.x$gamma),
    step = map_dbl(fit, ~.x$step0)
  ) |> select(n, p, gamma, step)

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

gg <- res |>
  bind_rows(times) |>
  ggplot(aes(x = n, y = val, col = method)) +
  facet_grid(metric~p, scales = 'free')+
  # geom_point() +
  geom_line(aes(linetype = as.factor(lab))) +
  theme_bw() +
  # scale_color_grey() +
  # scale_fill_grey()+
  scale_y_log10() +
  theme(legend.position = 'top') +
  labs(x = 'Number of observations', y = ' ', col = '', col = '')
plotly::ggplotly(gg, dynamickTicks = TRUE)


##########
# armb <- qs::qread(file = paste0(path,"armb.qs"))
# set.seed(1)
# tmpP <- 50
# S <- matrix(.2, tmpP, tmpP)
# diag(S) <- 1
# X <- mvtnorm::rmvnorm(5000, mean = rep(0,tmpP), sigma = S)
# MU <- resp$linkinv(X%*%rep(1,tmpP)/sqrt(tmpP))
# Y <- sapply(MU, function(mui) rpois(1, lambda = mui))
# summary(Y)
#
# runif(tmpP, -1, 1)
#
# armb <- qs::qread(file = paste0(path,"armb.qs"))
#
# dttmp <- genData_experiment(rep(1,500), N=5000, P=500, SEED=1)
# summary(dttmp$Y)
# theta_tmp <- armb$mle[[1]]#rep(1, 50)
# tst <- armb::test_glm(as.numeric(dttmp$Y), dttmp$X, resp$family, resp$link, theta_tmp)
# resp$linkinv(dttmp$X%*%theta_tmp)
# testthat::expect_equal(tst$mu, as.numeric(resp$linkinv(dttmp$X%*%theta_tmp)))
