# This evaluates ARM-B on the APS dataset
library(tidyverse) |> suppressPackageStartupMessages()
library(readr) |> suppressPackageStartupMessages()
library(armb) |> suppressPackageStartupMessages()
library(patchwork) |> suppressPackageStartupMessages()
library(ROSE) |> suppressPackageStartupMessages()

#import
raw_data <- read_csv(
  "data/aps/aps_failure_training_set.csv",
  na = "na",
  skip = 19,
  show_col_types = FALSE
) |>
  as_tibble() |>
  mutate(class = as.numeric(factor(class)) - 1)
# mean(is.na(raw_data))

# drop poor rows and cols
nabyrow <- apply(raw_data, MARGIN = 1, FUN = function(x) mean(is.na(x)))
drop_data <- raw_data[nabyrow < .4, ]
drop_data <- drop_data[, sapply(raw_data, function(x) mean(is.na(x))) < .25]
drop_data <- drop_data[, c(
  TRUE,
  !(sapply(drop_data[, -1], function(x) length(unique(x))) <= 2)
)]

# dim(drop_data)
# mean(is.na(drop_data))
# table(drop_data$class)

# simple imputation
imputed_data <- drop_data |>
  mutate(across(c(-class), .fns = ~ replace_na(., median(., na.rm = TRUE)))) |>
  mutate(across(
    c(-class),
    .fns = ~ {
      as.numeric(scale(log(.x + .01)))
    }
  ))

# resample the original dataset
data.rose <- ROSE::ROSE(class ~ ., data = imputed_data, seed = 123)$data
# table(data.rose$class)

# compute the mle
n <- nrow(data.rose)
p <- ncol(data.rose)
fam = binomial(link = "logit")
mle <- glm(formula = "class~.-1", data = data.rose, family = fam)


# compute ROSE-resampled Bootstrap
ncores <- 4
tictoc::tic()
bt_list <- pbapply::pblapply(
  1:100,
  function(x) {
    DROSE <- ROSE(class ~ ., data = imputed_data, seed = x)$data
    fit <- biglm::bigglm(
      formula = terms(class ~ . - 1, data = DROSE),
      data = DROSE,
      family = fam
    )

    return(coef(fit))
  },
  cl = ncores
)
tictoc::toc()

btsd <- apply(Reduce(rbind, bt_list), 2, function(x) {
  sd(x)
})

# compute ROSE-resampled ARM-B
{
  tictoc::tic()
  step0 <- 1
  arm_ctrl <- list(
    BURN = 2 * n,
    STEPSIZE0 = step0,
    TRIM = .2 * n,
    VERBOSE = FALSE,
    SEED = 123,
    CONV_CHECK = TRUE,
    CONV_WINDOW = 2,
    TOL = 1e-2
  )
  tune_ctrl <- list(
    LENGTH = 1 * n,
    BURN = 0,
    MAXA = 20,
    SCALE = .5,
    AUTO = TRUE,
    VERBOSE = FALSE
  )

  set.seed(321)
  resample_idx <- sample(1:n, n, replace = TRUE)
  tune <- armb::tune_armGLM3(
    Y = data.rose$class[resample_idx],
    X = as.matrix(data.rose[resample_idx, -1]),
    FAMILY = fam$family,
    LINK = fam$link,
    THETA0 = coef(mle), #
    LENGTH = n,
    BURN = 0,
    STEPSIZE0 = arm_ctrl$STEPSIZE0,
    SCALE = tune_ctrl$SCALE,
    MAXA = tune_ctrl$MAXA,
    AUTO_STOP = tune_ctrl$AUTO,
    VERBOSE = tune_ctrl$VERBOSE,
    SEED = arm_ctrl$SEED,
    CONV_CHECK = arm_ctrl$CONV_CHECK,
    CONV_WINDOW = arm_ctrl$CONV_WINDOW,
    TOL = arm_ctrl$TOL
  )
  arm_ctrl$STEPSIZE0 <- tune$stepsizes[which.min(tune$nlls)]
  arm_ctrl$STEPSIZE0
  chains <- pbapply::pblapply(
    1:100,
    function(id) {
      DROSE <- ROSE::ROSE(class ~ ., data = imputed_data, seed = id)$data
      set.seed(id)
      chain.id <- armb::armGLM3(
        Y = as.numeric(DROSE$class),
        X = as.matrix(DROSE[, -1]),
        FAMILY = fam$family,
        LINK = fam$link,
        THETA0 = coef(mle),
        LENGTH = n,
        BURN = arm_ctrl$BURN,
        STEPSIZE0 = arm_ctrl$STEPSIZE0,
        TRIM = arm_ctrl$TRIM,
        VERBOSE = arm_ctrl$VERBOSE,
        SEED = arm_ctrl$SEED,
        CONV_CHECK = arm_ctrl$CONV_CHECK,
        CONV_WINDOW = arm_ctrl$CONV_WINDOW,
        TOL = arm_ctrl$TOL
      )
      return(chain.id)
    },
    cl = ncores
  )
  tictoc::toc()
}
armbtsd <- apply(
  purrr::reduce(purrr::map(chains, ~ .$avdelta), rbind),
  2,
  function(x) {
    sd(x)
  }
)


# Plot results
res <- tibble(
  var_lab = names(coef(mle)),
  mle = coef(mle),
  sd_btstrp = btsd,
  sd_armb = armbtsd
)

gg1 <- res |>
  ggplot(aes(x = sd_btstrp, y = sd_armb, col = mle)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "grey80") +
  geom_point() +
  theme_bw() +
  scale_x_continuous(limits = c(0.002, .015), breaks = c(0.005, 0.01, 0.015)) +
  scale_y_continuous(limits = c(0.002, .015), breaks = c(0.005, 0.01, 0.015)) +
  labs(
    x = "Bootstrap",
    y = "ARM-B",
    col = "MLE",
    subtitle = "Estimated Standard Deviations"
  ) +
  theme(legend.position = "bottom", plot.subtitle = element_text(hjust = 0.5)) +
  scale_colour_gradientn(colours = c("black", "grey80"))


gg <- res |>
  select(var_lab, mle, sd_btstrp, sd_armb) |>
  pivot_longer(
    cols = starts_with("sd"),
    names_to = "method",
    values_to = "sd"
  ) |>
  mutate(
    method = factor(
      method,
      levels = c("sd_armb", "sd_btstrp"),
      labels = c("ARM-B", "Oracle")
    )
  ) |>
  ggplot(aes(x = var_lab, y = mle)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "lightgrey") +
  geom_linerange(
    aes(color = method, ymin = mle - 1.96 * sd, ymax = mle + 1.96 * sd),
    position = position_dodge(width = 1),
    linewidth = 1
  ) +
  scale_color_grey(start = .8, end = .2) +
  labs(
    y = "95% Confidence Intervals",
    x = "Regression coefficients",
    col = "",
    subtitle = ""
  )
gg
ggsave(gg, filename = "output/aps.pdf", width = 10, height = 6)
