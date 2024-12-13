library(tidyverse) |> suppressPackageStartupMessages()
library(readr) |> suppressPackageStartupMessages()
library(armb) |> suppressPackageStartupMessages()
library(patchwork) |> suppressPackageStartupMessages()
library(ROSE) |> suppressPackageStartupMessages()

#import
raw_data <- read_csv("data/aps_failure_training_set.csv",
                     na = "na", skip = 19, show_col_types = FALSE) |>  as_tibble()  |>
  mutate(class = as.numeric(factor(class))-1)
# mean(is.na(raw_data))

# drop poor rows and cols
nabyrow <- apply(raw_data, MARGIN = 1,FUN =  function(x) mean(is.na(x)))
drop_data <- raw_data[nabyrow<.3,]
drop_data <- drop_data[,sapply(raw_data, function(x) mean(is.na(x)))<.7]
drop_data <- drop_data[,c(TRUE,!(sapply(drop_data[,-1], function(x) length(unique(x)))<=2))]

# dim(drop_data)
# mean(is.na(drop_data))
# table(drop_data$class)

# simple imputation
imputed_data <- drop_data |>
  mutate(across(c(-class), .fns = ~replace_na(.,median(., na.rm=TRUE)))) |>
  mutate(across(c(-class), .fns = ~ {
    as.numeric(scale(log(.x+1)))
  }))

# resample the original dataset
data.rose <- ROSE(class ~ ., data = imputed_data, seed = 123)$data
# table(data.rose$class)

# compute the mle
n <- nrow(data.rose)
p <- ncol(data.rose)
fam=binomial(link="logit")
mle <- glm(formula = "class~.-1", data = data.rose, family = fam)


# compute ROSE-resampled Bootstrap
ncores <- 6
tictoc::tic()
bt_list <- pbapply::pblapply(1:500, function(x){
  DROSE <- ROSE(class ~ ., data = imputed_data, seed = x)$data
  fit <- biglm::bigglm(formula = terms(class~.-1, data = DROSE), data = DROSE,
                       family = fam)

  return(coef(fit))
}, cl = ncores)
tictoc::toc()

btsd <- apply(Reduce(rbind, bt_list), 2, function(x) {sd(x)})

# compute ROSE-resampled ARM-B
tictoc::tic()
{
  brn <- 1*n
  arm_ctrl <- list(
    MAXT           = brn+n,
    BURN           = brn,
    BATCH          = 1,
    STEPSIZE0      = 1e-3,
    PAR1           = 1,
    PAR2           = 1e-4,
    PAR3           = .75,
    PATH_WINDOW    = n,
    VERBOSE_WINDOW = 1,
    VERBOSE        = F,
    SEED           = 123
  )

  tune_ctrl <- list(
    MAXT    = arm_ctrl$MAXT,
    MAXA    = 200,
    SCALE   = .3,
    AUTO    = TRUE,
    VERBOSE = FALSE,
    CONV_WINDOW = 100,
    CONV_TOL = 1e-5
  )

  tune <- tune_armGLM(
    Y = data.rose$class,
    X = as.matrix(data.rose[,-1]),
    FAMILY = fam$family,
    LINK = fam$link,
    THETA0 = coef(mle), #
    MAXT = arm_ctrl$MAXT,
    BURN = arm_ctrl$BURN,
    BATCH = arm_ctrl$BATCH,
    STEPSIZE0 = arm_ctrl$STEPSIZE0,
    SCALE = tune_ctrl$SCALE,
    MAXA = tune_ctrl$MAXA,
    PAR1 = arm_ctrl$PAR1,
    PAR2 = arm_ctrl$PAR2,
    PAR3 = arm_ctrl$PAR3,
    AUTO_STOP = tune_ctrl$AUTO,
    VERBOSE = tune_ctrl$VERBOSE,
    SEED = arm_ctrl$SEED,
    SKIP_PRINT = 0L
  )
  arm_ctrl$STEPSIZE0 <- tune$stepsizes[which.min(tune$devresids)]
  armb_list <- pbapply::pblapply(1:100, function(x){
    DROSE <- ROSE(class ~ ., data = imputed_data, seed = x)$data
    fit <- armGLM(Y = as.numeric(DROSE$class),
                  X = as.matrix(DROSE[,-1]),
                  FAMILY = fam$family,
                  LINK = fam$link,
                  THETA0 = coef(mle),
                  MAXT = arm_ctrl$MAXT,
                  BURN = arm_ctrl$BURN,
                  BATCH = arm_ctrl$BATCH,
                  STEPSIZE0 = arm_ctrl$STEPSIZE0,
                  PAR1 = arm_ctrl$PAR1,
                  PAR2 = arm_ctrl$PAR2,
                  PAR3 = arm_ctrl$PAR3,
                  PATH_WINDOW = arm_ctrl$PATH_WINDOW,
                  VERBOSE_WINDOW = arm_ctrl$VERBOSE_WINDOW,
                  VERBOSE = arm_ctrl$VERBOSE,
                  SEED = x
    )
    return(fit$avtheta)
  }, cl = ncores)
}
tictoc::toc()
armbtsd <- apply(Reduce(rbind, armb_list), 2, function(x) {sd(x)})



# Plot results
res <- tibble(var_lab = names(coef(mle)), mle = coef(mle), sd_btstrp = btsd, sd_armb = armbtsd)
qs::qsave(res, "output/aps.qs")

# average ci length ratio
res |>
  mutate(ci_length_ratio = sd_btstrp/sd_armb) |>
  summarise_at("ci_length_ratio", mean)
gg1 <- res |> ggplot(aes(x=sd_btstrp, y = sd_armb, col = mle))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", col ="grey80")+
  geom_point() +
  theme_bw()+
  # scale_x_continuous(limits = c(0.002, .015), breaks = c(0.005, 0.01, 0.015))+
  # scale_y_continuous(limits = c(0.002, .015), breaks = c(0.005, 0.01, 0.015))+
  labs(x = "Oracle", y= "ARM-B", col = "MLE", subtitle = "Estimated Standard Deviations")+
  theme(legend.position = "bottom", plot.subtitle = element_text(hjust = 0.5))+
  scale_colour_gradientn(colours = c("black", "grey80"), breaks=c(-.1,0, .15), labels=c(-.1,0, .15))


gg2 <- res |>
  select(var_lab, mle, sd_btstrp, sd_armb) |>
  pivot_longer(cols=starts_with("sd"), names_to = "method", values_to = "sd") |>
  mutate(method=factor(method, levels=c("sd_armb", "sd_btstrp"), labels=c("ARM-B", "Oracle"))) |>
  ggplot(aes(x=mle, y = var_lab)) +
  geom_vline(xintercept = 0, linetype="dashed", col="lightgrey")+
  geom_segment(aes(x = mle-1.96*sd, xend = mle+1.96*sd, col = method),
               position = position_dodge(width = 0.1))+
  theme_bw()+
  theme(
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank(),  #remove y axis ticks
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",plot.subtitle = element_text(hjust = 0.5)
  )+
  scale_color_grey(start = .8, end =.2 )+
  labs(x="MLE", y="Regression coefficients", col="", subtitle = "95% Confidence Intervals")
gg2

ggsave(gg, filename="output/aps_failure.png", width = 10, height = 6)


gg2 <- res |>
  select(var_lab, mle, sd_btstrp, sd_armb) |>
  pivot_longer(cols=starts_with("sd"), names_to = "method", values_to = "sd") |>
  mutate(method=factor(method, levels=c("sd_armb", "sd_btstrp"), labels=c("ARM-B", "Oracle"))) |>
  ggplot(aes(x=var_lab, y = mle))+
  theme_bw()+
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",plot.subtitle = element_text(hjust = 0.5)
  )  +
  geom_hline(yintercept = 0, linetype="dashed", col="lightgrey") +
  geom_linerange(aes(color = method,
                     ymin = mle-1.96*sd,
                     ymax = mle+1.96*sd),
                 position = position_dodge(width = 1), size = 1)+
  scale_color_grey(start = .8, end =.2 )+
  labs(y="MLE", x="Regression coefficients", col="", subtitle = "95% Confidence Intervals")
gg <- gg1+gg2+plot_layout(widths = c(1,2))
# gg
ggsave(gg, filename="output/aps_failure.pdf", width = 10, height = 6)
ggsave(gg2, filename="output/aps_failure_ci.pdf", width = 10, height = 6)
# gg2



