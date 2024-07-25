
#########################################################
## title: predicting redlist or pop change status

# -------------------------------------------------------
# package loading
library(tidyverse)
library(survey)
library(mgcv)
library(dplyr)

#########################################################
## load in data
load("C:/Data1/UK/tolerance_output_with_conservation.RData")

d <- tolerance_output_with_conservation |>
  mutate(redamber = ifelse(uk_conservation_status < 3, 1, ifelse(uk_conservation_status == 4, NA, 0))) |>
  mutate(red = ifelse(uk_conservation_status == 1, 1, ifelse(uk_conservation_status < 4, 0, NA)))


cor(d$mean_peak, d$mean_breadth_0.5)
cor(d$ci_range_peak, d$ci_range_breadth_0.5)

#########################################################
## try the models


mod1 <- glm(redamber ~ mean_breadth_0.5*mean_peak, family = binomial, data = d)
mod2 <- glm(redamber ~ mean_breadth_0.5 + mean_peak, family = binomial, data = d)
mod3 <- glm(redamber ~ mean_peak, family = binomial, data = d)
mod4 <- glm(redamber ~ mean_breadth_0.5, family = binomial, data = d)

mod1r <- glm(red ~ mean_breadth_0.5*mean_peak, family = binomial, data = d)
mod2r <- glm(red ~ mean_breadth_0.5 + mean_peak, family = binomial, data = d)
mod3r <- glm(red ~ mean_peak, family = binomial, data = d)
mod4r <- glm(red ~ mean_breadth_0.5, family = binomial, data = d)

#########################################################
## with weights: (need to use a different function)

# need to use function svydesign
# I tried to do this, but it is a bit confusing and 
# I couldn't figure it out in the time I have. :( 


#########################################################
## try some gams with no weights: 
library(mgcv)

gam0 <- gam(redamber ~ 1, family = binomial, data = d)
gam1 <- gam(redamber ~ s(mean_breadth_0.5, mean_peak), family = binomial, data = d, gamma = 1.4)
gam2 <- gam(redamber ~ s(mean_breadth_0.5) + s(mean_peak), family = binomial, data = d, gamma = 1.4)

par(mfrow = c(1,2))
vis.gam(gam1, too.far = 0.1, plot.type = "contour", type = "response", xlim = c(0, 50))
vis.gam(gam2, too.far = 0.1, plot.type = "contour", type = "response", xlim = c(0, 50))


#########################################################
## population trend. 

# transform so that a doubling is equivalent to a halving: 

d$log_population_change_10 <- log10(d$population_change_10/100 + 1)
d$log_population_change_27 <- log10(d$population_change_27/100 + 1)
d$wt <- 1/d$ci_range_breadth_0.5

# i don't think we should do it with the 1-year and 5-year trends. 

mod5 <- glm(log_population_change_10 ~ mean_breadth_0.5*mean_peak, data = d, weights = wt)
mod6 <- glm(log_population_change_10 ~ mean_breadth_0.5 + mean_peak, data = d, weights = wt)

mod7 <- glm(log_population_change_27 ~ mean_breadth_0.5*mean_peak, data = d, weights = wt)
mod8 <- glm(log_population_change_27 ~ mean_breadth_0.5 + mean_peak, data = d, weights = wt)

# larger percentage of deviance explained with models 5 and 6. 


#########################################################
## plot predictions

pred <- predict(mod5, d, type = "response")
plot(d$log_population_change_10, pred)

# not much predictive power here! 
# maybe due to no underlying association or small dataset and low power?


#########################################################
## gams for population change

gam5 <- gam(log_population_change_10 ~ s(mean_breadth_0.5, mean_peak), data = d, weights = wt, gamma = 1.4)
vis.gam(gam5, type = "response", plot.type = "contour", too.far = 0.1)

#########################################################

## plotting out mod 7.2

mod7.2 <- glm(log_population_change_27 ~ mean_breadth_0.5*mean_peak + mean_right_proportion_breadth_0.5, data = d, weights = wt)

nd1 <- data.frame(mean_breadth_0.5 = mean(d$mean_breadth_0.5),
                  mean_peak = mean(d$mean_peak),
                  mean_right_proportion_breadth_0.5 = seq(0, 1, by = 0.01))

pred <- predict(mod7.2, nd1, type = "link", se.fit = TRUE)

est <- pred$fit

lcl <- est - 1.96*pred$se.fit

ucl <- est + 1.96*pred$se.fit

plot(nd1$mean_right_proportion_breadth_0.5, est, type = "l", lwd = 2,
     xlab = "Proportion right of peak",
     ylab = "Estimated population trend",
     ylim = c(-3, 3))

polygon(x = c(nd1$mean_right_proportion_breadth_0.5, rev(nd1$mean_right_proportion_breadth_0.5), nd1$mean_right_proportion_breadth_0.5[1]),
        y = c(lcl, rev(ucl), lcl[1]),
        col = alpha("firebrick", 0.3), border = alpha("white", 0))

lines(nd1$mean_right_proportion_breadth_0.5, est, lwd = 2)

nd2 <- expand.grid(mean_breadth_0.5 = seq(min(d$mean_breadth_0.5), max(d$mean_breadth_0.5), length.out = 20),
                   mean_peak = seq(0, 50, length.out = 50),
                   mean_right_proportion_breadth_0.5 = 0.5)

pred <- predict(mod7.2, nd2, type = "link", se.fit = TRUE)

nd2$pred <- pred$fit

est <- pred$fit
lcl <- est - 1.96*pred$se.fit
ucl <- est + 1.96*pred$se.fit

nd2 <- arrange(nd2, mean_breadth_0.5, mean_peak)

library(ggplot2)

ggplot(nd2, aes(x = mean_breadth_0.5, y = mean_peak, fill = pred)) +
  geom_tile() +
  scale_fill_gradient2(low = "firebrick", mid = "white", high = "steelblue", midpoint = 0)

labs(title = "Heatmap", x = "Niche breadth", y = "HFI Peak of niche", fill = "Estimated population trend") +
  theme_minimal()
