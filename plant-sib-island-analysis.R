#### Set up: ####
rm(list = ls())

set.seed(123)

library(glmmTMB)
library(tidyverse)
library(MuMIn)
library(DHARMa)
library(MASS)
library(patchwork) 

dat <- read.csv("plant-sib-island-analysis-data.csv")

#### Fit model: ####
mod1 <- glmmTMB(rar_s ~ l.area.std * sqrt.wrack.std + l.area.std * soild15n.std + slope.std + dist_near.std + (1|node), data = dat, REML = FALSE)

# Look at summary:
summary(mod1)
r.squaredGLMM(mod1) # Marginal: 0.28, Conditional: 0.60

# Fit all possible subsets of the global model. Average across all of them to obtain model-averaged coefficients and RVIs.
mod_dredge <- dredge(mod1)
many_mods <- get.models(mod_dredge, subset = TRUE)
mod_avg <- model.avg(many_mods, fit = TRUE)

#### Check model diagnostics: ####
# Variance Inflation Factors: Refitting with lmer because this test doesn't work with glmmTMB models.
vif.test <- lme4::lmer(rar_s ~ l.area.std * sqrt.wrack.std + l.area.std * soild15n.std + slope.std + dist_near.std + (1|node),
                       data = dat,
                       REML = FALSE)
car::vif(vif.test) # These are fine.

# Model diagnostics (using DHARMa package):
# Simulate residuals
simulated_resids <- simulateResiduals(fittedModel = mod1)

# Residual check on simulated residuals
plot(simulated_resids) # This looks pretty good. No significant deviations.

# Check individual predictors
plotResiduals(simulated_resids, form = dat$l.area.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$dist_near.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$slope.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$sqrt.wrack.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$soild15n.std, quantreg = FALSE)

# Pearson residuals:
ggplot(data.frame(l.area.std = dat$l.area.std, pearson = residuals(mod1, type="pearson")),
       aes(x = l.area.std, y = pearson)) +
  geom_point() +
  theme_bw()

ggplot(data.frame(dist_near.std = dat$dist_near.std, pearson = residuals(mod1, type="pearson")),
       aes(x = dist_near.std, y = pearson)) +
  geom_point() +
  theme_bw()


ggplot(data.frame(slope.std = dat$slope.std, pearson = residuals(mod1, type="pearson")),
       aes(x = slope.std, y = pearson)) +
  geom_point() +
  theme_bw()

# This one looks a bit wonky but it is real data.
ggplot(data.frame(sqrt.wrack.std = dat$sqrt.wrack.std, pearson = residuals(mod1, type="pearson")),
       aes(x = sqrt.wrack.std, y = pearson)) +
  geom_point() +
  theme_bw()


ggplot(data.frame(soild15n.std = dat$soild15n.std, pearson = residuals(mod1, type="pearson")),
       aes(x = soild15n.std, y = pearson)) +
  geom_point() +
  theme_bw()


# Goodness of fit test
testUniformity(simulationOutput = simulated_resids) # Good

# Check for Outliers
testOutliers(simulationOutput = simulated_resids) # One outlier. Prob GS04 due to v high wrack.

# Check for overdispersion
testDispersion(simulated_resids) # No over/underdispersion detected

#### RVI: ####
MuMIn::importance(mod_avg)
summary(mod_avg)

#### Coefplot: ####
coefs <- coefTable(mod_avg, full = FALSE)
coefs <- as.data.frame(coefs)
coefs$parameters <- rownames(coefs)
coefs$df <- NULL
names(coefs) <- c("fit", "fit.se", "parameter")

params <- factor(c("cond(l.area.std)", "cond(slope.std)", "cond(dist_near.std)", "cond(sqrt.wrack.std)", "cond(soild15n.std)", "cond(l.area.std:sqrt.wrack.std)", "cond(l.area.std:soild15n.std)"))

coefs$parameter <- factor(coefs$parameter, levels = params)

coefs <- coefs %>%
  slice(match(params, parameter)) # This removes the intercept.

coefs$parameter <- factor(coefs$parameter,
                          levels = c("cond(slope.std)",
                                     "cond(sqrt.wrack.std)",
                                     "cond(dist_near.std)",
                                     "cond(soild15n.std)",
                                     "cond(l.area.std:soild15n.std)",
                                     "cond(l.area.std:sqrt.wrack.std)",
                                     "cond(l.area.std)"),
                          labels = c("Mean island slope",
                                     "Wrack biomass",
                                     "Distance to nearest landmass",
                                     "Soil d15N",
                                     "Area * d15n",
                                     "Area * wrack",
                                     "Island area"))

param.labs <- c("Mean island slope",
                "Wrack biomass",
                "Distance to nearest landmass",
                expression(paste("Forest-edge soil"~delta^{15}, "N")),
                expression(paste("Island area * forest-edge soil"~delta^{15}, "N")),
                "Island area * wrack biomass",
                "Island area")

coefplot <- ggplot(coefs, aes(x = parameter, y = fit)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_point(size = 2, col = "gray60") +
  geom_point(data = coefs[coefs$parameter == "Island area", ],
             col = "blue1",
             size = 2) +
  geom_point(data = coefs[coefs$parameter == "Mean island slope", ],
             col = "maroon",
             size = 2) +
  geom_errorbar(aes(ymin = (fit - 1.96 * fit.se),
                    ymax = (fit + 1.96 * fit.se),
                    width = 0),
                col = "grey60") +
  geom_errorbar(data = coefs[coefs$parameter == "Island area", ],
                col = "blue1",
                aes(ymin = (fit - 1.96 * fit.se),
                    ymax = (fit + 1.96 * fit.se),
                    width = 0)) +
  geom_errorbar(data = coefs[coefs$parameter == "Mean island slope", ],
                col = "maroon",
                aes(ymin = (fit - 1.96 * fit.se),
                    ymax = (fit + 1.96 * fit.se),
                    width = 0)) +
  geom_hline(yintercept = 0, col = "grey50", linetype = "dashed") +
  ylab("Coefficient estimate") +
  xlab("") +
  scale_x_discrete(labels = param.labs) +
  coord_flip();coefplot

# Figure:
# tiff("isl-coef-plot.tiff", units = "px", height = 2200, width = 3000, res = 600)
coefplot
# dev.off()

#### Biplots - area: ####
# Simulate the uncertainty in parameter estimates and generate 95% credible intervals (Code credit to Dr. Dan Greenberg)
coef_a <- c(mod1$fit$par[1:8]) # Extract the point estimates for coefficients from your model (eg. intercept and slope)

vcov_a <- vcov(mod1)[[1]] # Variance-covariance matrix of the parameter estimates

pars.resamp.a <- mvrnorm(500, mu = coef_a, Sigma = vcov_a) # Multivariate normal simulation of errors - gives the mean and the VCV matrix

x.a <- seq(min(dat$l.area.std), max(dat$l.area.std), length.out = 1000) # range of predictor

y.conf <- list() # Create an empty list for the estimate y values for each simulation run

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.a <-  expand.grid(x = x.a,
                             l.ci = NA,
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA)


# Calculate additional variance due to random effects & random noise:
summary(mod1)
sqrt(12.37 + 15.63) # 5.291503

for(z in 1:length(x.a)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[,2] * x.a[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500,pars.resamp.a[t,1] + pars.resamp.a[t,2]*x.a[z], 5.291503)   # add in model standard deviation

  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors=FALSE))
  pred.frame.a[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.a[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.a[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.a[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.a[z,6] <- coef_a[1] + coef_a[2] * x.a[z]
}

pred.frame.a <- pred.frame.a %>%
  rename(l.area.std = x, rar_s = m)

# write.csv(pred.frame.a, "isl-level-predictions-area.csv", row.names = FALSE)

xbreaks.area <- c(((log10(100)-mean(dat$l.area) / sd(dat$l.area))),
                  ((log10(1000)-mean(dat$l.area) / sd(dat$l.area))),
                  ((log10(10000)-mean(dat$l.area) / sd(dat$l.area))),
                  ((log10(100000)-mean(dat$l.area) / sd(dat$l.area))),
                  ((log10(1000000)-mean(dat$l.area) / sd(dat$l.area))))

# Plot with prediction interval:
area_plot <- ggplot(NULL, aes(y = rar_s, x = l.area.std)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  geom_ribbon(data = pred.frame.a, aes(ymin = l.ci, ymax = u.ci), alpha = 1, linetype = 0, fill = "grey60") +
  geom_ribbon(data = pred.frame.a, aes(ymin = l.pi, ymax = u.pi), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.a, aes(y = rar_s), colour = "black") +
  geom_point(data = dat, colour = "black", size = 2) +
  xlab(expression(~Island~area~(m^{"2"}))) +
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_continuous(breaks = xbreaks.area,
                     limits = c(-2.4, 2.5),
                     labels = c("100", "1 000", "10 000", "100 000", "1 000 000")) +
  ylab("Rarefied species richness");area_plot

# tiff("Output/isl-level-area-richness.tiff", units = "px", width = 4000, height = 2400, res = 600)
# area_plot
# dev.off()


#### Biplots - slope: ####
# Simulate the uncertainty in parameter estimates and generate 95% credible intervals
coef_a <- c(mod1$fit$par[1:8]) # Extract the point estimates for coefficients from your model (eg. intercept and slope)

vcov_a <- vcov(mod1)[[1]] # Variance-covariance matrix of the parameter estimates

pars.resamp.b <- mvrnorm(500, mu = coef_a, Sigma = vcov_a) # Multivariate normal simulation of errors - gives the mean and the VCV matrix

x.b <- seq(min(dat$slope.std), max(dat$slope.std), length.out = 1000) # range of predictor

y.conf <- list() # Reset empty list for the estimate y values for each simulation run

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.b <-  expand.grid(x = x.b,
                             l.ci = NA,
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA)


# Calculate additional variance due to random effects & random noise:
summary(mod1)
sqrt(12.37 + 15.63) # 5.291503

for(z in 1:length(x.b)){
  y.conf[[z]] <- pars.resamp.b[,1] + pars.resamp.b[,5] * x.b[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500, pars.resamp.b[t,1] + pars.resamp.b[t,5]*x.b[z], 5.291503)   # add in model standard deviation

  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors=FALSE))
  pred.frame.b[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.b[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.b[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.b[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.b[z,6] <- coef_a[1] + coef_a[5] * x.b[z]
}

pred.frame.b <- pred.frame.b %>%
  rename(slope.std = x, rar_s = m)

# write.csv(pred.frame.b, "isl-level-predictions-slope.csv", row.names = FALSE)

# Unstandardize slope:
unstandardize.slope <- function() {
  function(x) format(x*sd(dat$slope_mean) + mean(dat$slope_mean), digits = 1)
}

# Plot with PI:
slope_plot <- ggplot(NULL, aes(y = rar_s, x = slope.std)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  geom_ribbon(data = pred.frame.b, aes(ymin = l.ci, ymax = u.ci), alpha = 1, linetype = 0, fill = "grey60") +
  geom_ribbon(data = pred.frame.b, aes(ymin = l.pi, ymax = u.pi), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.b, aes(y = rar_s), colour = "black") +
  geom_point(data = dat, colour = "black", size = 2) +
  xlab("Mean island slope") +
  scale_x_continuous(labels = unstandardize.slope(),
                     breaks = c(-2, -1.2, -0.3, 0.6, 1.5, 2.4)) +
  scale_y_continuous(limits = c(0, 40)) +
  ylab("Rarefied species richness");slope_plot

# tiff("Output/isl-level-slope-richness.tiff", units = "px", width = 4000, height = 2400, res = 600)
# slope_plot
# dev.off()

# Island-level plot for publication using patchwork:
patchwork = coefplot + area_plot + slope_plot

# Remove title from third subplot, add label:
patchwork[[3]] = patchwork[[3]] + theme(axis.text.y = element_blank(),
                                        axis.ticks.y = element_blank(),
                                        axis.title.y = element_blank()) +
  annotate("text", x = -2.7, y = 39, label = "c", size = 6)

patchwork[[1]] = patchwork[[1]] + theme(axis.text = element_text(size = 10),
                                        axis.title = element_text(size = 14)) +
  annotate("text", x = 7, y = -2.9, label = "a", size = 6)

patchwork[[2]] = patchwork[[2]] + annotate("text", x = -2.3, y = 39, label = "b", size = 6)

patchwork

layout <- "
AABBCC
AABBCC
AABBCC
"
patchwork

# tiff("3-panel-isl-level-dec20.tiff", units = "px", width = 8000, height = 2400, res = 600)
# patchwork
# dev.off()

#### Predictions based on GLOBAL model:####
# Area:
newdat_avg_isl <- expand.grid(l.area.std = median(dat$l.area.std),
                              slope.std = mean(dat$slope.std),
                              soild15n.std = mean(dat$soild15n.std),
                              dist_near.std = mean(dat$dist_near.std),
                              sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                              node = NA)

pred_avg_isl <- predict(mod1, newdat_avg_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_avg_isl$se.fit * 1.96
# 17.4 +/- 2.5

median(10^dat$l.area) # 13108
median(dat$area) # 13108

newdat_large_isl <- expand.grid(l.area.std = median(dat$l.area.std) + sd(dat$l.area.std),
                                slope.std = mean(dat$slope.std),
                                soild15n.std = mean(dat$soild15n.std),
                                dist_near.std = mean(dat$dist_near.std),
                                sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                                node = NA)

pred_large_isl <- predict(mod1, newdat_large_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_large_isl[[2]] * 1.96
# 19.7 +/- 2.8

# Size:
10^(median(dat$l.area) + sd(dat$l.area)) # 126000

newdat_small_isl <- expand.grid(l.area.std = median(dat$l.area.std) - sd(dat$l.area.std),
                                slope.std = mean(dat$slope.std),
                                soild15n.std = mean(dat$soild15n.std),
                                dist_near.std = mean(dat$dist_near.std),
                                sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                                node = NA)

pred_small_isl <- predict(mod1, newdat_small_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_small_isl[[2]] * 1.96

# 15.1 +/- 2.8

10^(median(dat$l.area) - sd(dat$l.area)) # 1362

# Slope:
newdat_slope_isl <- expand.grid(l.area.std = median(dat$l.area.std),
                                slope.std = mean(dat$slope.std),
                                soild15n.std = mean(dat$soild15n.std),
                                dist_near.std = mean(dat$dist_near.std),
                                sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                                node = NA)

pred_avg_slope_isl <- predict(mod1, newdat_slope_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_avg_slope_isl$se.fit * 1.96
# 17.4 +/- 2.5

mean(dat$slope_mean) # 21.3
median(dat$slope_mean) # 21.7

hist(dat$slope_mean)

# 30 degree island

# This is ~ 30 degrees:
mean(dat$slope_mean) + 1.6*sd(dat$slope_mean) # 30.25

newdat_steep_isl <- expand.grid(l.area.std = median(dat$l.area.std),
                                slope.std = mean(dat$slope.std) + 1.6*sd(dat$slope.std),
                                soild15n.std = mean(dat$soild15n.std),
                                dist_near.std = mean(dat$dist_near.std),
                                sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                                node = NA)

pred_steep_isl <- predict(mod1, newdat_steep_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_steep_isl[[2]] * 1.96

# 14.3 +/- 3.3

# Okay, now an island ~ 15 degrees:
mean(dat$slope_mean) - 1.1*sd(dat$slope_mean) # 15.09

newdat_shallow_isl <- expand.grid(l.area.std = median(dat$l.area.std),
                                  slope.std = mean(dat$slope.std) - 1.1*sd(dat$slope.std),
                                  soild15n.std = mean(dat$soild15n.std),
                                  dist_near.std = mean(dat$dist_near.std),
                                  sqrt.wrack.std = mean(dat$sqrt.wrack.std),
                                  node = NA)

pred_shallow_isl <- predict(mod1, newdat_shallow_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_shallow_isl[[2]] * 1.96

# 19.5 +/- 3.0
