#### Set up: ####
rm(list = ls())
set.seed(123)

library(glmmTMB)
library(tidyverse)
library(MuMIn)
library(DHARMa)
library(MASS)
library(patchwork)

dat <- read.csv("plant-sib-plot-analysis-data.csv")
  
#### Fit model: ####
mod1 <- glmmTMB(s ~ l.area.std * shoredist.std + dist_near.std + sqrt.sum.std + fs_pc1.std + sm_av.std + slope.std + avg.n.std + avg.d15n.std + (1|node/island/unq_tran), 
                data = dat,
                family = poisson(),
                na.action = "na.fail",
                REML = FALSE)

summary(mod1)
r.squaredGLMM(mod1) 

# Note: This takes a long time.
mod_dredge <- dredge(mod1, trace = 2)
many_mods <- get.models(mod_dredge, subset = TRUE)
mod_avg <- model.avg(many_mods, fit = TRUE) 

#### Model diagnostics: ####
performance::check_collinearity(mod1) # These are fine.

# Simulate residuals
simulated_resids <- simulateResiduals(fittedModel = mod1)

# Residual check on simulated residuals 
plot(simulated_resids) # Not perfect but not the worst. Significant deviation on the low end but a ton of data. Looks fine overall.

# Check individual predictors
plotResiduals(simulated_resids, form = dat$l.area.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$dist_near.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$slope.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$sqrt.wrack.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$shoredist.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$sm_av.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$fs_pc1.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$avg.n.std, quantreg = FALSE)
plotResiduals(simulated_resids, form = dat$avg.d15n.std, quantreg = FALSE)

# Goodness of fit test
testUniformity(simulationOutput = simulated_resids) # Significant deviation on low end.

# Check for Outliers
testOutliers(simulationOutput = simulated_resids) # 3 outlier out of 1383 points. not affecting model.

# Check for overdispersion
testDispersion(simulated_resids) # SLIGHT overdispersion (p = 0.048)

#### Coefficient plot: ####
# Coefplot: 
coefs <- coefTable(mod_avg, full = FALSE)
coefs <- as.data.frame(coefs)
coefs$parameters <- rownames(coefs)
coefs$df <- NULL
names(coefs) <- c("fit", "fit.se", "parameter")

params <- factor(c("cond(l.area.std)", "cond(shoredist.std)", "cond(dist_near.std)", "cond(sqrt.sum.std)", "cond(fs_pc1.std)", "cond(sm_av.std)", "cond(slope.std)", "cond(avg.n.std)", "cond(avg.d15n.std)", "cond(l.area.std:shoredist.std)"))

coefs$parameter <- factor(coefs$parameter, levels = params)

coefs <- coefs %>%
  slice(match(params, parameter)) # This removes the intercept.

param.labs <- c("Plot d15N", "Forest Structure PC", "Island area", "Distance to shore", "Plot slope", "Soil moisture", "Total N")


coefs$parameter <- factor(coefs$parameter, 
                          levels = c("cond(shoredist.std)", 
                                     "cond(avg.d15n.std)", 
                                     "cond(avg.n.std)", 
                                     "cond(sqrt.sum.std)", 
                                     "cond(l.area.std:shoredist.std)",
                                     "cond(dist_near.std)", 
                                     "cond(l.area.std)", 
                                     "cond(sm_av.std)", 
                                     "cond(slope.std)", 
                                     "cond(fs_pc1.std)"),
                          labels = c("Distance to shore",
                                     "Average d15N",
                                     "Average total N",
                                     "Wrack biomass",
                                     "Island area * distance to shore",
                                     "Distance to nearest landmass",
                                     "Island area",
                                     "Soil moisture", 
                                     "Plot slope",
                                     "Forest openness"))

coefplot <- ggplot(coefs, aes(x = parameter, y = fit)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_errorbar(aes(ymin = (fit - 1.96 * fit.se), 
                    ymax = (fit + 1.96 * fit.se), 
                    width = 0)) +
  geom_hline(yintercept = 0, col = "grey50", linetype = "dashed", lwd = 1) +
  geom_point(size = 2) +  
  ylab("Coef. Est.") +
  xlab("") +
  scale_x_discrete(name="",
                   breaks = c("Distance to shore",
                              "Average d15N",
                              "Average total N",
                              "Wrack biomass",
                              "Island area * distance to shore",
                              "Distance to nearest landmass",
                              "Island area",
                              "Soil moisture",
                              "Plot slope",
                              "Forest openness"),
                   labels = c("Distance to shore",
                              expression(paste("\n", " ", "Average", " ", delta^{15}, "N")),
                              "Average %N",
                              "Wrack biomass",
                              "Island area * distance to shore",
                              "Distance to nearest landmass",
                              "Island area",
                              "Soil moisture",
                              "Plot slope",
                              "Forest openness")) +
  coord_flip();coefplot

coefplot_c <- ggplot(coefs, aes(x = parameter, y = fit)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_errorbar(aes(ymin = (fit - 1.96 * fit.se), 
                    ymax = (fit + 1.96 * fit.se), 
                    width = 0), 
                col = "grey60") +
  geom_hline(yintercept = 0, col = "grey50", linetype = "dashed", lwd = 1) +
  geom_point(size = 2, col = "grey60") +  
  geom_errorbar(data = coefs[coefs$parameter %in% c("Forest openness", "Plot slope", "Soil moisture", "Island area"), ],
                aes(ymin = (fit - 1.96 * fit.se), 
                    ymax = (fit + 1.96 * fit.se), 
                    width = 0), 
                col = "blue1") +
  geom_errorbar(data = coefs[coefs$parameter %in% c("Average d15N", "Distance to shore"), ],
                aes(ymin = (fit - 1.96 * fit.se), 
                    ymax = (fit + 1.96 * fit.se), 
                    width = 0), 
                col = "maroon") +
  geom_point(data = coefs[coefs$parameter %in% c("Forest openness", "Plot slope", "Soil moisture", "Island area"), ],
             col = "blue1", 
             size = 2) + 
  geom_point(data = coefs[coefs$parameter %in% c("Average d15N", "Distance to shore"), ],
             col = "maroon", 
             size = 2) + 
  ylab("Coefficient estimate") +
  xlab("") +
  scale_x_discrete(name="",
                   breaks = c("Distance to shore",
                              "Average d15N",
                              "Average total N",
                              "Wrack biomass",
                              "Island area * distance to shore",
                              "Distance to nearest landmass",
                              "Island area",
                              "Soil moisture",
                              "Plot slope",
                              "Forest openness"),
                   labels = c("Distance to shore",
                              expression(paste("\n", " ", "Average", " ", delta^{15}, "N")),
                              "Average %N",
                              "Wrack biomass",
                              "Island area * distance to shore",
                              "Distance to nearest landmass",
                              "Island area",
                              "Soil moisture",
                              "Plot slope",
                              "Forest openness")) +
  coord_flip();coefplot_c

tiff("Output/plot-level-coefplot-with-transect-level-nutrients.tiff", units = "px", height = 2200, width = 2400, res = 600)
coefplot
dev.off()

# 
#### Bootstrapping CIs and PIs: ####
# Pull out model coefficients for all: 
# Simulate the uncertainty in parameter estimates and generate 95% credible intervals
coef_a <- c(mod1$fit$par[1:11]) # Extract the point estimates for coefficients from your model (eg. intercept and slope)

vcov_a <- vcov(mod1)[[1]] # Variance-covariance matrix of the parameter estimates

pars.resamp.a <- mvrnorm(500, mu = coef_a, Sigma = vcov_a) # Multivariate normal simulation of errors - gives the mean and the VCV matrix

# Area:
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
sqrt(0.000163 + 0.0205259 + 0.0120953) # 0.1810641

for(z in 1:length(x.a)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[,2] * x.a[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500, pars.resamp.a[t,1] + pars.resamp.a[t,2]*x.a[z], 0.1810641)   # add in model standard deviation
    
  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors = FALSE))
  pred.frame.a[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.a[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.a[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.a[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.a[z,6] <- coef_a[1] + coef_a[2] * x.a[z]
}

pred.frame.a <- pred.frame.a %>% 
  rename(l.area.std = x, s = m)

# write.csv(pred.frame.a, "plot-level-predictions-area-apr0721.csv", row.names = FALSE)

# Soil moisture: 
x.b <- seq(min(dat$sm_av.std), max(dat$sm_av.std), length.out = 1000) # range of predictor 

y.conf <- list() # Reset for simulation:

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.b <-  expand.grid(x = x.b,
                             l.ci = NA, 
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA) 

# Calculate additional variance due to random effects & random noise: 
summary(mod1)
sqrt(0.000163 + 0.0205259 + 0.0120953) # 0.1810641

for(z in 1:length(x.b)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[, 7] * x.b[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500, pars.resamp.a[t,1] + pars.resamp.a[t, 7]*x.b[z], 0.1810641)   # add in model standard deviation
    
  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors = FALSE))
  pred.frame.b[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.b[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.b[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.b[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.b[z,6] <- coef_a[1] + coef_a[7] * x.b[z]
}

pred.frame.b <- pred.frame.b %>% 
  rename(sm_av.std = x, s = m)

# write.csv(pred.frame.b, "plot-level-predictions-soil-moisture-apr0721.csv", row.names = FALSE)

# Plot slope: 
x.c <- seq(min(dat$slope.std), max(dat$slope.std), length.out = 1000) # range of predictor 

y.conf <- list() # Create an empty list for the estimate y values for each simulation run

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.c <-  expand.grid(x = x.c,
                             l.ci = NA, 
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA) 


# Calculate additional variance due to random effects & random noise: 
summary(mod1)
sqrt(0.000163 + 0.0205259 + 0.0120953) # 0.1810641

for(z in 1:length(x.c)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[,8] * x.c[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500,pars.resamp.a[t,1] + pars.resamp.a[t,8]*x.c[z], 0.1810641)   # add in model standard deviation
    
  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors=FALSE))
  pred.frame.c[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.c[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.c[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.c[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.c[z,6] <- coef_a[1] + coef_a[8] * x.c[z]
}

pred.frame.c <- pred.frame.c %>% 
  rename(slope.std = x, s = m)

# write.csv(pred.frame.c, "plot-level-predictions-slope-apr0721.csv", row.names = FALSE)

# Distance to shore: 
x.d <- seq(min(dat$shoredist.std), max(dat$shoredist.std), length.out = 1000) # range of predictor 

y.conf <- list() # Create an empty list for the estimate y values for each simulation run

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.d <-  expand.grid(x = x.d,
                             l.ci = NA, 
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA) 


# Calculate additional variance due to random effects & random noise: 
summary(mod1)
sqrt(0.000163 + 0.0205259 + 0.0120953) # 0.1810641

for(z in 1:length(x.d)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[,3] * x.d[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500,pars.resamp.a[t,1] + pars.resamp.a[t,3]*x.d[z], 0.1810641)   # add in model standard deviation
    
  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors=FALSE))
  pred.frame.d[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.d[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.d[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.d[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.d[z,6] <- coef_a[1] + coef_a[3] * x.d[z]
}

pred.frame.d <- pred.frame.d %>% 
  rename(shoredist.std = x, s = m)

# write.csv(pred.frame.d, "plot-level-predictions-shoredist-apr0721.csv", row.names = FALSE)

# Forest structure: 
x.e <- seq(min(dat$fs_pc1.std), max(dat$fs_pc1.std), length.out = 1000) # range of predictor 

y.conf <- list() # Create an empty list for the estimate y values for each simulation run

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.e <-  expand.grid(x = x.e,
                             l.ci = NA, 
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA) 


# Calculate additional variance due to random effects & random noise: 
summary(mod1)
sqrt(0.000163 + 0.0205259 + 0.0120953) # 0.1810641

for(z in 1:length(x.e)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[,6] * x.e[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500,pars.resamp.a[t,1] + pars.resamp.a[t,6]*x.e[z], 0.1810641)   # add in model standard deviation
    
  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors=FALSE))
  pred.frame.e[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.e[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.e[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.e[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.e[z,6] <- coef_a[1] + coef_a[6] * x.e[z]
}

pred.frame.e <- pred.frame.e %>% 
  rename(fs_pc1.std = x, s = m)

# write.csv(pred.frame.e, "plot-level-predictions-fs-apr0721.csv", row.names = FALSE)

# Soil d15N:
x.f <- seq(min(dat$avg.d15n.std), max(dat$avg.d15n.std), length.out = 1000) # range of predictor 

y.conf <- list() # Create an empty list for the estimate y values for each simulation run

# Create dataframe for each value of X, the limits on y from the simulation
pred.frame.f <-  expand.grid(x = x.f,
                             l.ci = NA, 
                             u.ci = NA,
                             l.pi = NA,
                             u.pi = NA,
                             m = NA) 


# Calculate additional variance due to random effects & random noise: 
summary(mod1)
sqrt(0.000163 + 0.0205259 + 0.0120953) # 0.1810641

for(z in 1:length(x.e)){
  y.conf[[z]] <- pars.resamp.a[,1] + pars.resamp.a[,10] * x.f[z]  # fit the estimate based on the sampled intercept and slopes
  pred.temp <- list()
  for(t in 1:500){
    pred.temp[[t]]<- rnorm(500,pars.resamp.a[t,1] + pars.resamp.a[t,10]*x.f[z], 0.1810641)   # add in model standard deviation
    
  }
  pred.temp.all<- do.call(rbind, lapply(pred.temp, data.frame, stringsAsFactors=FALSE))
  pred.frame.f[z,2]<- quantile(y.conf[[z]], 0.025)
  pred.frame.f[z,3]<- quantile(y.conf[[z]], 0.975)
  pred.frame.f[z,4]<- quantile(pred.temp.all[,1], 0.025)
  pred.frame.f[z,5]<- quantile(pred.temp.all[,1], 0.975)
  pred.frame.f[z,6] <- coef_a[1] + coef_a[10] * x.f[z]
}

pred.frame.f <- pred.frame.f %>% 
  rename(avg.d15n.std = x, s = m)

# write.csv(pred.frame.f, "plot-level-predictions-d15n-apr0721.csv", row.names = FALSE)

#### Bi-plots: ####
xbreaks.area <- c(log10(100),
                  log10(1000),
                  log10(10000),
                  log10(100000),
                  log10(1000000))

pred.frame.a$l.area <- pred.frame.a$l.area.std * sd(dat$l.area) + mean(dat$l.area)
pred.frame.a$area <- 10^pred.frame.a$l.area
dat$area <- 10^dat$l.area

area_plot <- ggplot(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10)) +
  geom_ribbon(data = pred.frame.a, aes(y = exp(s), x = l.area, ymin = exp(l.ci), ymax = exp(u.ci)), alpha = 1, linetype = 0, fill = "grey60") + 
  geom_ribbon(data = pred.frame.a, aes(y = exp(s), x = l.area, ymin = exp(l.pi), ymax = exp(u.pi)), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.a, aes(y = exp(s), x = l.area), colour = "black") +
  geom_point(data = dat, aes(y = s, x = l.area), colour = "black", alpha = 0.1) + 
  xlab(expression(~Island~area~(m^{"2"}))) +
  scale_x_continuous(breaks = xbreaks.area,
                     labels = c("100", "1 000", "10 000", "100 000", "1 000 000")) +
  ylab("Species richness");area_plot

pred.frame.b$sm_av <- pred.frame.b$sm_av.std*sd(dat$sm_av) + mean(dat$sm_av)

sm_plot <- ggplot(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10)) +
  geom_ribbon(data = pred.frame.b, aes(y = exp(s), x = sm_av, ymin = exp(l.ci), ymax = exp(u.ci)), alpha = 1, linetype = 0, fill = "grey60") + 
  geom_ribbon(data = pred.frame.b, aes(y = exp(s), x = sm_av, ymin = exp(l.pi), ymax = exp(u.pi)), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.b, aes(y = exp(s), x = sm_av), colour = "black") +
  geom_point(data = dat, aes(y = s, x = sm_av), colour = "black", alpha = 0.1) + 
  xlab("Soil moisture (%)") +
  ylab("Species richness") +
  xlim(c(0, 82));sm_plot


pred.frame.c$slope <- pred.frame.c$slope.std * sd(dat$slope) + mean(dat$slope)

slope_plot <- ggplot(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10)) +
  geom_ribbon(data = pred.frame.c, aes(y = exp(s), x = slope, ymin = exp(l.ci), ymax = exp(u.ci)), alpha = 1, linetype = 0, fill = "grey60") + 
  geom_ribbon(data = pred.frame.c, aes(y = exp(s), x = slope, ymin = exp(l.pi), ymax = exp(u.pi)), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.c, aes(y = exp(s), x = slope), colour = "black") +
  geom_point(data = plot_analysis, aes(y = s, x = slope), colour = "black", alpha = 0.1) +
  xlab("Plot slope (%)") +
  ylab("Species richness");slope_plot

pred.frame.d$shoredist <- pred.frame.d$shoredist.std * sd(dat$shore_dist) + mean(dat$shore_dist)

shoredist_plot <- ggplot(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10)) +
  geom_ribbon(data = pred.frame.d, aes(y = exp(s), x = shoredist, ymin = exp(l.ci), ymax = exp(u.ci)), alpha = 1, linetype = 0, fill = "grey60") + 
  geom_ribbon(data = pred.frame.d, aes(y = exp(s), x = shoredist, ymin = exp(l.pi), ymax = exp(u.pi)), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.d, aes(y = exp(s), x = shoredist), colour = "black") +
  geom_point(data = dat, aes(y = s, x = shore_dist), colour = "black", alpha = 0.1) + 
  xlab("Distance to shore (m)") + 
  ylab("Species richness") + geom_jitter(data = dat, aes(y = s, x = shore_dist), width = 1.4, alpha = 0.1 );shoredist_plot

pred.frame.e$fs_pc1 <- pred.frame.e$fs_pc1.std * sd(dat$fs_pc1) + mean(dat$fs_pc1)

fs_plot <- ggplot(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10)) +
  geom_ribbon(data = pred.frame.e, aes(y = exp(s), x = fs_pc1, ymin = exp(l.ci), ymax = exp(u.ci)), alpha = 1, linetype = 0, fill = "grey60") + 
  geom_ribbon(data = pred.frame.e, aes(y = exp(s), x = fs_pc1, ymin = exp(l.pi), ymax = exp(u.pi)), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.e, aes(y = exp(s), x = fs_pc1), colour = "black") +
  geom_point(data = dat, aes(y = s, x = fs_pc1), colour = "black", alpha = 0.1) +
  xlab("Forest openness") +
  ylab("Species richness");fs_plot

pred.frame.f$avg.d15n <- pred.frame.f$avg.d15n.std * sd(dat$avg.d15n) + mean(dat$avg.d15n)

d15n_plot <- ggplot(NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 10)) +
  geom_ribbon(data = pred.frame.f, aes(y = exp(s), x = avg.d15n, ymin = exp(l.ci), ymax = exp(u.ci)), alpha = 1, linetype = 0, fill = "grey60") + 
  geom_ribbon(data = pred.frame.f, aes(y = exp(s), x = avg.d15n, ymin = exp(l.pi), ymax = exp(u.pi)), alpha = .1, linetype = 0, fill = "grey30") +
  geom_line(data = pred.frame.f, aes(y = exp(s), x = avg.d15n), colour = "black") +
  geom_point(data = dat, aes(y = s, x = avg.d15n), colour = "black", alpha = 0.1) +
  xlab((expression(paste("\n", "Average soil", " ", delta^{15}, "N")))) +
  ylab("Species richness");d15n_plot

#### Final plot for publication: ####
area_plot <- area_plot + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y = element_blank()) 

shoredist_plot <- shoredist_plot + theme(axis.text.y = element_blank(),
                                         axis.ticks.y = element_blank(),
                                         axis.title.y = element_blank())

slope_plot <- slope_plot + theme(axis.text.y = element_blank(),
                                 axis.ticks.y = element_blank(),
                                 axis.title.y = element_blank())

patchwork = (area_plot / fs_plot) | (d15n_plot / sm_plot)| (shoredist_plot / slope_plot)

tiff("Output/plot-level-figs.tiff", units = "px", width = 5000, height = 3200, res = 600)
patchwork
dev.off()

coefplot_c2 = coefplot_c + annotate("text", x = 10.2, y = -0.1, label = "a", size = 5)
fs_plot2 = fs_plot + annotate("text", x = -6, y = 20, label = "b", size = 5)
slope_plot2 = slope_plot + annotate("text", x = 5, y = 20, label = "c", size = 5)
sm_plot2 = sm_plot + annotate("text", x = 3, y = 20, label = "d", size = 5)
area_plot2 = area_plot + annotate("text", x = 2.2, y = 20, label = "e", size = 5)
d15n_plot2 = d15n_plot + annotate("text", x = -1, y = 20, label = "f", size = 5)
shoredist_plot2 = shoredist_plot + annotate("text", x = 1.5, y = 20, label = "g", size = 5)

patchwork2 = coefplot_c2 | (fs_plot2 / area_plot2) | (slope_plot2 / d15n_plot2)| (sm_plot2 / shoredist_plot2)


tiff("Output/plot-level-figs-coefplot-dec2021.tiff", units = "px", width = 7500, height = 3000, res = 600)
patchwork2
dev.off()

patchwork3 = coefplot_c2 | (fs_plot2 / sm_plot2 / d15n_plot2) | (slope_plot2 / area_plot2 / shoredist_plot2)

A <- coefplot_c2
B <- fs_plot2
C <- slope_plot2
D <- sm_plot2 + ylab("Species richness")
E <- area_plot2
G <- d15n_plot2
H <- shoredist_plot2

layout <- "
ABBEE
ACCGG
#DDHH
"

tiff("Output/plot-level-figs-coefplot-dec2021-long.tiff", units = "px", width = 6000, height = 4500, res = 600)
# A + B + C + D + E + H + G + plot_layout(design = "layout")
A / plot_spacer() | (B / D / G) | (C / E / H)
dev.off()

# Colour plots: 
sm_plot_c <- sm_plot_c + theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.title.y = element_blank()) 

d15n_plot_c <- d15n_plot_c + theme(axis.text.y = element_blank(),
                                   axis.ticks.y = element_blank(),
                                   axis.title.y = element_blank())


shoredist_plot_c <- shoredist_plot_c + theme(axis.text.y = element_blank(),
                                             axis.ticks.y = element_blank(),
                                             axis.title.y = element_blank())

slope_plot_c <- slope_plot_c + theme(axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     axis.title.y = element_blank())

patchwork3 = coefplot_c | (fs_plot_c / area_plot_c) | (slope_plot_c / d15n_plot_c)| (sm_plot_c / shoredist_plot_c)

tiff("Output/plot-level-figs-coefplot_c3.tiff", units = "px", width = 7500, height = 3000, res = 600)
patchwork3
dev.off()


#### Predictions: ####
# Island area: 
newdat_large_isl <- expand.grid(l.area.std = median(dat$l.area.std) + sd(dat$l.area.std),
                                dist_near.std = mean(dat$dist_near.std),
                                slope.std = mean(dat$slope.std),
                                sm_av.std = mean(dat$sm_av.std),
                                fs_pc1.std = mean(dat$fs_pc1.std),
                                sqrt.sum.std = mean(dat$sqrt.sum.std),
                                shoredist.std = mean(dat$shoredist.std),
                                avg.n.std = mean(dat$avg.n.std),
                                avg.d15n.std = mean(dat$avg.d15n.std),
                                node = NA,
                                island = NA,
                                unq_tran = NA)

pred_large_isl <- predict(mod1, newdat_large_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_large_isl[[2]] * 1.96
# 6.3 +/- 0.6

10^(median(dat$l.area) + sd(dat$l.area)) # 179292.9

newdat_small_isl <- expand.grid(l.area.std = median(dat$l.area.std) - sd(dat$l.area.std),
                                dist_near.std = mean(dat$dist_near.std),
                                slope.std = mean(dat$slope.std),
                                sm_av.std = mean(dat$sm_av.std),
                                fs_pc1.std = mean(dat$fs_pc1.std),
                                sqrt.sum.std = mean(dat$sqrt.sum.std),
                                shoredist.std = mean(dat$shoredist.std),
                                avg.n.std = mean(dat$avg.n.std),
                                avg.d15n.std = mean(dat$avg.d15n.std),
                                node = NA,
                                island = NA,
                                unq_tran = NA)

pred_small_isl <- predict(mod1, newdat_small_isl, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_small_isl[[2]] * 1.96

# 5.6 +/- 0.5

# Shoredist: 
newdat_shoreside <- expand.grid(shoredist.std = min(dat$shoredist.std),
                              dist_near.std = mean(dat$dist_near.std),
                              slope.std = mean(dat$slope.std),
                              sm_av.std = mean(dat$sm_av.std),
                              fs_pc1.std = mean(dat$fs_pc1.std),
                              sqrt.sum.std = mean(dat$sqrt.sum.std),
                              l.area.std = mean(dat$l.area.std),
                              avg.n.std = mean(dat$avg.n.std),
                              avg.d15n.std = mean(dat$avg.d15n.std),
                              node = NA,
                              island = NA,
                              unq_tran = NA)

pred_shoreside <- predict(mod1, newdat_shoreside, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_shoreside$se.fit * 1.96
# 6.7 +/- 0.6

newdat_inland <- expand.grid(shoredist.std = max(dat$shoredist.std) + sd(dat$l.area.std),
                                dist_near.std = mean(dat$dist_near.std),
                                slope.std = mean(dat$slope.std),
                                sm_av.std = mean(dat$sm_av.std),
                                fs_pc1.std = mean(dat$fs_pc1.std),
                                sqrt.sum.std = mean(dat$sqrt.sum.std),
                                l.area.std = mean(dat$l.area.std),
                                avg.n.std = mean(dat$avg.n.std),
                                avg.d15n.std = mean(dat$avg.d15n.std),
                                node = NA,
                                island = NA,
                                unq_tran = NA)

pred_inland <- predict(mod1, newdat_inland, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_inland[[2]] * 1.96
# 4.7 +/- 0.5

# d15N: 
newdat_low_d15N <- expand.grid(avg.d15n.std = median(dat$avg.d15n.std) - sd(dat$avg.d15n.std),
                                dist_near.std = mean(dat$dist_near.std),
                                slope.std = mean(dat$slope.std),
                                sm_av.std = mean(dat$sm_av.std),
                                fs_pc1.std = mean(dat$fs_pc1.std),
                                sqrt.sum.std = mean(dat$sqrt.sum.std),
                                shoredist.std = mean(dat$shoredist.std),
                                avg.n.std = mean(dat$avg.n.std),
                                l.area.std = mean(dat$l.area.std),
                                node = NA,
                                island = NA,
                                unq_tran = NA)

pred_low_d15N <- predict(mod1, newdat_low_d15N, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_low_d15N[[2]] * 1.96

# 6.3 +/- 0.6

newdat_high_d15N <- expand.grid(avg.d15n.std = median(dat$avg.d15n.std) + sd(dat$avg.d15n.std),
                               dist_near.std = mean(dat$dist_near.std),
                               slope.std = mean(dat$slope.std),
                               sm_av.std = mean(dat$sm_av.std),
                               fs_pc1.std = mean(dat$fs_pc1.std),
                               sqrt.sum.std = mean(dat$sqrt.sum.std),
                               shoredist.std = mean(dat$shoredist.std),
                               avg.n.std = mean(dat$avg.n.std),
                               l.area.std = mean(dat$l.area.std),
                               node = NA,
                               island = NA,
                               unq_tran = NA)

pred_high_d15N <- predict(mod1, newdat_high_d15N, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_high_d15N[[2]] * 1.96

# 5.6 +/- 0.5
