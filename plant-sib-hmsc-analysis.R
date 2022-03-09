#### Set up: ####
rm(list = ls())

set.seed(123)

library(Hmsc)
library(tidyverse)
library(corrplot)
library(knitr)
library(ape)
library(MASS)
library(fields)
library(greekLetters)

# Read in data: 
plants <- read.csv("plant-sib-percent-cover-long.csv")
env <- read.csv("plant-sib-plot-analysis-data.csv")
xycoords <- read.csv("plant-sib-xycoords2.csv")
names <- read.csv("plant-sib-species-names.csv")

plants2 <- read.csv("site-by-species-test.csv", row.names = 1)


#### Set up: env matrix ####
# Environmental matrix: 
env <- env %>% 
  dplyr::select(node, island, unq_tran, unq_plot, l.area.std, dist_near.std, slope.std, sm_av.std, fs_pc1.std, sqrt.sum.std, shoredist.std, avg.n.std, avg.d15n.std)

# Making the row name the plot id: 
row.names(env) <- env$unq_plot

# Call this X to follow suit with HMSC terminology.
X <- env %>% 
  dplyr::select(l.area.std, dist_near.std, slope.std, sm_av.std, fs_pc1.std, sqrt.sum.std, shoredist.std, avg.n.std, avg.d15n.std)


#### Set up species matrix: ####
names <- names %>% 
  dplyr::select(species, species_bin)

plants_key <- xycoords %>% 
  dplyr::select(old_unq_tran, unq_tran, unq_plot) %>% 
  distinct()

plants2 <- merge(plants, plants_key, by = "unq_plot")

# This results in 7992 observations.
plants2 <- plants2 %>% 
  dplyr::select(unq_isl, old_unq_tran, unq_tran, unq_plot, species, pcover) %>% 
  dplyr::filter(!species %in% c("l", "doje")) %>% # l is litter. Doje is not in the list of species and it's not common enough to make it in to the analysis anyways. 
  dplyr::group_by(unq_isl, unq_tran, unq_plot, species) %>% 
  dplyr::summarize(n = n())

# This results in 7148 observations. 
plants2 <- plants2 %>% 
  dplyr::filter(unq_plot %in% env$unq_plot) # Only keep plots that have accompanying environmental data.

# How many plant species?
length(unique(plants2$species)) # 94

# Figure out which species are present in over 5% of plots to keep for analysis.
# Number of plots: 1381
length(unique(plants2$unq_plot))
length(unique(plants2$unq_plot)) * 0.05 # ~70 plots. Plants need to be in at least 70 plots.

species_to_keep <- plants2 %>% 
  dplyr::group_by(unq_plot, species) %>% 
  dplyr::select(unq_plot, species)

# This tells you how many plots a species occurred in.
n.plots <- species_to_keep %>% 
  dplyr::group_by(species) %>% 
  distinct() %>% 
  dplyr::summarize(n = n())

# Put these in order from smallest number to largest:
n.plots <- n.plots[order(n.plots$n), ]

# Species to use in analysis need to be present in at least 70 plots.
species_to_keep <- n.plots %>% 
  dplyr::group_by(species) %>% 
  dplyr::filter(n >= 70)

# This leaves 18 species for analysis.

# Only include plants present in at least 70 quadrats: 
plants2 <- plants2 %>% 
  dplyr::filter(species %in% species_to_keep$species)

# lichen, m.
setdiff(unique(plants2$species), unique(names$species))

# Make these names consistent
names <- rbind(names, c("m", "moss"))
names <- rbind(names, c("lichen", "lichen"))

# Change names from 4 letter codes to full latin name: 
plants2 <- merge(plants2, names, by.x = "species", by.y = "species") # 6141 obs

plants2$species <- NULL
plants2$unq_isl <- NULL
plants2$unq_tran <- NULL

length(unique(plants2$species_bin))

plants2$species_bin <- as.factor(plants2$species_bin)

capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

plants2$species_bin <- capFirst(plants2$species_bin)

sum(plants2$n) # 6141 plants in 1326 plots, 18 species

plants3 <- plants2 %>% 
  dplyr::group_by(species_bin) %>% 
  dplyr::summarize(n = n())


# Change structure to one row per island, one species per column
reshape_plants <- plants2 %>% spread(species_bin, n)

# Change all NAs to 0s
reshape_plants[is.na(reshape_plants)] <- 0

# Making the row names the names of the islands
row.names(reshape_plants) <- reshape_plants$unq_plot

# Removing the column of island names
reshape_plants$unq_plot <- NULL


Y <- as.matrix(reshape_plants)


#### XY coordinates: ####
# To build spatial explicit JSDM, must include a matrix of xy coordinates for additional random effect:
xycoords <- xycoords %>% 
  dplyr::select(unq_plot, easting, northing) %>% 
  dplyr::filter(unq_plot %in% env$unq_plot) 

rownames(xycoords) <- xycoords$unq_plot

xycoords <- xycoords %>% 
  dplyr::select(easting, northing)

xycoords <- as.matrix(xycoords)
#### Run HMSC ####
XData <- X

XFormula <- ~l.area.std * shoredist.std + dist_near.std + sqrt.sum.std + fs_pc1.std + sm_av.std + slope.std + avg.n.std + avg.d15n.std

studyDesign <- env %>% 
  dplyr::filter(unq_plot %in% plants$unq_plot) %>% 
  dplyr::select(unq_plot, unq_tran, island)

# Not sure why this doesn't work without converting to a character and then back to a factor after but otherwise I get an error saying: Error in checkForRemoteErrors(val) : 
# 2 nodes produced errors; first error: subscript out of bounds
studyDesign$unq_plot <- as.character(studyDesign$unq_plot)
studyDesign$unq_tran <- as.character(studyDesign$unq_tran)
studyDesign$island <- as.character(studyDesign$island)
rownames(studyDesign) <- NULL

studyDesign <- studyDesign %>% 
  dplyr::select(island, unq_tran, unq_plot) 

studyDesign$unq_plot <- as.factor(studyDesign$unq_plot)
studyDesign$unq_tran <- as.factor(studyDesign$unq_tran)
studyDesign$island <- as.factor(studyDesign$island)
studyDesign$spatial <- as.factor(studyDesign$unq_plot)

n.plots <- studyDesign %>% 
  dplyr::select(unq_plot) %>% 
  dplyr::summarize(n = n());n.plots # 1381 plots

n.transects <- studyDesign %>% 
  dplyr::select(unq_tran) %>% 
  distinct();n.transects # 347 transects

n.islands <- studyDesign %>% 
  dplyr::select(island) %>% 
  distinct();n.islands # 90 islands

rL1 = HmscRandomLevel(units = unique(env$unq_plot))
rL2 = HmscRandomLevel(units = unique(env$unq_tran))
rL3 = HmscRandomLevel(units = unique(env$island))
# rL.spatial = HmscRandomLevel(sData = xycoords) # This option includes all all.
rL.spatial = HmscRandomLevel(sData = xycoords, sMethod = 'NNGP', nNeighbours = 10)
rL.spatial = setPriors(rL.spatial, nfMin = 1, nfMax = 1)

mod <- Hmsc(Y = Y,
            XData = X, 
            XFormula = XFormula,
            studyDesign = studyDesign,
            ranLevels = list(unq_plot = rL1, 
                             unq_tran = rL2, 
                             island = rL3,
                             spatial = rL.spatial),
            distr = "normal")

nChains <- 2

test.run = TRUE
if (test.run){
  # with this option, the vignette runs fast but results are not reliable 
  thin = 5
  nParallel = 3
  samples = 10
  transient = 0
  verbose = 1
} else {
  # with this option, the vignette evaluates slow but it reproduces the results of 
  # the .pdf version
  thin = 10
  samples = 25000
  transient = 1000
  verbose = 100
}

mod <- sampleMcmc(mod,
                  thin = thin,
                  samples = samples,
                  transient = transient,
                  nChains = nChains,
                  verbose = verbose,
                  updater = list(GammaEta = FALSE))

#### Step 2: Evaluate MCMC convergence: ####
m <- mod
mpost <- convertToCodaObject(m)
summary(mpost$Beta)

# From the HMSC 3.0 Univariate vignette: 
# Here, we want the two chains to yield approximately the same results. Second, they should mix well (go up and down without apparent autocorrelation). Finally, look to see that they seem to have reached a stationary distribution, as e.g. the first half of the recorded iterations looks statistically identical to the second half.
plot(mpost$Beta)

# More quantitatively, we can look at the following: 
effectiveSize(mpost$Beta) 

# If the effective sample sizes are close to theoretical value of the number of samples, it means there is little autocorrelation among consecutive samples.

# Evaluate potential scale reductor factors (psrfs)
gelman.diag(mpost$Beta,multivariate = FALSE)$psrf

# If potential scale reductor fators give numbers close to 1, it means the two chains give consistent results.

# Look at effective sample size and PSRF of beta estimates (environmental matrix)
ess.beta = effectiveSize(mpost$Beta)
psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf

hist(ess.beta, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (beta)", main = NULL)
hist(psrf.beta, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (beta)", main = NULL)

# Omega estimates (species interactions)
sppairs = matrix(sample(x = 1:18^2, size = 100)) # 25 species ^2
tmp = mpost$Omega[[1]]
for (chain in 1:length(tmp)){
  tmp[[chain]] = tmp[[chain]][,sppairs]
}
ess.omega = effectiveSize(tmp)
psrf.omega = gelman.diag(tmp, multivariate = FALSE)$psrf

hist(ess.omega, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (omega)", main = NULL)
hist(psrf.omega, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (gamma)", main = NULL)

par(mfrow = c(2,2))

hist(ess.beta, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (beta)", main = NULL)
hist(psrf.beta, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (beta)", main = NULL)

hist(ess.omega, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "ESS (omega)", main = NULL)
hist(psrf.omega, breaks = 30, cex.lab = 1.5, cex.axis = 1.5, xlab = "PSRF (omega)", main = NULL)

#### Step 3: Evaluate model fit (explanatory power): ####
# Evaluate the explanatory power of the model: 
# Computes predicted values by taking the posterior medians.
preds <- computePredictedValues(m, expected = FALSE)

MF <- evaluateModelFit(hM = m, predY = preds)

# From the vignette: 
# RMSE = Root-mean-square-error between the predicted and observed values times the sign of the corelation. 
# SR2 = pseudo-R2 for poisson models. This is computed as the squared spearman correlation between observed and predicted values times the sign of the correlation.
# O.RMSE, O.AUC, O.Tjur2 = observed and predicted data are also truncated to occurrences (presence/absences). 
# C.RMSE and C.SR2 are conditional on presence.

par(mfrow = c(1,2))

hist(MF$RMSE, 
     main = paste0("Mean = ", round(mean(MF$RMSE), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$R2, main = paste0("Mean = ", round(mean(MF$R2), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$O.AUC, main = paste0("Mean = ", round(mean(MF$O.AUC), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$O.TjurR2, main = paste0("Mean = ", round(mean(MF$O.TjurR2), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$O.RMSE, main = paste0("Mean = ", round(mean(MF$O.RMSE), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$C.SR2, main = paste0("Mean = ", round(mean(MF$C.SR2), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)
hist(MF$C.RMSE, main = paste0("Mean = ", round(mean(MF$C.RMSE), 2)),
     breaks = 10, 
     cex.lab = 1.5, 
     cex.axis = 1.5)

#### Results: ####
postBeta = getPostEstimate(m, parName = "Beta")
mbeta = postBeta$mean
betaP = postBeta$support
supportLevel = 0.95

# Means: 
# Mean coefficients: 
toPlot = mbeta
toPlot = toPlot * ((betaP > supportLevel) + (betaP < (1 - supportLevel)) > 0)
betaMat = matrix(toPlot, nrow = m$nc, ncol = ncol(m$Y))

colnames(betaMat) <- colnames(Y)

rownames(betaMat) <- c("Intercept", "Island area", "Distance to land", "Plot slope", "Soil moisture", "Forest openness", "Wrack biomass", "Distance to shore", "Average soil %N", "Average soil d15N", "Island area:Distance to shore")

rownames(betaMat)[rownames(betaMat) == "Average soil d15N"] <- paste("Average soil", paste(greeks("delta"), "15N", sep = ""))
rownames(betaMat)[rownames(betaMat) == "Average soil d15N"] <- paste("Average soil", paste(":delta", "15N", sep = ""))


# Sort the matrix by row of interest:
betaMat2 <- betaMat[,order(betaMat[10,]), drop=F]

# Drop the intercept:
betaMat3 <- betaMat2[-1, ]


# Plot:
corrplot(betaMat3,
         method = "color",
         mar = c(0, 0, 0, 0),
         tl.col = "black",
         tl.cex = 0.75,
         tl.srt = 45,
         cl.pos = "r",
         cl.ratio = 0.2,
         cl.length = 3,
         col.lim = c(-0.11, 0.11),
         font = 3,
         yaxt = "none",
         bg = "white",
         addgrid.col = "grey",
         is.corr = FALSE)

# Extract species-to-species associations: 
# Plot-level:
OmegaCor = computeAssociations(m) # This converts covariances to the scale of correlations (-1 to 1)
supportLevel = 0

toPlot = ((OmegaCor[[1]]$support > supportLevel) 
          + (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) * OmegaCor[[1]]$mean
par(mfrow = c(1,1))

corrplot(toPlot, 
         method = "color",
         title = "Plot-level species interactions", 
         mar=c(0,0,1,0),
         order = "FPC",
         type = "lower",
         bg = "white", 
         tl.col = "black", 
         tl.cex = 1,
         tl.srt = 45,
         font = 3)

# Island-level (transect level is 2):
toPlot3 = ((OmegaCor[[3]]$support > supportLevel) 
           + (OmegaCor[[3]]$support < (1 - supportLevel)) > 0) * OmegaCor[[3]]$mean

par(mfrow = c(1,1))

corrplot(toPlot3, 
         method = "color",
         title = "Island-level co-occurrences", 
         mar=c(0,0,1,0),
         order = "FPC",
         type = "lower",
         bg = "white", 
         tl.col = "black", 
         tl.cex = 1,
         tl.srt = 45,
         font = 3)

# Multipanel plot:
par(mfrow = c(1, 2))

corrplot(toPlot, 
         method = "color",
        # title = "Plot-level co-occurrences",
         mar=c(0,0,0,0),
         order = "FPC",
         type = "lower",
         bg = "white", 
         tl.col = "black", 
         tl.cex = 0.65,
         tl.srt = 45,
         cl.cex = 0.8,
         cl.length = 3,
         font = 3)

mtext("a", side = 3, at = -4, cex = 1.25)
mtext("Plot-level co-occurrences", at = 9.5, padj = 25)

corrplot(toPlot3, 
         method = "color",
        # title = "Island-level co-occurrences", 
         mar=c(0,0,0,0),
         order = "FPC",
         type = "lower",
         bg = "white", 
         tl.col = "black", 
         tl.cex = 0.65,
         tl.srt = 45,
         cl.cex = 0.8,
         cl.length = 3,
         font = 3)
mtext("b", side = 3, at = -4, cex = 1.25)
mtext("Island-level co-occurrences", at = 10, padj = 25)

dev.off()

# Multipanel: species interactions plot: 
layout.matrix <- matrix(c(1, 2, 0, 3), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(1, 1)) # Widths of the two columns

layout.show(3)

layout.matrix <- matrix(c(1, 2, 1, 3), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(1.8, 2), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns

# par(mfrow = c(3, 1))

corrplot(betaMat3,
         method = "color",
         mar = c(0, 0, 0, 0),
         tl.col = "black",
         tl.cex = 1,
         tl.srt = 45,
         cl.pos = "r",
         cl.ratio = 0.2,
         cl.length = 3,
         col.lim = c(-0.11, 0.11),
         font = 3,
         yaxt = "none",
         bg = "white",
         addgrid.col = "grey",
         is.corr = FALSE)
mtext("a", side = 3, at = -4, cex = 1.25) 

  
corrplot(toPlot, 
         method = "color",
         # title = "Plot-level co-occurrences",
         mar=c(0,0,0,0),
         order = "FPC",
         type = "lower",
         bg = "white", 
         tl.col = "black", 
         tl.cex = 1,
         tl.srt = 45,
         cl.cex = 1,
         cl.length = 3,
         font = 3) 

mtext("b", side = 3, at = -4, cex = 1.25)
mtext("Plot-level co-occurrences", at = 9.5, padj = 28)

corrplot(toPlot3, 
         method = "color",
         # title = "Island-level co-occurrences", 
         mar=c(0,0,0,0),
         order = "FPC",
         type = "lower",
         bg = "white", 
         tl.col = "black", 
         tl.cex = 1,
         tl.srt = 45,
         cl.cex = 1,
         cl.length = 3,
         font = 3)
mtext("c", side = 3, at = -4, cex = 1.25)
mtext("Island-level co-occurrences", at = 10, padj = 28)

#### Variance partitioning: ####
VP = computeVariancePartitioning(m, 
                                 group = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3), 
                                 groupnames = c("Island biogeography","Local environment", "Marine subsidies"))

library(PNWColors)
pnw_colors <- pnw_palette(name = "Lake", n = 7)

plotVariancePartitioning(m, 
                         VP = VP,
                         cols = pnw_colors, 
                         cex.names = 0.8, 
                         las = 2,
                         args.legend = list(x = 29, 
                                            y = 0.25, 
                                            cex = 0.7))
