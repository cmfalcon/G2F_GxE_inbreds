##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-linear regression
## Date: 2018-09-21
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(broman)

rowMedian <- function(x, na.rm = FALSE) apply(x, 1, median, na.rm = na.rm)
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

#### Read in, format, and organize data ####
selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")


#### Linear regression models and plots (slope and MSE) ####

selectedPhenoData_standardized <- selectedPhenoData
selectedPhenoData_standardized <- selectedPhenoData_standardized %>% transmute(
  Anthesis_st = (`Anthesis (GDU)` - mean(`Anthesis (GDU)`, na.rm = T)) / sd(`Anthesis (GDU)`, na.rm = T),
  Silking_st = (`Silking (GDU)` - mean(`Silking (GDU)`, na.rm = T)) / sd(`Silking (GDU)`, na.rm = T),
  PlantHeight_st = (`Plant height` - mean(`Plant height`, na.rm = T)) / sd(`Plant height`, na.rm = T),
  EarHeight_st = (`Ear height` - mean(`Ear height`, na.rm = T)) / sd(`Ear height`, na.rm = T),
  PlotWeight_st = (`Plot grain weight` - mean(`Plot grain weight`, na.rm = T)) / sd(`Plot grain weight`, na.rm = T),
  EarLength_st = (`Ear length` - mean(`Ear length`, na.rm = T)) / sd(`Ear length`, na.rm = T),
  EarWidth_st = (`Ear width` - mean(`Ear width`, na.rm = T)) / sd(`Ear width`, na.rm = T),
  KernelsPerRow_st = (`Kernels per row` - mean(`Kernels per row`, na.rm = T)) / sd(`Kernels per row`, na.rm = T),
  KRN_st = (`Kernel row number` - mean(`Kernel row number`, na.rm = T)) / sd(`Kernel row number`, na.rm = T),
  KernelWeight_st = (`Kernel weight` - mean(`Kernel weight`, na.rm = T)) / sd(`Kernel weight`, na.rm = T),
  KernelArea_st = (`Kernel area` - mean(`Kernel area`, na.rm = T)) / sd(`Kernel area`, na.rm = T),
  KernelLength_st = (`Kernel length` - mean(`Kernel length`, na.rm = T)) / sd(`Kernel length`, na.rm = T),
  KernelWidth_st = (`Kernel width` - mean(`Kernel width`, na.rm = T)) / sd(`Kernel width`, na.rm = T),
  KernelThickness_st = (`Kernel thickness` - mean(`Kernel thickness`, na.rm = T)) / sd(`Kernel thickness`, na.rm = T)
)  

selectedPhenoData_standardized <- cbind(selectedPhenoData[,1:4], selectedPhenoData_standardized)

#### Calculate environment means
for (i in 5:ncol(selectedPhenoData_standardized)) {
  # i = 5
  means_by_exp <- aggregate(selectedPhenoData_standardized,
                            by = list(selectedPhenoData_standardized$Enviro), FUN = mean, na.rm = T)
  means_by_exp$Enviro <- means_by_exp$Group.1
  means_by_exp <- means_by_exp[,-1]
}

####

names <- colnames(means_by_exp[ ,5:ncol(means_by_exp)])

index <- 1

inbreds <- unique(selectedPhenoData_standardized$Genotype)
traits <- colnames(selectedPhenoData_standardized[ ,5:ncol(selectedPhenoData_standardized)]) 

reg <- matrix(nrow = length(inbreds), ncol = 6, 
              dimnames = list(inbreds, c("beta0", "beta1", "MSE", "nobs", "nenv", "var_bw_locs")))
allReg <- matrix(nrow = length(inbreds), ncol = 0)
fred2 <- character(0)

slopes <- matrix(nrow = length(inbreds), ncol = 1, dimnames = list(inbreds, "slope"))
MSE <- matrix(nrow = length(inbreds), ncol = 1, dimnames = list(inbreds, "MSE"))
Rsquared_of_model <- matrix(nrow = length(inbreds), ncol = 1, dimnames = list(inbreds, "R^2"))

allSlopes <- matrix(nrow = (length(inbreds)), ncol = 0)
allMSE <- matrix(nrow = length(inbreds), ncol = 0)
all_Rsquared_of_model <- matrix(nrow = length(inbreds), ncol = 0)

for (trait in traits) {
  # trait = "Anthesis_st"
  dir.create(paste0("./standardized", trait))
  traitName <- names[index]
  
  expMeans <- means_by_exp[order(means_by_exp[ ,trait]), ]
  expMeans$rank <- 1:nrow(expMeans)
  expMeans <- expMeans[complete.cases(expMeans[,trait]), ]
  
  for (inbred in inbreds) {
    # inbred = "2369"
    inbredDat <- selectedPhenoData_standardized[selectedPhenoData_standardized$Genotype == inbred, c("Enviro", trait)]
    
    if (all(is.na(inbredDat[,trait]))) {
      reg[inbred, ] <- NA
      next
    }
    inbredDat <- inbredDat[!is.na(inbredDat[,trait]), ]
    
    nenv <- length(unique(inbredDat$Enviro))
    
    inb_loc_means <- aggregate(inbredDat, by = list(inbredDat$Enviro), FUN = mean, na.rm = T)
    inb_mean <- mean(inbredDat[,trait])
    inb_var <- var(inb_loc_means[,trait], na.rm = T)
    
    regression_dat <- merge(expMeans[,c("Enviro", trait)],
                            inbredDat,
                            by = "Enviro",
                            all = T)
    if (all(is.na(regression_dat[,3]))) {
      reg[inbred, "beta0"] <- NA
      reg[inbred, "beta1"] <- NA
      reg[inbred, "MSE"] <- NA
    } else {
      trait_regression <- lm(regression_dat[,3] ~ regression_dat[,2])
      reg[inbred, "beta0"] <- trait_regression$coefficients[1]
      reg[inbred, "beta1"] <- trait_regression$coefficients[2]
      reg[inbred, "MSE"] <- sum(trait_regression$residuals^2) / trait_regression$df.residual
      slopes[inbred, "slope"] <- trait_regression$coefficients[2]
      MSE[inbred, "MSE"] <- sum(trait_regression$residuals^2) / trait_regression$df.residual
      Rsquared_of_model[inbred,"R^2"] <- summary(trait_regression)$r.squared
      }
    
    reg[inbred, "var_bw_locs"] <- inb_var
    reg[inbred, "nobs"] <- nrow(inbredDat)
    reg[inbred, "nenv"] <- nenv
  }
  
  write.csv(reg, paste0("./standardized", trait, "/regression_", trait, ".csv"), sep = "\t", col.names = NA, row.names = T)
  
  allReg <- cbind(allReg, reg)
  fred <- c(paste0("beta0_", trait), paste0("Slope - ", trait),
            paste0("MSE_", trait), paste0("nobs_", trait),
            paste0("nenv_", trait), paste0("var_bw_locs_", trait))
  fred2 <- c(fred2,fred)
  
  allSlopes <- cbind(allSlopes, slopes)
  allMSE <- cbind(allMSE, MSE)
  all_Rsquared_of_model <- cbind(all_Rsquared_of_model, Rsquared_of_model)
  
  pdf(paste0("./standardized", trait, "/stability_histograms_",trait,".pdf"), width = 6, height = 6)
  par(mfrow = c(3,2))
  hist(reg[,1], breaks = 20, col = "cadetblue", main = "Intercept", ylab = "Count", xlab = "")
  hist(reg[,2], breaks = 20, col = "cadetblue", main = "Slope", ylab = "Count", xlab = "")
  hist(reg[,3], breaks = 20, col = "cadetblue", main = "MSE", ylab = "Count", xlab = "")
  hist(reg[,6], breaks = 20, col = "cadetblue", main = "Variance Between Enviros", ylab = "Count", xlab = "")
  hist(reg[,4], breaks = 20, col = "cadetblue", main = "Number of Obs", ylab = "Count", xlab = "")
  hist(reg[,5], breaks = 20, col = "cadetblue", main = "Number of Environments", ylab = "Count", xlab = "")
  dev.off()
  
  pdf(paste0("./standardized", trait, "/stability_regressions_", trait, ".pdf"), width = 12, height = 6)
  spread <- range(selectedPhenoData_standardized[,trait], na.rm = T)
  plot(1, xlim = c((min(spread)*.9),(max(spread)*1.1)), ylim = c((min(spread)*.9), (max(spread)*1.1)),
       main = "Inbred Regression \nAcross Environment Means",
       xlab = paste0("Environmental ", traitName, " Mean"),
       ylab = paste0(traitName, " by Inbred"))
  for (i in 1:nrow(reg)) {
    if (!is.na(reg[i,2])) {
      abline(a = reg[i,1], b = reg[i,2])
    }}
  abline(a = 0, b = 1, col = "red")
  dev.off()
  
  
  index <- index + 1
}

colnames(allReg) <- c(fred2)
colnames(allSlopes) <- names
colnames(allMSE) <- names
colnames(all_Rsquared_of_model) <- names

write.csv(allReg, "selectedTraitsRegression_traitsStandardized.csv")
write.csv(allSlopes, "selectedTraitsSlopes_traitsStandardized.csv")
write.csv(allMSE, "selectedTraitsMSE_traitsStandardized.csv")
write.csv(all_Rsquared_of_model, "selectedTraits_Rsquaredofmodel_traitsStandardized.csv")

### Create slope summary
colMins <- as.data.frame(apply(X = allSlopes[,1:ncol(allSlopes)], MARGIN = 2, min, na.rm = T))
colMaxs <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, max, na.rm = T))
colRanges <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, range, na.rm = T))
colMeans <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, mean, na.rm = T))
colMedians <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, median, na.rm = T))
colVars <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, var, na.rm = T))
traits <- as.data.frame(rownames(colMins))

slopeSummary <- cbind(traits,colMins,colMedians,colMaxs, colVars)
colnames(slopeSummary) <- c("Trait", "Min", "Median", "Max", "Var")
rownames(slopeSummary) <- seq(1:nrow(slopeSummary))

write.csv(slopeSummary, "slope_summary.csv", row.names = F)

rm(colMaxs, colMeans, colMedians, colMins, colRanges, colVars)

allSlopesWithMedian <- as.data.frame(allSlopes)
allSlopesWithMedian$Anthesis_st <- as.numeric(myround(allSlopesWithMedian$Anthesis_st, 2))
allSlopesWithMedian$Silking_st <- as.numeric(myround(allSlopesWithMedian$Silking_st, 2))
allSlopesWithMedian$PlantHeight_st <- as.numeric(myround(allSlopesWithMedian$PlantHeight_st, 2))
allSlopesWithMedian$EarHeight_st <- as.numeric(myround(allSlopesWithMedian$EarHeight_st, 2))
allSlopesWithMedian$PlotWeight_st <- as.numeric(myround(allSlopesWithMedian$PlotWeight_st, 2))
# allSlopesWithMedian$CupWeight_st <- as.numeric(myround(allSlopesWithMedian$CupWeight_st, 2))
allSlopesWithMedian$KernelWidth_st <- as.numeric(myround(allSlopesWithMedian$KernelWidth_st, 2))
allSlopesWithMedian$KernelLength_st <- as.numeric(myround(allSlopesWithMedian$KernelLength_st, 2))
allSlopesWithMedian$KernelArea_st <- as.numeric(myround(allSlopesWithMedian$KernelArea_st, 2))
allSlopesWithMedian$KernelWeight_st <- as.numeric(myround(allSlopesWithMedian$KernelWeight_st, 2))
allSlopesWithMedian$KRN_st <- as.numeric(myround(allSlopesWithMedian$KRN_st, 2))
allSlopesWithMedian$KernelThickness_st <- as.numeric(myround(allSlopesWithMedian$KernelThickness_st, 2))
allSlopesWithMedian$KernelsPerRow_st <- as.numeric(myround(allSlopesWithMedian$KernelsPerRow_st, 2))
allSlopesWithMedian$EarLength_st <- as.numeric(myround(allSlopesWithMedian$EarLength_st, 2))
allSlopesWithMedian$EarWidth_st <- as.numeric(myround(allSlopesWithMedian$EarWidth_st, 2))
allSlopesWithMedian$medianSlope <- as.numeric(myround(rowMedian(allSlopesWithMedian[,1:ncol(allSlopesWithMedian)], na.rm = T), 2))

write.csv(allSlopesWithMedian, "all_slopes_with_median.csv")

### Create MSE summary
colMins <- as.data.frame(apply(X = allMSE[,1:ncol(allMSE)], MARGIN = 2, min, na.rm = T))
colMaxs <-  as.data.frame(apply(allMSE[,1:ncol(allMSE)], 2, max, na.rm = T))
colRanges <-  as.data.frame(apply(allMSE[,1:ncol(allMSE)], 2, range, na.rm = T))
colMeans <-  as.data.frame(apply(allMSE[,1:ncol(allMSE)], 2, mean, na.rm = T))
colMedians <-  as.data.frame(apply(allMSE[,1:ncol(allMSE)], 2, median, na.rm = T))
colVars <- as.data.frame(apply(allMSE[,1:ncol(allMSE)], 2, var, na.rm = T))
traits <- as.data.frame(rownames(colMins))

MSEsummary <- cbind(traits,colMins,colMedians,colMaxs,colVars)
colnames(MSEsummary) <- c("Trait", "Min", "Median", "Max", "Var")
rownames(MSEsummary) <- seq(1:nrow(slopeSummary))

write.csv(MSEsummary, "MSE_summary.csv", row.names = F)

rm(colMaxs, colMeans, colMins, colRanges, colVars)

allMSEWithMedian <- as.data.frame(allMSE)
allMSEWithMedian$Anthesis_st <- as.numeric(myround(allMSEWithMedian$Anthesis_st, 2))
allMSEWithMedian$Silking_st <- as.numeric(myround(allMSEWithMedian$Silking_st, 2))
allMSEWithMedian$PlantHeight_st <- as.numeric(myround(allMSEWithMedian$PlantHeight_st, 2))
allMSEWithMedian$EarHeight_st <- as.numeric(myround(allMSEWithMedian$EarHeight_st, 2))
allMSEWithMedian$PlotWeight_st <- as.numeric(myround(allMSEWithMedian$PlotWeight_st, 2))
# allMSEWithMedian$CupWeight_st <- as.numeric(myround(allMSEWithMedian$CupWeight_st, 2))
allMSEWithMedian$KernelWidth_st <- as.numeric(myround(allMSEWithMedian$KernelWidth_st, 2))
allMSEWithMedian$KernelLength_st <- as.numeric(myround(allMSEWithMedian$KernelLength_st, 2))
allMSEWithMedian$KernelArea_st <- as.numeric(myround(allMSEWithMedian$KernelArea_st, 2))
allMSEWithMedian$KernelWeight_st <- as.numeric(myround(allMSEWithMedian$KernelWeight_st, 2))
allMSEWithMedian$KRN_st <- as.numeric(myround(allMSEWithMedian$KRN_st, 2))
allMSEWithMedian$KernelThickness_st <- as.numeric(myround(allMSEWithMedian$KernelThickness_st, 2))
allMSEWithMedian$KernelsPerRow_st <- as.numeric(myround(allMSEWithMedian$KernelsPerRow_st, 2))
allMSEWithMedian$EarLength_st <- as.numeric(myround(allMSEWithMedian$EarLength_st, 2))
allMSEWithMedian$EarWidth_st <- as.numeric(myround(allMSEWithMedian$EarWidth_st, 2))
allMSEWithMedian$medianMSE <- as.numeric(myround(rowMedian(allMSEWithMedian[,1:ncol(allMSEWithMedian)], na.rm = T), 2))

write.csv(allMSEWithMedian, "all_MSE_with_median.csv")

#### Plot slope and MSE for traits and genotypes ####

### Melt slope and MSE matrices for making plots
colnames(allSlopes) <- c("Anthesis", "Silking", 
                         "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                         "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                         "Kernel length", "Kernel width", "Kernel thickness")
meltedSlope <- melt(allSlopes, id = "Slope")
meltedSlope$value <- as.numeric(meltedSlope$value)
colnames(meltedSlope) <- c("Genotype", "Trait", "Slope")

colnames(allMSE) <- c("Anthesis", "Silking", 
                      "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                      "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                      "Kernel length", "Kernel width", "Kernel thickness")
meltedMSE <- melt(allMSE, id = "MSE")
meltedMSE$value <- as.numeric(meltedMSE$value)
colnames(meltedMSE) <- c("Genotype", "Trait", "StandMSE")
meltedMSE$Trait <- gsub("_st", "", meltedMSE$Trait)

### Make boxplots
slopeBoxplot <-
  ggplot(meltedSlope, aes(x = reorder(Trait, X = Slope, FUN = median), y = Slope)) +
  geom_jitter(color = "gray40", shape = 1, width = 0.15) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  labs(x = "", y = "Slope") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 70, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')
ggsave(slopeBoxplot, filename = "boxplot_slopeByTrait_gxe20142015Inbreds.tif", device = "tiff",
       height = 7, width = 7, units = "in")

MSEBoxplot <-
  ggplot(meltedMSE, aes(x = reorder(Trait, X = StandMSE, FUN = median), y = StandMSE)) +
  geom_jitter(color = "gray40", shape = 1, width = 0.15) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  labs(x = "Trait", y = "MSE") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 70, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none') 
ggsave(MSEBoxplot, filename = "boxplot_standMSEByTrait_gxe20142015Inbreds.tif", device = "tiff",
       height = 7, width = 7, units = "in")

MSEBoxplot_limit_y_axis <-
  ggplot(meltedMSE, aes(x = reorder(Trait, X = StandMSE, FUN = median), y = StandMSE)) +
  geom_jitter(color = "gray40", shape = 1, width = 0.15) +
  geom_boxplot(fill = NA, outlier.shape = NA) +
  labs(x = "Trait", y = "MSE") +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.25)) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 70, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none') 
ggsave(MSEBoxplot_limit_y_axis, filename = "boxplot_standMSEByTrait_limit_y_axis_gxe20142015Inbreds.tif", device = "tiff",
       height = 7, width = 7, units = "in")
combinedSlopeMSEboxplots <- ggarrange(slopeBoxplot, MSEBoxplot_limit_y_axis, 
                                      labels = c("A", "B"),
                                      ncol = 1, nrow = 2)

ggsave(combinedSlopeMSEboxplots, filename = "combined_standardized_slope_and_MSE_boxplots.tif", device = "tiff",
       height = 11, width = 8.5, units = "in")


