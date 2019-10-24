##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-Stability analysis in 
##                 Midwest locations
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(dplyr)
library(plyr)
library(broman)
library(GGEBiplots)
library(ggpubr)
library(data.table)
library(forcats)
library(corrplot)
library(ggcorrplot)
library(agricolae)

rowMedian <- function(x, na.rm = FALSE) apply(x, 1, median, na.rm = na.rm)
mean_ <- function(...) mean(..., na.rm = T)
median_ <- function(...) median(..., na.rm = T)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#### Read in, format, and organize data ####

selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")

allBlups <- as.data.frame(read.csv("./gxe20142015_allTraits_BLUPs.csv", header = T))
colnames(allBlups) <- c("Genotype", "Anthesis", "Silking", 
                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                        "Kernel length", "Kernel width", "Kernel thickness")

allSlopesWithMedian <- read.csv("all_slopes_with_median.csv", header = T)
colnames(allSlopesWithMedian) <- c("Genotype", "Anthesis", "Silking", 
                                   "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                   "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                   "Kernel length", "Kernel width", "Kernel thickness", "Median")

allMSEWithMedian <- read.csv("all_MSE_with_median.csv", header = T)
colnames(allMSEWithMedian) <- c("Genotype", "Anthesis", "Silking", 
                                "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                "Kernel length", "Kernel width", "Kernel thickness", "Median")

stabilityValue_allTraits <- read.csv("stability_value_all_traits.csv", header = T)
colnames(stabilityValue_allTraits) <- c("Genotype", "Anthesis", "Silking", 
                                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                        "Kernel length", "Kernel width", "Kernel thickness")
stabilityValue_allTraits_melt <- melt(stabilityValue_allTraits)

stabilityValue_allTraits_withMed <- read.csv("stability_value_all_traits_with_median.csv", header = T)
colnames(stabilityValue_allTraits_withMed) <- c("Genotype", "Anthesis", "Silking", 
                                                "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                                "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                                "Kernel length", "Kernel width", "Kernel thickness", "Median")

#### Stability analysis in Midwest locations ####
# dir.create("./stability_in_midwest")
setwd("./stability_in_midwest/")
midwest_pheno_data <- selectedPhenoData %>%
  filter(Enviro == "IA1_14" |
           Enviro == "IA1_15" |
           Enviro == "IA2_14" |
           Enviro == "IA2_15" |
           Enviro == "IA3_14" |
           Enviro == "IA3_15" |
           Enviro == "IA4_15" |
           Enviro == "IL1_14" |
           Enviro == "IL1_15" |
           Enviro == "IN1_14" |
           Enviro == "IN1_15" |
           Enviro == "MN1_15" |
           Enviro == "MN2_14" |
           Enviro == "MO1_15" |
           Enviro == "MO2_15" |
           Enviro == "WI1_14" |
           Enviro == "WI1_15" |
           Enviro == "WI2_15")

#### Linear regression models and plots (slope and MSE) --- TRAITS STANDARDIZED ####
midwest_pheno_data_standardized <- midwest_pheno_data
midwest_pheno_data_standardized <- midwest_pheno_data_standardized %>% transmute(
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

midwest_pheno_data_standardized <- cbind(midwest_pheno_data[,1:4], midwest_pheno_data_standardized)


#### Calculate environment means
for (i in 5:ncol(midwest_pheno_data_standardized)) {
  # i = 5
  means_by_exp <- aggregate(midwest_pheno_data_standardized,
                            by = list(midwest_pheno_data_standardized$Enviro), FUN = mean, na.rm = T)
  means_by_exp$Enviro <- means_by_exp$Group.1
  means_by_exp <- means_by_exp[,-1]
}

####

names <- colnames(means_by_exp[ ,5:ncol(means_by_exp)])

index <- 1

inbreds <- unique(midwest_pheno_data_standardized$Genotype)
traits <- colnames(midwest_pheno_data_standardized[ ,5:ncol(midwest_pheno_data_standardized)]) 

reg <- matrix(nrow = length(inbreds), ncol = 6, 
              dimnames = list(inbreds, c("beta0", "beta1", "MSE", "nobs", "nenv", "var_bw_locs")))
allReg <- matrix(nrow = length(inbreds), ncol = 0)
fred2 <- character(0)

slopes <- matrix(nrow = length(inbreds), ncol = 1, dimnames = list(inbreds, "slope"))
MSE <- matrix(nrow = length(inbreds), ncol = 1, dimnames = list(inbreds, "MSE"))

allSlopes <- matrix(nrow = (length(inbreds)), ncol = 0)
allMSE <- matrix(nrow = length(inbreds), ncol = 0)

for (trait in traits) {
  # trait = "Anthesis_st"
  dir.create(paste0("./standardized_", trait))
  traitName <- names[index]
  
  expMeans <- means_by_exp[order(means_by_exp[ ,trait]), ]
  expMeans$rank <- 1:nrow(expMeans)
  expMeans <- expMeans[complete.cases(expMeans[,trait]), ]
  
  for (inbred in inbreds) {
    # inbred = "2369"
    inbredDat <- midwest_pheno_data_standardized[midwest_pheno_data_standardized$Genotype == inbred, c("Enviro", trait)]
    
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
    }
    
    reg[inbred, "var_bw_locs"] <- inb_var
    reg[inbred, "nobs"] <- nrow(inbredDat)
    reg[inbred, "nenv"] <- nenv
  }
  
  write.csv(reg, paste0("./standardized_", trait, "/regression_", trait, ".csv"), sep = "\t", col.names = NA, row.names = T)
  
  allReg <- cbind(allReg, reg)
  fred <- c(paste0("beta0_", trait), paste0("Slope - ", trait),
            paste0("MSE_", trait), paste0("nobs_", trait),
            paste0("nenv_", trait), paste0("var_bw_locs_", trait))
  fred2 <- c(fred2,fred)
  
  allSlopes <- cbind(allSlopes, slopes)
  allMSE <- cbind(allMSE, MSE)
  
  pdf(paste0("./standardized_", trait, "/stability_histograms_",trait,".pdf"), width = 6, height = 6)
  par(mfrow = c(3,2))
  hist(reg[,1], breaks = 20, col = "cadetblue", main = "Intercept", ylab = "Count", xlab = "")
  hist(reg[,2], breaks = 20, col = "cadetblue", main = "Slope", ylab = "Count", xlab = "")
  hist(reg[,3], breaks = 20, col = "cadetblue", main = "MSE", ylab = "Count", xlab = "")
  hist(reg[,6], breaks = 20, col = "cadetblue", main = "Variance Between Enviros", ylab = "Count", xlab = "")
  hist(reg[,4], breaks = 20, col = "cadetblue", main = "Number of Obs", ylab = "Count", xlab = "")
  hist(reg[,5], breaks = 20, col = "cadetblue", main = "Number of Environments", ylab = "Count", xlab = "")
  dev.off()
  
  pdf(paste0("./standardized_", trait, "/stability_regressions_", trait, ".pdf"), width = 12, height = 6)
  spread <- range(midwest_pheno_data_standardized[,trait], na.rm = T)
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

write.csv(allReg, "selectedTraitsRegression_traitsStandardized.csv")
write.csv(allSlopes, "selectedTraitsSlopes_traitsStandardized.csv")
write.csv(allMSE, "selectedTraitsMSE_traitsStandardized.csv")

### Create slope summary
colMins <- as.data.frame(apply(X = allSlopes[,1:ncol(allSlopes)], MARGIN = 2, min, na.rm = T))
colMaxs <-  as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, max, na.rm = T))
colRanges <-  as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, range, na.rm = T))
colMeans <-  as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, mean, na.rm = T))
colMedians <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, median, na.rm = T))
colVars <- as.data.frame(apply(allSlopes[,1:ncol(allSlopes)], 2, var, na.rm = T))
traits <- as.data.frame(rownames(colMins))

slopeSummary <- cbind(traits,colMins,colMedians,colMaxs, colVars)
colnames(slopeSummary) <- c("Trait", "Min", "Median", "Max", "Var")
rownames(slopeSummary) <- seq(1:nrow(slopeSummary))

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

#### GGE biplot analysis ####
dir.create("./GGEbiplots_standardized")
dir.create("./boxplotsByEnviro_standardized_")

midwest_pheno_data_standardized$Plot = as.factor(midwest_pheno_data_standardized$Plot)
midwest_pheno_data_standardized$Enviro = as.factor(midwest_pheno_data_standardized$Enviro)
midwest_pheno_data_standardized$Genotype = as.factor(midwest_pheno_data_standardized$Genotype)
midwest_pheno_data_standardized$Rep = as.factor(midwest_pheno_data_standardized$Rep)

midwest_pheno_data_standardized$Enviro <- gsub("_20", "_", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("AZI", "AZ", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("DEI", "DE", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("GAI", "GA", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("IAI", "IA", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("ILI", "IL", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("INI", "IN", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("KSI", "KS", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("MNI", "MN", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("MOI", "MO", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("NCI", "NC", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("NYI", "NY", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("PAI", "PA", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("SDI", "SD", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("TXI", "TX", midwest_pheno_data_standardized$Enviro, fixed = T)
midwest_pheno_data_standardized$Enviro <- gsub("WII", "WI", midwest_pheno_data_standardized$Enviro, fixed = T)

#### create enviro-centered data ####
midwest_pheno_data_standardized$Plot <- as.factor(midwest_pheno_data_standardized$Plot)
midwest_pheno_data_standardized$Enviro <- as.factor(midwest_pheno_data_standardized$Enviro)
midwest_pheno_data_standardized$Rep <- as.factor(midwest_pheno_data_standardized$Rep)
midwest_pheno_data_standardized$Genotype <- as.factor(midwest_pheno_data_standardized$Genotype)

traitMean <- summarise_all(midwest_pheno_data_standardized, .funs = mean, na.rm = T)

enviroCenteredData <- midwest_pheno_data_standardized %>%
  group_by(Enviro) %>%
  mutate_all(funs(. - mean_(.)))

enviroCenteredData <- cbind(midwest_pheno_data_standardized[,1:4], enviroCenteredData[,5:18])

# dir.create("./GGEbiplots_standardized")

enviroCenteredData$Plot = as.factor(enviroCenteredData$Plot)
enviroCenteredData$Enviro = as.factor(enviroCenteredData$Enviro)
enviroCenteredData$Genotype = as.factor(enviroCenteredData$Genotype)
enviroCenteredData$Rep = as.factor(enviroCenteredData$Rep)

#### GGE model analysis ####
enviros <- sort(unique(enviroCenteredData$Enviro), decreasing = F)
genos <- sort(unique(enviroCenteredData$Genotype), decreasing = F)

discrimRank_allTraits <- matrix(enviros, nrow = length(unique(enviroCenteredData$Enviro)), ncol = 1)
stabilityRank_allTraits <- matrix(genos, nrow = length(unique(enviroCenteredData$Genotype)), ncol = 1)

colnames(discrimRank_allTraits) <- "enviros"
colnames(stabilityRank_allTraits) <- "genos"

GGE_PC_varExpl <- matrix(NA, nrow = 0, ncol = 4)
colnames(GGE_PC_varExpl) <- c("Trait", "PC1", "PC2", "PC1+PC2")

for (t in c(5:ncol(enviroCenteredData))) {
  # t = 6
  
  #### calculate means across reps and create a matrix of Gs by Es for each trait ####
  trait = enviroCenteredData[,t]
  df <- enviroCenteredData[!is.na(trait),]
  df <- droplevels(df)
  
  df2 <- daply(.data = df,
               .variables = c("Enviro", "Genotype"),
               .fun = function(df) mean(df[,t]))
  data <- as.data.frame(t(df2))
  
  #### remove enviros and then genos with >30% missing data ####
  dataRmHighNAEnviro <- data[,colMeans(is.na(data)) < 0.30]
  dataRmHighNAEnviroAndGeno <- dataRmHighNAEnviro[rowMeans(is.na(dataRmHighNAEnviro)) < 0.30,]
  
  data <- dataRmHighNAEnviroAndGeno 
  remove_NaN <- function(x) ifelse(is.nan(x), NA, x) 
  data <- apply(data, 2, remove_NaN)
  
  dataWithNAs <- data
  
  newData <- data
  
  #### estimate missing data as row/column mean ####
  for (r in 1:nrow(data)) {
    for (c in 1:ncol(data)) {
      # r = 9
      # c = 24
      if (is.na(data[r,c])) {
        geno_mean <- mean(data[r,], na.rm = T)
        enviro_mean <- mean(data[,c], na.rm = T)
        data[r,c] <- mean_(geno_mean, enviro_mean)
      }
    }
  }
  
  #### create models ####
  GGE1 <- GGEModel(data, centering = "tester", SVP = "column") # Model for which won where plot & discriminability vs representativeness plot
  GGE2 <- GGEModel(data, centering = "tester", SVP = "row") # Model for mean vs stability plot
  
  #### create which won where and discriminability vs representativenss plots ####
  whichWonWhere <- GGEPlot(GGE1, type = 6) # which won where
  ggsave(plot = whichWonWhere, filename = paste0("./GGEbiplots_standardized/WWW-", colnames(enviroCenteredData[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
  discrimVsRepresent <- GGEPlot(GGE1, type = 7, sizeGen = 2, sizeEnv = 3) # discrimination vs. representativeness
  ggsave(plot = discrimVsRepresent, filename = paste0("./GGEbiplots_standardized/DvR-", colnames(enviroCenteredData[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
  #### discriminability value output ####
  nenv <- ncol(data)
  ngen <- nrow(data)
  
  vectorLengths <- matrix(nrow = nenv, ncol = 1, dimnames = list(colnames(data), "vectorLength"))
  for (v in (ngen + 1):(ngen + nenv)) {
    # v = 65
    vectorLengths[(v - ngen),1] <- sqrt(discrimVsRepresent$data[v,1]^2 + (discrimVsRepresent$data[v,2])^2)
  }
  
  rankVectorLengths <- transform(vectorLengths, rank = ave(vectorLength,
                                                           FUN = function(x) rank(-x, ties.method = "min")))
  rankVectorLengths <- setDT(rankVectorLengths, keep.rownames = T)[]
  colnames(rankVectorLengths)[1] <- "enviros"
  colnames(rankVectorLengths)[2] <- paste0("vectorLength_", colnames(enviroCenteredData[t]))
  colnames(rankVectorLengths)[3] <- paste0("rank_", colnames(enviroCenteredData[t]))
  
  discrimRank_allTraits <- merge(discrimRank_allTraits, rankVectorLengths, by = "enviros", all = T)
  
  #### create mean vs stability plots ####
  meanVsStability <- GGEPlot(GGE2, type = 9, sizeEnv = 2, sizeGen = 3) # mean vs. stability
  ggsave(plot = meanVsStability, filename = paste0("./GGEbiplots_standardized/MvS-", colnames(enviroCenteredData[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
  #### stability value output ####
  avgEnvX <- mean(meanVsStability$data[(ngen + 1):(ngen + nenv),1])
  avgEnvY <- mean(meanVsStability$data[(ngen + 1):(ngen + nenv),2])
  
  m = abs(avgEnvY/avgEnvX)
  
  sinTheta <- m/(sqrt(1 + (m)^2))
  cosTheta <- 1/(sqrt(1 + (m)^2))
  
  genoLengths <- matrix(nrow = ngen, ncol = 1, dimnames = list(rownames(data), "genoLength"))
  for (g in 1:ngen) {
    # g = 1
    genoLengths[g, 1] <- abs((sinTheta * (meanVsStability$data[g,1])) + (cosTheta * (meanVsStability$data[g,2])))
  }
  
  rankGenoLengths <- transform(genoLengths, rank = ave(genoLengths,
                                                       FUN = function(x) rank(x, ties.method = "min")))
  rankGenoLengths <- setDT(rankGenoLengths, keep.rownames = T)[]
  colnames(rankGenoLengths)[1] <- "genos"
  colnames(rankGenoLengths)[2] <- paste0("genoLength_", colnames(enviroCenteredData[t]))
  colnames(rankGenoLengths)[3] <- paste0("rank_", colnames(enviroCenteredData[t]))
  
  stabilityRank_allTraits <- merge(stabilityRank_allTraits, rankGenoLengths, by = "genos", all = T)
  
  #### Output var explained by PC1, PC2, and PC1+PC2
  PC1_varExpl <- GGE1$varexpl[1]
  PC2_varExpl <- GGE1$varexpl[2]
  sumPC1PC2_varExpl <- sum(PC1_varExpl, PC2_varExpl)
  
  row <- c(colnames(enviroCenteredData)[t], PC1_varExpl, PC2_varExpl, sumPC1PC2_varExpl)
  
  GGE_PC_varExpl <- rbind(GGE_PC_varExpl, row)
}

write.csv(discrimRank_allTraits, "discrim_rank_all_traits.csv", row.names = F)
write.csv(stabilityRank_allTraits, "stability_rank_all_traits.csv", row.names = F)

#### Stability ####
stabilityRank_allTraits2 <- stabilityRank_allTraits[,seq(1, ncol(stabilityRank_allTraits), 2)]
# stabilityRank_allTraits2$meanRank <- as.numeric(myround(rowMeans(stabilityRank_allTraits2[,2:ncol(stabilityRank_allTraits2)], na.rm = T), 2))
stabilityRank_allTraits2$medianRank <- as.numeric(rowMedian(stabilityRank_allTraits2[,2:ncol(stabilityRank_allTraits2)], na.rm = T))
stabilityRank_allTraits2[stabilityRank_allTraits2 == "NaN"] <- NA_character_
write.csv(stabilityRank_allTraits2, "stabilityRank_allTraits2.csv", row.names = F)

stabilityRank_allTraits2_melt <- melt(stabilityRank_allTraits2)

stability <- ggplot(stabilityRank_allTraits2_melt, mapping = aes(x = reorder(genos, -value, median, na.rm = T), y = value)) +
  geom_boxplot() +
  # geom_jitter() +
  labs(x = "Genotype", y = "Stability rank") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = c(0.15,0.9), legend.title = element_text(" "))
ggsave(plot = stability, filename = "medianStability.tif", device = "tiff",
       width = 6, height = 7.75, units = "in")


stabilityValue_allTraits <- stabilityRank_allTraits[,c(1, seq(2, ncol(stabilityRank_allTraits), 2))]
# stabilityValue_allTraits$median <- as.numeric(rowMedian(stabilityValue_allTraits[,2:ncol(stabilityValue_allTraits)], na.rm = T))
stabilityValue_allTraits[stabilityValue_allTraits == "NaN"] <- NA_character_
traits <- colnames(selectedPhenoData[5:ncol(selectedPhenoData)])
colnames(stabilityValue_allTraits) <- c("Genotype", traits)
stabilityValue_allTraits <- subset(stabilityValue_allTraits, Genotype != "PHG83")
write.csv(stabilityValue_allTraits, "GGE_stability.csv", row.names = F)

stabilityValue_allTraits_withMed <- stabilityValue_allTraits
stabilityValue_allTraits_withMed$median <- as.numeric(rowMedian(stabilityValue_allTraits_withMed[,2:ncol(stabilityValue_allTraits_withMed)], na.rm = T))

stabilityValue_allTraits_melt <- melt(stabilityValue_allTraits)

stabilityValue <- ggplot(stabilityValue_allTraits_melt, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot(outlier.color = "white") +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_flip() +
  labs(x = "Genotypes", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue, filename = "medianStabilityValue.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")

stabilityValue2 <- ggplot(stabilityValue_allTraits_melt, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  # scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  # scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue2, filename = "medianStabilityValue_noJitter.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")


stabilityValue_allTraits_melt_traitType <- stabilityValue_allTraits_melt
stabilityValue_allTraits_melt_traitType$traitType <- NA
for (i in 1:nrow(stabilityValue_allTraits_melt_traitType)) {
  if (stabilityValue_allTraits_melt_traitType[i, "variable"] == "Anthesis (GDU)" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Silking (GDU)" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Plant height" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Ear height") {
    stabilityValue_allTraits_melt_traitType[i, "traitType"] <- "Non"
  } else {
    stabilityValue_allTraits_melt_traitType[i, "traitType"] <- "Yield-component"
  }
}

stabilityValue_allTraits_melt_traitType$traitType2 <- NA
for (i in 1:nrow(stabilityValue_allTraits_melt_traitType)) {
  if (stabilityValue_allTraits_melt_traitType[i, "variable"] == "Anthesis (GDU)" |
      stabilityValue_allTraits_melt_traitType[i, "variable"] == "Silking (GDU)") {
    stabilityValue_allTraits_melt_traitType[i, "traitType2"] <- "Flowering date"
  } else if (stabilityValue_allTraits_melt_traitType[i, "variable"] == "Plant height" |
             stabilityValue_allTraits_melt_traitType[i, "variable"] == "Ear height") {
    stabilityValue_allTraits_melt_traitType[i, "traitType2"] <- "Height"
  } else {
    stabilityValue_allTraits_melt_traitType[i, "traitType2"] <- "Yield-component"
  }
}

stabilityValue_traitType <- ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = traitType, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:14, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_traitType, filename = "medianStabilityValue_traitType.tif", device = "tiff")

stabilityValue_traitType2 <- ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = traitType2, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:14, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_traitType2, filename = "medianStabilityValue_traitType2.tif", device = "tiff")

stabilityValue_melt_yieldComponents <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType == "Yield-component")
stabilityValue_melt_nonYieldComponents <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType == "Non")

stabilityValue_yieldComponents <- ggplot(stabilityValue_melt_yieldComponents, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_yieldComponents, filename = "medianStabilityValue_yieldComponents.tif", device = "tiff")

stabilityValue_nonYieldComponents <- ggplot(stabilityValue_melt_nonYieldComponents, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_nonYieldComponents, filename = "medianstabilityValue_nonYieldComponents.tif", device = "tiff")

stabilityValue_melt_flowering <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType2 == "Flowering date")
stabilityValue_melt_height <- subset(stabilityValue_allTraits_melt_traitType, stabilityValue_allTraits_melt_traitType$traitType2 == "Height")

stabilityValue_flowering <- ggplot(stabilityValue_melt_flowering, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_flowering, filename = "medianStabilityValue_flowering.tif", device = "tiff")

stabilityValue_height <- ggplot(stabilityValue_melt_height, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_height, filename = "medianStabilityValue_height.tif", device = "tiff")

genos_ordered <- reorder(x = stabilityValue_melt_yieldComponents$Genotype, X = -stabilityValue_melt_yieldComponents$value, FUN = median_)
genos_ordered <- levels(genos_ordered)
genos_ordered <- genos_ordered[1:31]
genos_ordered2 <- c(genos_ordered, genos_ordered)

stabilityValue_flowering <- ggplot(stabilityValue_melt_flowering, mapping = aes(x = fct_inorder(genos_ordered2), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("red", "goldenrod"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  scale_x_discrete(limits = genos_ordered) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "stabilityinabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_flowering, filename = "medianStabilityValue_flowering_orderedYC.tif", device = "tiff")

stabilityValue_height <- ggplot(stabilityValue_melt_height, mapping = aes(x = fct_inorder(genos_ordered2), y = value)) +
  geom_boxplot() +
  geom_jitter(mapping = aes(color = variable, shape = variable)) +
  scale_color_manual(values = c(rep(c("darkgreen", "blue"), 5))) +
  scale_shape_manual(values = c(rep(0:6, 3))) +
  coord_cartesian(xlim = c(0,6)) +
  coord_flip() +
  labs(x = "Environment", y = "stabilityinabilty value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        # legend.position = c(0.15,0.9), 
        legend.title = element_text(" "))
ggsave(plot = stabilityValue_height, filename = "medianStabilityValue_height_orderedYC.tif", device = "tiff")

genos_ordered3 <- c(genos_ordered, genos_ordered, genos_ordered, genos_ordered, 
                    genos_ordered, genos_ordered, genos_ordered, genos_ordered,
                    genos_ordered, genos_ordered, genos_ordered, genos_ordered,
                    genos_ordered, genos_ordered, genos_ordered)

stabilityValue_allTraits_melt_traitType$value2 <- NA
for (i in 1:nrow(stabilityValue_allTraits_melt_traitType)) {
  if (stabilityValue_allTraits_melt_traitType[i,"traitType2"] == "Yield-component") {
    stabilityValue_allTraits_melt_traitType[i,"value2"] <- stabilityValue_allTraits_melt_traitType[i,"value"]
  }
}

stabilityValue_panels <- 
  ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(Genotype, value2, median_), y = value)) +
  geom_boxplot() +
  # geom_jitter(mapping = aes(color = variable, shape = variable)) +
  # scale_color_manual(values = c("red", "goldenrod", rep(c("darkgreen", "blue"), 20))) +
  # scale_shape_manual(values = c(rep(0:6, 20))) +
  coord_flip() +
  facet_grid(~ traitType2) +
  labs(x = "Genotype", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = stabilityValue_panels, filename = "medianStabilityValue_orderedYC_panels.tif", device = "tiff",
       width = 7.5, height = 6.5, units = "in")

stabilityValue_panels_eachTrait <- 
  ggplot(stabilityValue_allTraits_melt_traitType, mapping = aes(x = reorder(Genotype, value2, median_), y = value)) +
  geom_col() +
  coord_flip() +
  facet_grid(~ variable) +
  labs(x = "Genotype", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12, angle = 90), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = stabilityValue_panels_eachTrait, filename = "medianStabilityValue_panelForEachTrait.tif", device = "tiff",
       width = 9, height = 6.5, units = "in")

stabilityValue3 <- ggplot(stabilityValue_allTraits_melt, mapping = aes(x = reorder(x = Genotype, X = -value, FUN = median_), y = value)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = "Environment", y = "Stability value") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.title = element_text(" "))
ggsave(plot = stabilityValue3, filename = "medianStabilityValue3.tif", device = "tiff",
       width = 5, height = 7.75, units = "in")


combinedStabilityBoxplots <- ggarrange(stabilityValue3, stabilityValue_panels,
                                       widths = c(1,1.95),
                                       labels = c("A", "B"),
                                       ncol = 2, nrow = 1)

ggsave(combinedStabilityBoxplots, filename = "combined_stability_boxplots.tif", device = "tiff",
       height = 8.5, width = 11, units = "in")

## stability by release year
lineMetadata <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/gxeInbredLines-Metadata.csv", header = T)
lineMetadata$Year <- as.factor(lineMetadata$Year)
colnames(stabilityRank_allTraits2)[1] <- "Genotype"
stabilityRank_andMetadata <- merge(stabilityRank_allTraits2, lineMetadata, by = "Genotype", )

stabilityRank_andMetadata_melt <- melt(stabilityRank_andMetadata)
stabilityRank_andMetadata_melt$Year <- as.numeric.factor(stabilityRank_andMetadata_melt$Year)

stability2 <- ggplot(stabilityRank_andMetadata_melt, mapping = aes(x = reorder(stabilityRank_andMetadata_melt$Genotype, stabilityRank_andMetadata_melt$Year), y = value)) +
  geom_boxplot() +
  # geom_jitter() +
  labs(x = "Genotype", y = "Stability rank") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 14, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = c(0.15,0.9), legend.title = element_text(" "))
ggsave(plot = stability2, filename = "medianStability2.tif", device = "tiff")

#### Rank correlation of stability of different types of traits/different traits ####
data_stability <- stabilityValue_allTraits[,-1]
cormat <- cor(data_stability, use = "p", method = "pearson")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(data_stability)
head(p.mat[, 1:5])

tiff("stability-correlationPlot.tif")
corrplot(cormat, method = "color", 
         # addCoef.col = "gray40",
         tl.col = "black",  
         type = "upper", diag = F, p.mat = p.mat, 
         mar = c(0,3,0,0), insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")
dev.off()

ggcorrplot_spearman <- ggcorrplot(cor(data_stability, use = "complete.obs", method = "spearman"), p.mat = p.mat, type = "upper",
                                  colors = c("darkred", "white", "steelblue"),
                                  ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_spearman, filename = "stability-ggCorrPlot-spearman.tif", device = "tiff", width = 6, height = 4.5, units = "in")

ggcorrplot_pearson <- ggcorrplot(cor(data_stability, use = "complete.obs", method = "pearson"), p.mat = p.mat, type = "upper",
                                 colors = c("darkred", "white", "steelblue"),
                                 ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_pearson, filename = "stability-ggCorrPlot-pearson.tif", device = "tiff", width = 6, height = 4.5, units = "in")


combinedStabilityBoxplotAndCorrplot <- ggarrange(stabilityValue3, ggcorrplot_spearman,
                                                 widths = c(1,1.8),
                                                 labels = c("A", "B"),
                                                 ncol = 2, nrow = 1)

ggsave(combinedStabilityBoxplotAndCorrplot, filename = "combined_stability_boxplots_and_correlation.tif", device = "tiff",
       height = 8, width = 11, units = "in")

#### Combine BLUPs with line metadata and do boxplots by category ####
lineMetadata <- read.csv("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/gxeInbredLines-Metadata.csv", header = T)
allMetadataAndBlups <- merge(lineMetadata, allBlups, by = "Genotype")
allMetadataAndBlups$X <- NULL; allMetadataAndBlups$X.1 <- NULL; allMetadataAndBlups$X.2 <- NULL; allMetadataAndBlups$X.3 <- NULL; allMetadataAndBlups$X.4 <- NULL; allMetadataAndBlups$X.5 <- NULL
rm(lineMetadata)

boxplot_theme <- theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), 
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())

## Boxplots Silking
boxplot_silking_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = Silking)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(1300, 1600)) +
  boxplot_theme

boxplot_silking_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = Silking)) +
  geom_boxplot() +
  labs(y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(1300, 1600)) +
  theme(axis.text.y = element_blank())

boxplot_silking_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = Silking)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(1300, 1600)) +
  theme(axis.text.y = element_blank())

boxplot_silking_combined <- ggarrange(boxplot_silking_allInbreds, boxplot_silking_group, boxplot_silking_exPVP, 
                                      nrow = 1, ncol = 3)

ggsave(boxplot_silking_combined, filename = "boxplot-silking-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

## Boxplots Plant height
boxplot_plant_height_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = `Plant height`)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(120, 191)) +
  boxplot_theme

boxplot_plant_height_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = `Plant height`)) +
  geom_boxplot() +
  labs(y = "") +
  coord_cartesian(ylim = c(120, 191)) +
  boxplot_theme +
  theme(axis.text.y = element_blank())

boxplot_plant_height_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = `Plant height`)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  coord_cartesian(ylim = c(120, 191)) +
  boxplot_theme +
  theme(axis.text.y = element_blank())


boxplot_plant_height_combined <- ggarrange(boxplot_plant_height_allInbreds, boxplot_plant_height_group, boxplot_plant_height_exPVP, 
                                           nrow = 1, ncol = 3)

ggsave(boxplot_plant_height_combined, filename = "boxplot-plant_height-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

## Boxplots Kernel row number
boxplot_kernel_row_number_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = `Kernel row number`)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(10.5, 17.4)) +
  boxplot_theme

boxplot_kernel_row_number_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = `Kernel row number`)) +
  geom_boxplot() +
  labs(y = "") +
  coord_cartesian(ylim = c(10.5, 17.4)) +
  boxplot_theme +
  theme(axis.text.y = element_blank())

boxplot_kernel_row_number_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = `Kernel row number`)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  coord_cartesian(ylim = c(10.5, 17.4)) +
  boxplot_theme +
  theme(axis.text.y = element_blank()) 


boxplot_kernel_row_number_combined <- ggarrange(boxplot_kernel_row_number_allInbreds, boxplot_kernel_row_number_group, boxplot_kernel_row_number_exPVP, 
                                                nrow = 1, ncol = 3)

ggsave(boxplot_kernel_row_number_combined, filename = "boxplot-kernel_row_number-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

## Boxplots Kernel weight
boxplot_kernel_weight_allInbreds <- ggplot(allMetadataAndBlups, mapping = aes(x = "", y = `Kernel weight`)) +
  geom_boxplot() +
  labs(x = "All inbreds") +
  coord_cartesian(ylim = c(0.185, 0.31)) +
  boxplot_theme

boxplot_kernel_weight_group <- ggplot(allMetadataAndBlups, mapping = aes(x = Group, y = `Kernel weight`)) +
  geom_boxplot() +
  labs(y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(0.185, 0.31)) +
  theme(axis.text.y = element_blank())

boxplot_kernel_weight_exPVP <- ggplot(allMetadataAndBlups, mapping = aes(x = factor(Ex.PVP,levels = c("Non", "Ex-PVP")), y = `Kernel weight`)) +
  geom_boxplot() +
  labs(x = "Ex-PVP", y = "") +
  boxplot_theme +
  coord_cartesian(ylim = c(0.185, 0.31)) +
  theme(axis.text.y = element_blank())

boxplot_kernel_weight_combined <- ggarrange(boxplot_kernel_weight_allInbreds, boxplot_kernel_weight_group, boxplot_kernel_weight_exPVP, 
                                            nrow = 1, ncol = 3)

ggsave(boxplot_kernel_weight_combined, filename = "boxplot-kernel_weight-by_categories.tif", device = "tiff",
       height = 3.5, width = 4.5, units = "in")

boxplot_all_traits_combined <- ggarrange(boxplot_silking_combined, boxplot_plant_height_combined, boxplot_kernel_row_number_combined, boxplot_kernel_weight_combined,
                                         nrow = 4, ncol = 1,
                                         labels = c("A", "B", "C", "D"))
ggsave(boxplot_all_traits_combined, filename = "boxplot-all_traits-combined.tif", device = "tiff",
       height = 12, width = 5, units = "in")


boxplot_all_traits_combined_no_categories <- ggarrange(boxplot_silking_allInbreds, boxplot_plant_height_allInbreds, boxplot_kernel_row_number_allInbreds, boxplot_kernel_weight_allInbreds,
                                                       nrow = 2, ncol = 2,
                                                       labels = c("A", "B", "C", "D"))
ggsave(boxplot_all_traits_combined_no_categories, filename = "boxplot-all_traits-combined_no_categories.tif", device = "tiff",
       height = 5, width = 4, units = "in")


#### Combine slope, MSE, and GGE stability with line metadata ####
lineMetadata <- read.csv("../gxeInbredLines-Metadata.csv", header = T)

allSlopesWithMedian <- read.csv("../all_slopes_with_median.csv", header = T)
colnames(allSlopesWithMedian) <- c("Genotype", "Anthesis", "Silking", 
                                   "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                   "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                   "Kernel length", "Kernel width", "Kernel thickness", "Median")

allMSEWithMedian <- read.csv("../all_MSE_with_median.csv", header = T)
colnames(allMSEWithMedian) <- c("Genotype", "Anthesis", "Silking", 
                                "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                "Kernel length", "Kernel width", "Kernel thickness", "Median")

stabilityValue_allTraits <- read.csv("../stability_value_all_traits.csv", header = T)
colnames(stabilityValue_allTraits) <- c("Genotype", "Anthesis", "Silking", 
                                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                        "Kernel length", "Kernel width", "Kernel thickness")


slopeAndMetadata <- merge(lineMetadata, allSlopesWithMedian, by = "Genotype")

MSEandMetadata <- merge(lineMetadata, allMSEWithMedian, by = "Genotype")

GGEstabilityAndMetadata <- merge(lineMetadata, stabilityValue_allTraits, by = "Genotype")

#### Compare NSS vs SSS ####
# dir.create("./NSSvsSSS")
# dir.create("./NSSvsSSS/boxplotSlope")
# dir.create("./NSSvsSSS/boxplotTypeII")
# dir.create("./NSSvsSSS/boxplotMSE")
# dir.create("./NSSvsSSS/boxplotGGEstability")

#### Slope
sigDiffSlope_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffSlope_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "Slope", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotSlope/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(slopeAndMetadata) - 1)) {
  for (j in 6) {
    # i = 20
    # j = 6
    
    subpopMeans <- aggregate(slopeAndMetadata[,i],by = list(slopeAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(slopeAndMetadata$Genotype,by = list(slopeAndMetadata[,j]),FUN = length)
    
    NSS <- slopeAndMetadata[which(slopeAndMetadata$Group == "NSS"),]
    SSS <- slopeAndMetadata[which(slopeAndMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(slopeAndMetadata[,i]) + (0.00000001 * mean(slopeAndMetadata[,i]))
    
    myBoxplot(slopeAndMetadata, slopeAndMetadata[,j], slopeAndMetadata[,i])
    
    sigDiffSlope_NSSvsSSS <- rbind(sigDiffSlope_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffSlope_NSSvsSSS) <- colnames(slopeAndMetadata[,11:(ncol(slopeAndMetadata) - 1)])


#### Type II Stability (abs(slope-1))
typeIIandMetadata <- slopeAndMetadata %>%
  mutate(Anthesis_type2 = abs(Anthesis - 1)) %>%
  mutate(Silking_type2 = abs(Silking - 1)) %>%
  mutate(PlantHeight_type2 = abs(`Plant height` - 1)) %>%
  mutate(EarHeight_type2 = abs(`Ear height` - 1)) %>%
  mutate(PlotWeight_type2 = abs(`Plot grain weight` - 1)) %>%
  mutate(EarLength_type2 = abs(`Ear length` - 1)) %>%
  mutate(EarWidth_type2 = abs(`Ear width` - 1)) %>%
  mutate(KernelsPerRow_type2 = abs(`Kernels per row` - 1)) %>%
  mutate(KRN_type2 = abs(`Kernel row number` - 1)) %>%
  mutate(KernelWeight_type2 = abs(`Kernel weight` - 1)) %>%
  mutate(KernelArea_type2 = abs(`Kernel area` - 1)) %>%
  mutate(KernelLength_type2 = abs(`Kernel length` - 1)) %>%
  mutate(KernelWidth_type2 = abs(`Kernel width` - 1)) %>%
  mutate(KernelThickness_type2 = abs(`Kernel thickness` - 1)) %>%
  mutate(median_type2 = abs(Median - 1)) %>%
  select(Genotype, Pedigree, Derivation, Ex.PVP, Origin, Group, Company, Year, Release, Accession_Number.Origin,
         Anthesis_type2, Silking_type2, PlantHeight_type2, EarHeight_type2, PlotWeight_type2,
         EarLength_type2, EarWidth_type2, KernelsPerRow_type2, KRN_type2, KernelWeight_type2,
         KernelArea_type2, KernelLength_type2, KernelWidth_type2, KernelThickness_type2, median_type2)

sigDiffTypeII_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffTypeII_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "|slope-1|", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotTypeII/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(typeIIandMetadata) - 1)) {
  for (j in 6) {
    # i = 20
    # j = 6
    
    subpopMeans <- aggregate(typeIIandMetadata[,i],by = list(typeIIandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(typeIIandMetadata$Genotype,by = list(typeIIandMetadata[,j]),FUN = length)
    
    NSS <- typeIIandMetadata[which(typeIIandMetadata$Group == "NSS"),]
    SSS <- typeIIandMetadata[which(typeIIandMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(typeIIandMetadata[,i]) + (0.00000001 * mean(typeIIandMetadata[,i]))
    
    myBoxplot(typeIIandMetadata, typeIIandMetadata[,j], typeIIandMetadata[,i])
    
    sigDiffTypeII_NSSvsSSS <- rbind(sigDiffTypeII_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffTypeII_NSSvsSSS) <- colnames(typeIIandMetadata[,11:(ncol(typeIIandMetadata) - 1)])


#### MSE
sigDiffMSE_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffMSE_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "MSE", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotMSE/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(MSEandMetadata) - 1)) {
  for (j in 6) {
    # i = 12
    # j = 6
    subpopMeans <- aggregate(MSEandMetadata[,i],by = list(MSEandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(MSEandMetadata$Genotype,by = list(MSEandMetadata[,j]),FUN = length)
    
    NSS <- MSEandMetadata[which(MSEandMetadata$Group == "NSS"),]
    SSS <- MSEandMetadata[which(MSEandMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(MSEandMetadata[,i]) + (0.00000001 * mean(MSEandMetadata[,i]))
    
    myBoxplot(MSEandMetadata, MSEandMetadata[,j], MSEandMetadata[,i])
    
    sigDiffMSE_NSSvsSSS <- rbind(sigDiffMSE_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffMSE_NSSvsSSS) <- colnames(MSEandMetadata[,11:(ncol(MSEandMetadata) - 1)])

#### GGE Stability
sigDiffGGEstability_NSSvsSSS <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffGGEstability_NSSvsSSS) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "GGE Stability Value", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./NSSvsSSS/boxplotGGEstability/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:ncol(GGEstabilityAndMetadata)) {
  for (j in 6) {
    # i = 12
    # j = 6
    subpopMeans <- aggregate(GGEstabilityAndMetadata[,i],by = list(GGEstabilityAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(GGEstabilityAndMetadata$Genotype,by = list(GGEstabilityAndMetadata[,j]),FUN = length)
    
    NSS <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Group == "NSS"),]
    SSS <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Group == "SSS"),]
    tTest <- t.test(NSS[,i], SSS[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(GGEstabilityAndMetadata[,i]) + (0.00000001 * mean(GGEstabilityAndMetadata[,i]))
    
    myBoxplot(GGEstabilityAndMetadata, GGEstabilityAndMetadata[,j], GGEstabilityAndMetadata[,i])
    
    sigDiffGGEstability_NSSvsSSS <- rbind(sigDiffGGEstability_NSSvsSSS, myround(tTest$p.value,2))
  }
}

rownames(sigDiffGGEstability_NSSvsSSS) <- colnames(GGEstabilityAndMetadata[,11:ncol(GGEstabilityAndMetadata)])

#### Compare exPVP vs Non ####
# dir.create("./exPVPvsNon")
# dir.create("./exPVPvsNon/boxplotSlope")
# dir.create("./exPVPvsNon/boxplotTypeII")
# dir.create("./exPVPvsNon/boxplotMSE")
# dir.create("./exPVPvsNon/boxplotGGEstability")

#### Slope
sigDiffSlope_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffSlope_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "Slope", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotSlope/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(slopeAndMetadata) - 1)) {
  for (j in 4) {
    # i = 11
    # j = 4
    
    subpopMeans <- aggregate(slopeAndMetadata[,i],by = list(slopeAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(slopeAndMetadata$Genotype,by = list(slopeAndMetadata[,j]),FUN = length)
    
    exPVP <- slopeAndMetadata[which(slopeAndMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- slopeAndMetadata[which(slopeAndMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(slopeAndMetadata[,i]) + (0.00000001 * mean(slopeAndMetadata[,i]))
    
    myBoxplot(slopeAndMetadata, slopeAndMetadata[,j], slopeAndMetadata[,i])
    
    sigDiffSlope_exPVPvsNon <- rbind(sigDiffSlope_exPVPvsNon, myround(tTest$p.value,2))
  }
}

rownames(sigDiffSlope_exPVPvsNon) <- colnames(slopeAndMetadata[,11:(ncol(slopeAndMetadata) - 1)])

#### Type II Stability (abs(slope-1))
sigDiffTypeII_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffTypeII_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "|slope-1|", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotTypeII/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(typeIIandMetadata) - 1)) {
  for (j in 4) {
    # i = 11
    # j = 4
    
    subpopMeans <- aggregate(typeIIandMetadata[,i],by = list(typeIIandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(typeIIandMetadata$Genotype,by = list(typeIIandMetadata[,j]),FUN = length)
    
    exPVP <- typeIIandMetadata[which(typeIIandMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- typeIIandMetadata[which(typeIIandMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value,2))
    yLabPos <- max(typeIIandMetadata[,i]) + (0.00000001 * mean(typeIIandMetadata[,i]))
    
    myBoxplot(typeIIandMetadata, typeIIandMetadata[,j], typeIIandMetadata[,i])
    
    sigDiffTypeII_exPVPvsNon <- rbind(sigDiffTypeII_exPVPvsNon, myround(tTest$p.value,2))
  }
}
rownames(sigDiffTypeII_exPVPvsNon) <- colnames(typeIIandMetadata[,11:(ncol(typeIIandMetadata) - 1)])


#### MSE
sigDiffMSE_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffMSE_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "MSE", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotMSE/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:(ncol(MSEandMetadata) - 1)) {
  for (j in 4) {
    # i = 12
    # j = 4
    subpopMeans <- aggregate(MSEandMetadata[,i],by = list(MSEandMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(MSEandMetadata$Genotype,by = list(MSEandMetadata[,j]),FUN = length)
    
    exPVP <- MSEandMetadata[which(MSEandMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- MSEandMetadata[which(MSEandMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(MSEandMetadata[,i]) + (0.00000001 * mean(MSEandMetadata[,i]))
    
    myBoxplot(MSEandMetadata, MSEandMetadata[,j], MSEandMetadata[,i])
    
    sigDiffMSE_exPVPvsNon <- rbind(sigDiffMSE_exPVPvsNon, myround(tTest$p.value,2))
  }
}

rownames(sigDiffMSE_exPVPvsNon) <- colnames(MSEandMetadata[,11:(ncol(MSEandMetadata) - 1)])

#### GGE Stability
sigDiffGGEstability_exPVPvsNon <- matrix(NA, nrow = 0, ncol = 1)
colnames(sigDiffGGEstability_exPVPvsNon) <- "pvalue"

myBoxplot <- function(data, category, trait) {
  plot <- ggplot(data = data, aes(factor(category), trait)) + 
    geom_boxplot() +
    geom_text(x = 2.2, y = yLabPos, label = pval, size = 8) +
    labs(x = "Subpopulation", y = "GGE Stability Value", title = colnames(data)[i]) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 20), 
          axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 20),
          plot.title = element_text(size = 20),
          panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank())
  ggsave(plot = plot, filename = paste("./exPVPvsNon/boxplotGGEstability/",colnames(data)[i],".tif",sep = ""), device = "tiff")
}

for (i in 11:ncol(GGEstabilityAndMetadata)) {
  for (j in 4) {
    # i = 12
    # j = 4
    subpopMeans <- aggregate(GGEstabilityAndMetadata[,i],by = list(GGEstabilityAndMetadata[,j]),mean,na.rm = T)
    subpopCount <- aggregate(GGEstabilityAndMetadata$Genotype,by = list(GGEstabilityAndMetadata[,j]),FUN = length)
    
    exPVP <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Ex.PVP == "Ex-PVP"),]
    Non <- GGEstabilityAndMetadata[which(GGEstabilityAndMetadata$Ex.PVP == "Non"),]
    tTest <- t.test(exPVP[,i], Non[,i])
    pval <- paste0("p = ", myround(tTest$p.value, 2))
    yLabPos <- max(GGEstabilityAndMetadata[,i]) + (0.00000001 * mean(GGEstabilityAndMetadata[,i]))
    
    myBoxplot(GGEstabilityAndMetadata, GGEstabilityAndMetadata[,j], GGEstabilityAndMetadata[,i])
    
    sigDiffGGEstability_exPVPvsNon <- rbind(sigDiffGGEstability_exPVPvsNon, myround(tTest$p.value,2))
  }
}

rownames(sigDiffGGEstability_exPVPvsNon) <- colnames(GGEstabilityAndMetadata[,11:ncol(GGEstabilityAndMetadata)])


sigDiffSlope_NSSvsSSS <- as.data.frame(sigDiffSlope_NSSvsSSS)
sigDiffTypeII_NSSvsSSS <- as.data.frame(sigDiffTypeII_NSSvsSSS)
sigDiffMSE_NSSvsSSS <- as.data.frame(sigDiffMSE_NSSvsSSS)
sigDiffGGEstability_NSSvsSSS <- as.data.frame(sigDiffGGEstability_NSSvsSSS)

sigDiffSlope_exPVPvsNon <- as.data.frame(sigDiffSlope_exPVPvsNon)
sigDiffTypeII_exPVPvsNon <- as.data.frame(sigDiffTypeII_exPVPvsNon)
sigDiffMSE_exPVPvsNon <- as.data.frame(sigDiffMSE_exPVPvsNon)
sigDiffGGEstability_exPVPvsNon <- as.data.frame(sigDiffGGEstability_exPVPvsNon)

sigDiff <- cbind(sigDiffSlope_NSSvsSSS, sigDiffTypeII_NSSvsSSS, sigDiffMSE_NSSvsSSS, sigDiffGGEstability_NSSvsSSS,
                 sigDiffSlope_exPVPvsNon, sigDiffTypeII_exPVPvsNon, sigDiffMSE_exPVPvsNon, sigDiffGGEstability_exPVPvsNon)

colnames(sigDiff) <- c("slope_NSSvsSSS", "typeII_NSSvsSSS", "MSE_NSSvsSSS", "GGE_NSSvsSSS",
                       "slope_exPVPvsNon", "typeII_exPVPvsNon", "MSE_exPVPvsNon", "GGE_exPVPvsNon")

sigDiff[,1] <- as.numeric.factor(sigDiff[,1])
sigDiff[,2] <- as.numeric.factor(sigDiff[,2])
sigDiff[,3] <- as.numeric.factor(sigDiff[,3])
sigDiff[,4] <- as.numeric.factor(sigDiff[,4])
sigDiff[,5] <- as.numeric.factor(sigDiff[,5])
sigDiff[,6] <- as.numeric.factor(sigDiff[,6])
sigDiff[,7] <- as.numeric.factor(sigDiff[,7])
sigDiff[,8] <- as.numeric.factor(sigDiff[,8])

write.csv(sigDiff, "sigDiffs_stability_genotypeCategories.csv")

#### Rank correlation of measures of stability ####
colnames(allSlopesWithMedian)[2:ncol(allSlopesWithMedian)] <- paste0("slope_", colnames(allSlopesWithMedian)[2:ncol(allSlopesWithMedian)])
colnames(allMSEWithMedian)[2:ncol(allMSEWithMedian)] <- paste0("MSE_", colnames(allMSEWithMedian)[2:ncol(allMSEWithMedian)])
colnames(stabilityValue_allTraits_withMed)[2:ncol(stabilityValue_allTraits_withMed)] <- paste0("GGE_", colnames(stabilityValue_allTraits_withMed)[2:ncol(stabilityValue_allTraits_withMed)])

typeIIStability <- typeIIandMetadata[,c(1,11:25)]

stabilityMeasures <- merge(allSlopesWithMedian, allMSEWithMedian, by = "Genotype")
stabilityMeasures <- merge(stabilityMeasures, stabilityValue_allTraits_withMed, by = "Genotype")
stabilityMeasures <- merge(stabilityMeasures, typeIIStability, by = "Genotype")

corr_slope_MSE <- cor(stabilityMeasures$slope_Median, stabilityMeasures$MSE_Median, method = "spearman")
corr_slope_GGE <- cor(stabilityMeasures$slope_Median, stabilityMeasures$GGE_median, method = "spearman")

corr_typeII_MSE <- cor(stabilityMeasures$median_type2, stabilityMeasures$MSE_Median, method = "spearman")
corr_typeII_GGE <- cor(stabilityMeasures$median_type2, stabilityMeasures$GGE_median, method = "spearman")

corr_MSE_GGE <- cor(stabilityMeasures$MSE_Median, stabilityMeasures$GGE_median, method = "spearman")


cormat <- cor(stabilityMeasures[,-1], use = "p", method = "spearman")
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(stabilityMeasures[,-1])

ggcorrplot <- 
  ggcorrplot(cor(stabilityMeasures[,-1]), p.mat = cor.mtest(stabilityMeasures[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot, filename = "ggCorrPlot-stabilityMeasures.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

allSlopes <- allSlopesWithMedian[,-16]
colnames(allSlopes) <- gsub(pattern = "slope_", replacement = "", x = colnames(allSlopes))
colnames(allSlopes) <- gsub(pattern = "_st", replacement = "", x = colnames(allSlopes))

ggcorrplot_slopeAmongTraits <- 
  ggcorrplot(cor(allSlopes[,-1]), p.mat = cor.mtest(allSlopes[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_slopeAmongTraits, filename = "ggCorrPlot-slopeAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

typeIIStability <- typeIIStability[,-ncol(typeIIStability)]
colnames(typeIIStability) <- gsub(pattern = "_type2", replacement = "", x = colnames(typeIIStability))

ggcorrplot_typeIIAmongTraits <- 
  ggcorrplot(cor(typeIIStability[,-1]), p.mat = cor.mtest(typeIIStability[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_typeIIAmongTraits, filename = "ggCorrPlot-typeIIAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

MSEstability <- allMSEWithMedian[,-16]
colnames(MSEstability) <- gsub(pattern = "MSE_", replacement = "", x = colnames(MSEstability))
colnames(MSEstability) <- gsub(pattern = "_st", replacement = "", x = colnames(MSEstability))

ggcorrplot_MSEAmongTraits <- 
  ggcorrplot(cor(MSEstability[,-1]), p.mat = cor.mtest(MSEstability[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_MSEAmongTraits, filename = "ggCorrPlot-MSEAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")

GGEstability <- stabilityValue_allTraits_withMed[-16]
colnames(GGEstability) <- gsub(pattern = "GGE_", replacement = "", x = colnames(GGEstability))
colnames(GGEstability) <- gsub(pattern = "\\(GDU)", replacement = "", x = colnames(GGEstability))

ggcorrplot_GGEAmongTraits <- 
  ggcorrplot(cor(GGEstability[,-1]), p.mat = cor.mtest(GGEstability[,-1]), type = "upper",
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 8),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_GGEAmongTraits, filename = "ggCorrPlot-GGEAmongTraits.tif", device = "tiff",
       width = 6, height = 4.5, units = "in")


combinedStabilityCorrAmongTraits <- ggarrange(ggcorrplot_typeIIAmongTraits, ggcorrplot_MSEAmongTraits,
                                              nrow = 1, ncol = 2,
                                              labels = c("A", "B"))
ggsave(combinedStabilityCorrAmongTraits, filename = "combinedStabilityCorrAmongTraits.tif", device = "tiff",
       width = 9, height = 3.5, units = "in")


#### Plot of genotypes that are in 5 most or 5 least stable for slope, MSE, GGE ####
traits <- c("Anthesis", "Silking", 
            "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
            "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
            "Kernel length", "Kernel width", "Kernel thickness")

colnames(stabilityMeasures) <- c("Genotype", "slope_Anthesis", "slope_Silking", "slope_Plant height", "slope_Ear height", 
                                 "slope_Plot grain weight", "slope_Ear length", "slope_Ear width", "slope_Kernels per row", "slope_Kernel row number", 
                                 "slope_Kernel weight", "slope_Kernel area", "slope_Kernel length", "slope_Kernel width", "slope_Kernel thickness",
                                 "slope_medianSlope", "MSE_Anthesis", "MSE_Silking", "MSE_Plant height", "MSE_Ear height",     
                                 "MSE_Plot grain weight", "MSE_Ear length", "MSE_Ear width", "MSE_Kernels per row", "MSE_Kernel row number",           
                                 "MSE_Kernel weight", "MSE_Kernel area", "MSE_Kernel length", "MSE_Kernel width", "MSE_Kernel thickness",  
                                 "MSE_medianMSE", "GGE_Anthesis", "GGE_Silking", "GGE_Plant height", "GGE_Ear height",          
                                 "GGE_Plot grain weight", "GGE_Ear length", "GGE_Ear width", "GGE_Kernels per row", "GGE_Kernel row number",   
                                 "GGE_Kernel weight", "GGE_Kernel area", "GGE_Kernel length", "GGE_Kernel width", "GGE_Kernel thickness",    
                                 "GGE_median", "type2_Anthesis", "type2_Silking", "type2_Plant height", "type2_Ear height", 
                                 "type2_Plot grain weight", "type2_Ear length", "type2_Ear width", "type2_Kernels per row", "type2_Kernel row number", 
                                 "type2_Kernel weight", "type2_Kernel area", "type2_Kernel length", "type2_Kernel width", "type2_Kernel thickness",
                                 "type2_median")
stabilityMeasures$MSE_medianMSE <- NULL
stabilityMeasures$GGE_median <- NULL
stabilityMeasures$slope_medianSlope <- NULL
stabilityMeasures$type2_median <- NULL

stabilityMeasures_melt <- melt(stabilityMeasures, id.vars = "Genotype")
stabilityMeasures_melt$variable1 <- gsub(pattern = "_.*", replacement = "", stabilityMeasures_melt$variable)
stabilityMeasures_melt$variable2 <- gsub(pattern = "slope_", replacement = "", stabilityMeasures_melt$variable)
stabilityMeasures_melt$variable2 <- gsub(pattern = "MSE_", replacement = "", stabilityMeasures_melt$variable2)
stabilityMeasures_melt$variable2 <- gsub(pattern = "GGE_", replacement = "", stabilityMeasures_melt$variable2)
stabilityMeasures_melt$variable <- NULL
stabilityMeasures_melt <- stabilityMeasures_melt[,c(1,3,4,2)]
colnames(stabilityMeasures_melt) <- c("Genotype", "stability measure", "trait", "value")

stabilityMeasures_melt$`stability measure` <- as.factor(stabilityMeasures_melt$`stability measure`)
stabilityMeasures_melt$trait <- as.factor(stabilityMeasures_melt$trait)

####

least_and_most_stable_genos <- matrix(NA, nrow = 0, ncol = 4)
colnames(least_and_most_stable_genos) <- c("measure", "trait", "least", "most")

m = "GGE"
measureStability <- stabilityMeasures_melt %>%
  filter(`stability measure` == m)

for (t in 1:length(traits)) {
  # t = 1
  
  traitStability <- measureStability %>%
    filter(measureStability$trait == traits[t])
  
  # most stable = smallest distance from AEC
  mostStable <- traitStability %>%
    top_n(-5) %>%
    .$Genotype
  
  # least stable = greatest distance from AEC
  leastStable <- traitStability %>%
    top_n(5) %>%
    .$Genotype
  
  genos_rows <- c(rep(m, 5), rep(traits[t], 5), as.character(leastStable), as.character(mostStable))
  genos_matrix <- matrix(genos_rows, nrow = 5, ncol = 4)
  least_and_most_stable_genos <- rbind(least_and_most_stable_genos, genos_matrix)
}

m = "MSE"
measureStability <- stabilityMeasures_melt %>%
  filter(`stability measure` == m)

for (t in 1:length(traits)) {
  # t = 1
  
  traitStability <- measureStability %>%
    filter(measureStability$trait == traits[t])
  
  # most stable = smallest distance from AEC
  mostStable <- traitStability %>%
    top_n(-5) %>%
    .$Genotype
  
  # least stable = greatest distance from AEC
  leastStable <- traitStability %>%
    top_n(5) %>%
    .$Genotype
  
  genos_rows <- c(rep(m, 5), rep(traits[t], 5), as.character(leastStable), as.character(mostStable))
  genos_matrix <- matrix(genos_rows, nrow = 5, ncol = 4)
  least_and_most_stable_genos <- rbind(least_and_most_stable_genos, genos_matrix)
}

m = "slope"
measureStability <- stabilityMeasures_melt %>%
  filter(`stability measure` == m)
measureStability$value2 <- abs(measureStability$value - 1)

for (t in 1:length(traits)) {
  # t = 1
  
  traitStability <- measureStability %>%
    filter(measureStability$trait == traits[t])
  
  # most stable = smallest distance from AEC
  mostStable <- traitStability %>%
    top_n(-5) %>%
    .$Genotype
  
  # least stable = greatest distance from AEC
  leastStable <- traitStability %>%
    top_n(5) %>%
    .$Genotype
  
  genos_rows <- c(rep(m, 5), rep(traits[t], 5), as.character(leastStable), as.character(mostStable))
  genos_matrix <- matrix(genos_rows, nrow = 5, ncol = 4)
  least_and_most_stable_genos <- rbind(least_and_most_stable_genos, genos_matrix)
}

least_and_most_stable_genos <- as.data.frame(least_and_most_stable_genos)
least_and_most_stable_genos_melt <- melt(least_and_most_stable_genos, id.var = c("measure", "trait"))

write.csv(least_and_most_stable_genos, "least_and_most_stable_genos.csv", row.names = F)

myTable <- table(least_and_most_stable_genos_melt$measure, least_and_most_stable_genos_melt$variable, least_and_most_stable_genos_melt$value)
myTable_melt <- melt(myTable)

for (i in 1:nrow(myTable_melt)) {
  if (myTable_melt[i,"Var2"] == "least") {
    myTable_melt[i,"value"] <- -myTable_melt[i,"value"]
  }
}

colnames(myTable_melt) <- c("measure", "stability", "Genotype", "value")

stable_genos <-
  ggplot(myTable_melt, aes(x = Genotype, y = value)) +
  geom_bar(stat = "identity", fill = "gray50") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos, filename = paste("./most_least_stable.tif",sep = ""),
       device = "tiff", width = 5, height = 8, units = "in")


myTable_melt_with_metadata <- merge(lineMetadata, myTable_melt)

stable_genos_byYear <-
  ggplot(myTable_melt_with_metadata, aes(x = reorder(Genotype, Year), y = value)) +
  geom_bar(stat = "identity", fill = "gray50") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = c(0.225,0.9), legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos_byYear, filename = paste("./most_least_stable_byYear.tif",sep = ""),
       device = "tiff", width = 5, height = 8, units = "in")

stable_genos_byExPVP <-
  ggplot(myTable_melt_with_metadata, aes(x = reorder(Genotype, Year), y = value, fill = Ex.PVP)) +
  geom_bar(stat = "identity") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos_byExPVP, filename = paste("./most_least_stable_byExPVP.tif",sep = ""),
       device = "tiff", width = 6, height = 8, units = "in")


myTable_melt_with_metadata$Genotype <- factor(myTable_melt_with_metadata$Genotype,
                                              levels = unique(myTable_melt_with_metadata[order(myTable_melt_with_metadata$Group),"Genotype"]))


stable_genos_byGroup <-
  ggplot(myTable_melt_with_metadata, aes(x = reorder(Genotype, Group), y = value, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(x = "Genotype", y = "") +
  scale_y_continuous(breaks = c(-8, -4, 0, 4),
                     labels = c("8", "4", "0", "4")) +
  geom_abline(slope = 0, intercept = 0, size = 1.2) +
  annotate("text", label = "Most stable", x = 16, y = 8, size = 4, family = "Courier") +
  annotate("text", label = "Least stable", x = 16, y = -9, size = 4, family = "Courier") +
  facet_grid(measure~.) +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 0), axis.text.x = element_text(size = 10, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 10),
        panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave(plot = stable_genos_byGroup, filename = paste("./most_least_stable_byGroup.tif",sep = ""),
       device = "tiff", width = 6, height = 8, units = "in")


#### Correlation of performance (BLUPs) and stability ####
BLUPs_stabilityMeasures <- merge(allBlups, stabilityMeasures, by = "Genotype")

i = 2
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_anthesis <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_anthesis$typeII <- abs(BLUPs_stabilityMeasures_anthesis$slope_Anthesis - 1)

i = 3
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_silking <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_silking$typeII <- abs(BLUPs_stabilityMeasures_silking$slope_Silking - 1)

i = 4
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_plantHeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_plantHeight$typeII <- abs(BLUPs_stabilityMeasures_plantHeight$`slope_Plant height` - 1)

i = 5
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_earHeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_earHeight$typeII <- abs(BLUPs_stabilityMeasures_earHeight$`slope_Ear height` - 1)

i = 6
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_plotWeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_plotWeight$typeII <- abs(BLUPs_stabilityMeasures_plotWeight$`slope_Plot grain weight` - 1)

i = 7
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_earLength <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_earLength$typeII <- abs(BLUPs_stabilityMeasures_earLength$`slope_Ear length` - 1)

i = 8
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_earWidth <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_earWidth$typeII <- abs(BLUPs_stabilityMeasures_earWidth$`slope_Ear width` - 1)

i = 9
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelsPerRow <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelsPerRow$typeII <- abs(BLUPs_stabilityMeasures_kernelsPerRow$`slope_Kernels per row` - 1)

i = 10
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelRowNumber <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelRowNumber$typeII <- abs(BLUPs_stabilityMeasures_kernelRowNumber$`slope_Kernel row number` - 1)

i = 11
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelWeight <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelWeight$typeII <- abs(BLUPs_stabilityMeasures_kernelWeight$`slope_Kernel weight` - 1)

i = 12
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelArea <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelArea$typeII <- abs(BLUPs_stabilityMeasures_kernelArea$`slope_Kernel area` - 1)

i = 13
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelLength <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelLength$typeII <- abs(BLUPs_stabilityMeasures_kernelLength$`slope_Kernel length` - 1)

i = 14
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelWidth <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelWidth$typeII <- abs(BLUPs_stabilityMeasures_kernelWidth$`slope_Kernel width` - 1)

i = 15
columns <- c(1, i, (i + 14), (i + 28), (i + 42))
BLUPs_stabilityMeasures_kernelThickness <- BLUPs_stabilityMeasures %>%
  select(columns)
BLUPs_stabilityMeasures_kernelThickness$typeII <- abs(BLUPs_stabilityMeasures_kernelThickness$`slope_Kernel thickness` - 1)

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

ggcorrplot_anthesis_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_anthesis[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_anthesis[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_anthesis_BLUPs_stability, filename = "ggcorrplot_anthesis_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_silking_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_silking[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_silking[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_silking_BLUPs_stability, filename = "ggcorrplot_silking_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_plantHeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_plantHeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_plantHeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_plantHeight_BLUPs_stability, filename = "ggcorrplot_plantHeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_earHeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_earHeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_earHeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_earHeight_BLUPs_stability, filename = "ggcorrplot_earHeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_plotWeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_plotWeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_plotWeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_plotWeight_BLUPs_stability, filename = "ggcorrplot_plotWeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_earLength_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_earLength[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_earLength[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_earLength_BLUPs_stability, filename = "ggcorrplot_earLength_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_earWidth_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_earWidth[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_earWidth[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_earWidth_BLUPs_stability, filename = "ggcorrplot_earWidth_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelsPerRow_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelsPerRow[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelsPerRow[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelsPerRow_BLUPs_stability, filename = "ggcorrplot_kernelsPerRow_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelRowNumber_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelRowNumber[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelRowNumber[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelRowNumber_BLUPs_stability, filename = "ggcorrplot_kernelRowNumber_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelWeight_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelWeight[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelWeight[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelWeight_BLUPs_stability, filename = "ggcorrplot_kernelWeight_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelArea_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelArea[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelArea[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelArea_BLUPs_stability, filename = "ggcorrplot_kernelArea_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")

ggcorrplot_kernelLength_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelLength[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelLength[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelLength_BLUPs_stability, filename = "ggcorrplot_kernelLength_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")


ggcorrplot_kernelWidth_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelWidth[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelWidth[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelWidth_BLUPs_stability, filename = "ggcorrplot_kernelWidth_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")


ggcorrplot_kernelThickness_BLUPs_stability <-
  ggcorrplot(cor(BLUPs_stabilityMeasures_kernelThickness[,-1]), p.mat = cor.mtest(BLUPs_stabilityMeasures_kernelThickness[,-1]), 
             type = "upper", lab = T,
             colors = c("darkred", "white", "steelblue"),
             ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot_kernelThickness_BLUPs_stability, filename = "ggcorrplot_kernelThickness_BLUPs_stability.tif", 
       device = "tiff", width = 5, height = 5, units = "in")




