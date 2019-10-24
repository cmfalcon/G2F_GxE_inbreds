##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-Stability analysis in 
##                 non-Midwest locations
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


#### Stability analysis in non-Midwest locations ####
# dir.create("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_non_midwest")
setwd("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability_in_non_midwest")
non_midwest_pheno_data <- selectedPhenoData %>%
  filter(Enviro == "DE1_14" |
           Enviro == "DE1_15" |
           Enviro == "GA1_15" |
           Enviro == "GA2_14" |
           Enviro == "KS1_15" |
           Enviro == "NC1_14" |
           Enviro == "NC1_15" |
           Enviro == "NE1_14" |
           Enviro == "NY1_14" |
           Enviro == "NY1_15" |
           Enviro == "NY2_15" |
           Enviro == "PA1_14" |
           Enviro == "SD1_15" |
           Enviro == "TX1_14" |
           Enviro == "TX1_15" |
           Enviro == "TX2_14" |
           Enviro == "TX2_15" |
           Enviro == "TX3_15")

#### Linear regression models and plots (slope and MSE) --- TRAITS STANDARDIZED ####
non_midwest_pheno_data_standardized <- non_midwest_pheno_data
non_midwest_pheno_data_standardized <- non_midwest_pheno_data_standardized %>% transmute(
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

non_midwest_pheno_data_standardized <- cbind(non_midwest_pheno_data[,1:4], non_midwest_pheno_data_standardized)


#### Calculate environment means
for (i in 5:ncol(non_midwest_pheno_data_standardized)) {
  # i = 5
  means_by_exp <- aggregate(non_midwest_pheno_data_standardized,
                            by = list(non_midwest_pheno_data_standardized$Enviro), FUN = mean, na.rm = T)
  means_by_exp$Enviro <- means_by_exp$Group.1
  means_by_exp <- means_by_exp[,-1]
}

####

names <- colnames(means_by_exp[ ,5:ncol(means_by_exp)])

index <- 1

inbreds <- unique(non_midwest_pheno_data_standardized$Genotype)
traits <- colnames(non_midwest_pheno_data_standardized[ ,5:ncol(non_midwest_pheno_data_standardized)]) 

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
    inbredDat <- non_midwest_pheno_data_standardized[non_midwest_pheno_data_standardized$Genotype == inbred, c("Enviro", trait)]
    
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
  spread <- range(non_midwest_pheno_data_standardized[,trait], na.rm = T)
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

#### GGE biplot analysis --- TRAITS STANDARDIZED ####
dir.create("./GGEbiplots_standardized")
dir.create("./boxplotsByEnviro_standardized_")

non_midwest_pheno_data_standardized$Plot = as.factor(non_midwest_pheno_data_standardized$Plot)
non_midwest_pheno_data_standardized$Enviro = as.factor(non_midwest_pheno_data_standardized$Enviro)
non_midwest_pheno_data_standardized$Genotype = as.factor(non_midwest_pheno_data_standardized$Genotype)
non_midwest_pheno_data_standardized$Rep = as.factor(non_midwest_pheno_data_standardized$Rep)

non_midwest_pheno_data_standardized$Enviro <- gsub("_20", "_", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("AZI", "AZ", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("DEI", "DE", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("GAI", "GA", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("IAI", "IA", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("ILI", "IL", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("INI", "IN", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("KSI", "KS", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("MNI", "MN", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("MOI", "MO", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("NCI", "NC", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("NYI", "NY", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("PAI", "PA", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("SDI", "SD", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("TXI", "TX", non_midwest_pheno_data_standardized$Enviro, fixed = T)
non_midwest_pheno_data_standardized$Enviro <- gsub("WII", "WI", non_midwest_pheno_data_standardized$Enviro, fixed = T)

#### GGE model analysis ####
enviros <- sort(unique(non_midwest_pheno_data_standardized$Enviro), decreasing = F)
genos <- sort(unique(non_midwest_pheno_data_standardized$Genotype), decreasing = F)

discrimRank_allTraits <- matrix(enviros, nrow = length(unique(non_midwest_pheno_data_standardized$Enviro)), ncol = 1)
stabilityRank_allTraits <- matrix(genos, nrow = length(unique(non_midwest_pheno_data_standardized$Genotype)), ncol = 1)

colnames(discrimRank_allTraits) <- "enviros"
colnames(stabilityRank_allTraits) <- "genos"

GGE_PC_varExpl <- matrix(NA, nrow = 0, ncol = 4)
colnames(GGE_PC_varExpl) <- c("Trait", "PC1", "PC2", "PC1+PC2")

for (t in c(5:ncol(non_midwest_pheno_data_standardized))) {
  # t = 18
  
  #### calculate means across reps and create a matrix of Gs by Es for each trait ####
  trait = non_midwest_pheno_data_standardized[,t]
  df <- non_midwest_pheno_data_standardized[!is.na(trait),]
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
  
  # # write.csv(data, file = paste0("gxe1415-", colnames(non_midwest_pheno_data_standardized[t]), ".csv"))
  #
  # 
  # #### estimate missing G, E data ####
  # a <- sum(is.na(dataWithNAs))
  # c <- (nrow(dataWithNAs) * ncol(dataWithNAs)) - a
  
  # data[is.na(data)] <- 0
  # data <- as.matrix(data)
  # 
  # test <- 1
  # while (test > 0.01) {
  #   DDt <- data %*% t(data)
  #   DtD <- t(data) %*% data
  #   
  #   eig.result <- eigen(DDt)
  #   str(eig.result)
  #   eig.vec.geno <- eig.result$vectors
  #   lambda <- eig.result$values
  #   
  #   eig.result2 <- eigen(DtD)
  #   str(eig.result2)
  #   eig.vec.enviro <- eig.result2$vectors
  #   lambda2 <- eig.result2$values
  #   
  #   S <- sqrt(lambda2)
  #   
  #   rankApprox <- matrix(0, nrow = nrow(DDt), ncol = nrow(DtD))
  #   for (i in 1:ncol(DDt)) {
  #     for (j in 1:ncol(DtD)) {
  #       rankApprox[i,j] <- (S[1] * eig.vec.geno[i,1] * eig.vec.enviro[j,1]) + (S[2] * eig.vec.geno[i,2] * eig.vec.enviro[j,2])
  #     }
  #   }
  #   
  #   for (m in 1:nrow(dataWithNAs)) {
  #     for (n in 1:ncol(dataWithNAs)) {
  #       if (is.na(dataWithNAs[m,n])) {
  #         newData[m,n] <- rankApprox[m,n]
  #       }
  #     }
  #   }
  #   
  #   b = 0
  #   for (m in 1:nrow(dataWithNAs)) {
  #     for (n in 1:ncol(dataWithNAs)) {
  #       if (is.na(dataWithNAs[m,n])) {
  #         b <- b + (newData[m,n] - data[m,n])^2
  #       }
  #     }
  #   }
  #   
  #   e = 0
  #   for (m in 1:nrow(dataWithNAs)) {
  #     for (n in 1:ncol(dataWithNAs)) {
  #       if (!is.na(dataWithNAs[m,n])) {
  #         e <- e + (dataWithNAs[m,n])^2
  #       }
  #     }
  #   }
  #   
  #   d <- sqrt((1/a) * b)
  #   
  #   yBar <- sqrt((1/c) * e)
  #   
  #   test <- d/yBar
  #   
  #   data <- as.matrix(newData)
  # }
  
  #### create plots ####
  GGE1 <- GGEModel(data) #Produces genotype plus genotype-by-environment model from a 2-way table of means
  # basic <- GGEPlot(GGE1, type = 1) # basic biplot
  # ggsave(plot = basic, filename = paste0("ggeBiplot-basic-", colnames(non_midwest_pheno_data_standardized[t]), ".tif"), device = "tiff")
  
  # GGEPlot(GGE1, type = 2, selectedE = "E1") # examine a particular environment
  # GGEPlot(GGE1, type = 3, selectedG = "G1") # examine a particular genotype
  # GGEPlot(GGE1, type = 4) # relationship among environments
  # GGEPlot(GGE1, type = 5, selectedG1 = "G1", selectedG2 = "G2") # compare 2 genotypes
  whichWonWhere <- GGEPlot(GGE1, type = 6) # which won where/what
  ggsave(plot = whichWonWhere, filename = paste0("./GGEbiplots_standardized/WWW-", colnames(non_midwest_pheno_data_standardized[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
  discrimVsRepresent <- GGEPlot(GGE1, type = 7, sizeGen = 2, sizeEnv = 3) # discrimination vs. representativeness
  ggsave(plot = discrimVsRepresent, filename = paste0("./GGEbiplots_standardized/DvR-", colnames(non_midwest_pheno_data_standardized[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
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
  colnames(rankVectorLengths)[2] <- paste0("vectorLength_", colnames(non_midwest_pheno_data_standardized[t]))
  colnames(rankVectorLengths)[3] <- paste0("rank_", colnames(non_midwest_pheno_data_standardized[t]))
  
  discrimRank_allTraits <- merge(discrimRank_allTraits, rankVectorLengths, by = "enviros", all = T)
  
  # GGEPlot(GGE1, type = 8) # ranking environments
  meanVsStability <- GGEPlot(GGE1, type = 9, sizeEnv = 2, sizeGen = 3) # mean vs. stability
  ggsave(plot = meanVsStability, filename = paste0("./GGEbiplots_standardized/MvS-", colnames(non_midwest_pheno_data_standardized[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
  
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
  colnames(rankGenoLengths)[2] <- paste0("genoLength_", colnames(non_midwest_pheno_data_standardized[t]))
  colnames(rankGenoLengths)[3] <- paste0("rank_", colnames(non_midwest_pheno_data_standardized[t]))
  
  stabilityRank_allTraits <- merge(stabilityRank_allTraits, rankGenoLengths, by = "genos", all = T)
  
  # GGEPlot(GGE1, type = 10) # ranking genotypes
  
  #### Pull out var explained by PC1, PC2, and PC1+PC2
  PC1_varExpl <- GGE1$varexpl[1]
  PC2_varExpl <- GGE1$varexpl[2]
  sumPC1PC2_varExpl <- sum(PC1_varExpl, PC2_varExpl)
  
  row <- c(colnames(non_midwest_pheno_data_standardized)[t], PC1_varExpl, PC2_varExpl, sumPC1PC2_varExpl)
  
  GGE_PC_varExpl <- rbind(GGE_PC_varExpl, row)
}

#### Stability  ####
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


