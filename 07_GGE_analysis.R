##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-GGE analysis
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(dplyr)
library(plyr)
library(GGEBiplots)
library(ggplot2)
library(data.table)
library(broman)

mean_   <- function(...) mean(..., na.rm = T)
rowMedian <- function(x, na.rm = FALSE) apply(x, 1, median, na.rm = na.rm)

#### Read in, format, and organize data ####
selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")

#### standardize data ####
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
colnames(selectedPhenoData_standardized) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")


#### GGE model analysis ####
enviros <- sort(unique(selectedPhenoData_standardized$Enviro), decreasing = F)
genos <- sort(unique(selectedPhenoData_standardized$Genotype), decreasing = F)

discrimRank_allTraits <- matrix(enviros, nrow = length(unique(selectedPhenoData_standardized$Enviro)), ncol = 1)
stabilityRank_allTraits <- matrix(genos, nrow = length(unique(selectedPhenoData_standardized$Genotype)), ncol = 1)

colnames(discrimRank_allTraits) <- "enviros"
colnames(stabilityRank_allTraits) <- "genos"

GGE_PC_varExpl <- matrix(NA, nrow = 0, ncol = 4)
colnames(GGE_PC_varExpl) <- c("Trait", "PC1", "PC2", "PC1+PC2")

for (t in c(5:ncol(selectedPhenoData_standardized))) {
  # t = 6
  
  #### calculate means across reps and create a matrix of Gs by Es for each trait ####
  trait = selectedPhenoData_standardized[,t]
  df <- selectedPhenoData_standardized[!is.na(trait),]
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
  ggsave(plot = whichWonWhere, filename = paste0("./GGEbiplots_standardized/WWW-", colnames(selectedPhenoData_standardized[t]), ".tif"),
         device = "tiff", width = 4, height = 4, units = "in")
  
  discrimVsRepresent <- GGEPlot(GGE1, type = 7, sizeGen = 2, sizeEnv = 3) # discrimination vs. representativeness
  ggsave(plot = discrimVsRepresent, filename = paste0("./GGEbiplots_standardized/DvR-", colnames(selectedPhenoData_standardized[t]), ".tif"),
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
  colnames(rankVectorLengths)[2] <- paste0("vectorLength_", colnames(selectedPhenoData_standardized[t]))
  colnames(rankVectorLengths)[3] <- paste0("rank_", colnames(selectedPhenoData_standardized[t]))
  
  discrimRank_allTraits <- merge(discrimRank_allTraits, rankVectorLengths, by = "enviros", all = T)
  
  #### create mean vs stability plots ####
  meanVsStability <- GGEPlot(GGE2, type = 9, sizeEnv = 2, sizeGen = 3) # mean vs. stability
  ggsave(plot = meanVsStability, filename = paste0("./GGEbiplots_standardized/MvS-", colnames(selectedPhenoData_standardized[t]), ".tif"),
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
  colnames(rankGenoLengths)[2] <- paste0("genoLength_", colnames(selectedPhenoData_standardized[t]))
  colnames(rankGenoLengths)[3] <- paste0("rank_", colnames(selectedPhenoData_standardized[t]))
  
  stabilityRank_allTraits <- merge(stabilityRank_allTraits, rankGenoLengths, by = "genos", all = T)
  
  #### Output var explained by PC1, PC2, and PC1+PC2
  PC1_varExpl <- GGE1$varexpl[1]
  PC2_varExpl <- GGE1$varexpl[2]
  sumPC1PC2_varExpl <- sum(PC1_varExpl, PC2_varExpl)
  
  row <- c(colnames(selectedPhenoData_standardized)[t], PC1_varExpl, PC2_varExpl, sumPC1PC2_varExpl)
  
  GGE_PC_varExpl <- rbind(GGE_PC_varExpl, row)
}

write.csv(discrimRank_allTraits, "discrim_rank_all_traits.csv", row.names = F)
write.csv(stabilityRank_allTraits, "stability_rank_all_traits.csv", row.names = F)

### Create GGE stability summary

allStability <- stabilityRank_allTraits %>%
  select(genos, `genoLength_Anthesis (GDU)`, `genoLength_Silking (GDU)`,
         `genoLength_Plant height`, `genoLength_Ear height`, 
         `genoLength_Plot grain weight`, `genoLength_Ear length`,
         `genoLength_Ear width`, `genoLength_Kernels per row`,
         `genoLength_Kernel row number`, `genoLength_Kernel weight`,
         `genoLength_Kernel area`, `genoLength_Kernel length`,
         `genoLength_Kernel width`, `genoLength_Kernel thickness`) 

colMins <- as.data.frame(apply(X = allStability[,2:ncol(allStability)], MARGIN = 2, min, na.rm = T))
colMaxs <-  as.data.frame(apply(allStability[,2:ncol(allStability)], 2, max, na.rm = T))
colRanges <-  as.data.frame(apply(allStability[,2:ncol(allStability)], 2, range, na.rm = T))
colMeans <-  as.data.frame(apply(allStability[,2:ncol(allStability)], 2, mean, na.rm = T))
colMedians <- as.data.frame(apply(allStability[,2:ncol(allStability)], 2, median, na.rm = T))
colVars <- as.data.frame(apply(allStability[,2:ncol(allStability)], 2, var, na.rm = T))
traits <- as.data.frame(colnames(selectedPhenoData[,5:18]))

ggeSummary <- cbind(traits,colMins,colMedians,colMaxs, colVars)
colnames(ggeSummary) <- c("Trait", "Min", "Median", "Max", "Var")
rownames(ggeSummary) <- seq(1:nrow(ggeSummary))

write.csv(ggeSummary, "gge_summary.csv", row.names = F)

allStability$medianGGEstability <- as.numeric(myround(rowMedian(allStability[,2:ncol(allStability)], na.rm = T), 2))

write.csv(allStability, "gge_stability_values_all_traits.csv", row.names = F)

