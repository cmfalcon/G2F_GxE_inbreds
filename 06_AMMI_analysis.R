##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: GxE-AMMI analysis
## Date: 2018-09-20
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(plyr)
library(agricolae)
library(ggplot2)
library(dplyr)
library(data.table)
library(broman)

get_stars = function(p) {
  stars = findInterval(p, c(0, 0.001, 0.01, 0.05, 0.1))
  codes = c("***" , "**","*", ".", " ")
  codes[stars]
}
rowMedian <- function(x, na.rm = FALSE) apply(x, 1, median, na.rm = na.rm)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

#### Read in, format, and organize data ####
selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")

varianceComponents <- read.csv("variance_components.csv", header = T)

#### AMMI biplot analysis ####
# dir.create("./AMMIbiplots/standardized")

varComponents <- matrix(, nrow = 5, ncol = 0)
significance <- matrix(, nrow = 5, ncol = 0)
rownames(significance) <- (c("Loc", "Genotype", "LocRep", "LocGenotype", "Residuals"))
varExplainedByModel <- matrix(, nrow = 0, ncol = 1)

enviros <- sort(unique(selectedPhenoData$Enviro), decreasing = F)
genos <- sort(unique(selectedPhenoData$Genotype), decreasing = F)

ASVrank_allTraits <- matrix(genos, nrow = length(unique(selectedPhenoData$Genotype)), ncol = 1)
colnames(ASVrank_allTraits) <- "Genotype"

meanEuclideanDist_allTraits <- matrix(, nrow = 0, ncol = 2)
medianEuclideanDist_allTraits <- matrix(, nrow = 0, ncol = 2)

AMMI_PC_varExpl <- matrix(NA, nrow = 0, ncol = 4)
colnames(AMMI_PC_varExpl) <- c("Trait", "PC1", "PC2", "PC1+PC2")

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


for (i in 5:ncol(selectedPhenoData_standardized)) {
  # i = 5
  trait = selectedPhenoData_standardized[,i]
  df <- selectedPhenoData_standardized[!is.na(trait),]
  df <- droplevels(df)
  
  df2 <- daply(.data = df,
               .variables = c("Enviro", "Genotype"),
               .fun = function(df) mean(df[,i]))
  df2 <- t(df2)
  
  ngen <- nrow(df2)
  nenv <- ncol(df2)
  
  enviros <- sort(unique(df$Enviro), decreasing = F)
  genos <- sort(unique(df$Genotype), decreasing = F)
  
  model <- AMMI(ENV = df$Enviro, GEN = df$Genotype, REP = df$Rep, Y = df[,i])
  
  ## % explained by PC1, PC2, and PC1+PC2
  PC1_varExpl <- model$analysis$percent[1]
  PC2_varExpl <- model$analysis$percent[2]
  sumPC1PC2_varExpl <- sum(PC1_varExpl, PC2_varExpl)
  
  row <- c(colnames(selectedPhenoData_standardized)[i], PC1_varExpl, PC2_varExpl, sumPC1PC2_varExpl)
  
  AMMI_PC_varExpl <- rbind(AMMI_PC_varExpl, row)
  
  anovaTable <- model$ANOVA
  
  sig = get_stars(anovaTable$`Pr(>F)`)
  significance <- cbind(significance, sig)
  
  pdf(paste0("./AMMIbiplots/standardized/gxe20142015Inbreds-biplotPC1PC2-", colnames(df[i]), ".pdf"))
  biplot <- plot(model)
  AMMI.contour(model,distance = 0.7, shape = 8, col = "red", lwd = 2, lty = 5)
  dev.off()
  
  pdf(paste0("./AMMIbiplots/standardized/gxe20142015Inbreds-biplotPC1vs", colnames(df[i]), ".pdf"))
  # biplot PC1 vs plotWeight
  biplot <- plot(model,0,1,angle = 20,ecol = "brown", xlab = colnames(df[i]))
  dev.off()
  
  #### Mean Euclidean distance among points ####
  biplotDataPoints <- model$biplot[,3:4]
  
  ## calculate Euclidean distance for genotypes
  biplotDataPoints_geno <- biplotDataPoints[1:ngen,]
  meanEuclideanDist_geno <- mean(dist(biplotDataPoints_geno, method = "euclidean"))
  medianEuclideanDist_geno <- median(dist(biplotDataPoints_geno, method = "euclidean"))
  
  ## calculate Euclidean distance for environments
  biplotDataPoints_enviro <- biplotDataPoints[(ngen + 1):(ngen + nenv),]
  meanEuclideanDist_enviro <- mean(dist(biplotDataPoints_enviro, method = "euclidean"))
  medianEuclideanDist_enviro <- median(dist(biplotDataPoints_enviro, method = "euclidean"))
  
  meanEuclideanDist <- cbind(meanEuclideanDist_geno, meanEuclideanDist_enviro)
  meanEuclideanDist_allTraits <- rbind(meanEuclideanDist_allTraits, meanEuclideanDist)
  
  medianEuclideanDist <- cbind(medianEuclideanDist_geno, medianEuclideanDist_enviro)
  medianEuclideanDist_allTraits <- rbind(medianEuclideanDist_allTraits, medianEuclideanDist)
  
  #### AMMI Stability Value ####
  ASVvalues <- matrix(nrow = ngen, ncol = 1, dimnames = list(rownames(df2), "ASVvalue"))
  for (j in 1:ngen) {
    # j = 1
    IPCA1SS <- model$analysis[1,4]
    IPCA2SS <- model$analysis[2,4]
    
    IPCA1score <- model$biplot[j,3]
    IPCA2score <- model$biplot[j,4]
    
    ASVvalues[j,1] <- sqrt(((IPCA1SS/IPCA2SS) * IPCA1score)^2 + (IPCA2score)^2)
  }
  
  rankASVvalues <- transform(ASVvalues, rank = ave(ASVvalues,
                                                   FUN = function(x) rank(x, ties.method = "min")))
  rankASVvalues <- setDT(rankASVvalues, keep.rownames = T)[]
  colnames(rankASVvalues)[1] <- "Genotype"
  colnames(rankASVvalues)[2] <- paste0("ASV_", colnames(selectedPhenoData_standardized[i]))
  colnames(rankASVvalues)[3] <- paste0("rank_", colnames(selectedPhenoData_standardized[i]))
  
  ASVrank_allTraits <- merge(ASVrank_allTraits, rankASVvalues, by = "Genotype", all = T)
  
  source("/Users/celestefalcon/Box Sync/Celeste Falcon/GxE/stability.par.cfalcon.R")
  YSi <- stability.par.cfalcon(df2, rep = 2, MSerror = 4.405999, alpha = 0.05, main = "Genotype", console = T)
}

ASVrank_allTraits2 <- ASVrank_allTraits[,seq(1, ncol(ASVrank_allTraits), 2)]
ASVrank_allTraits2$medianRank <- as.numeric(myround(rowMedian(ASVrank_allTraits2[,2:ncol(ASVrank_allTraits2)], na.rm = T),2))

ASV_value_allTraits <- ASVrank_allTraits[,c(1, seq(2, ncol(ASVrank_allTraits), 2))]
ASV_value_allTraits[ASV_value_allTraits == "NaN"] <- NA_character_
traits <- colnames(selectedPhenoData[5:ncol(selectedPhenoData)])
colnames(ASV_value_allTraits) <- c("Genotype", traits)

write.csv(ASV_value_allTraits, "ASV_value_all_traits.csv", row.names = F)

ASV_value_allTraits_withMed <- ASV_value_allTraits
ASV_value_allTraits_withMed$median <- as.numeric(rowMedian(ASV_value_allTraits_withMed[,2:ncol(ASV_value_allTraits_withMed)], na.rm = T))

write.csv(ASV_value_allTraits_withMed, "ASV_value_all_traits_with_median.csv", row.names = F)

ASV_value_allTraits_melt <- melt(ASV_value_allTraits)

rownames(meanEuclideanDist_allTraits) <- colnames(selectedPhenoData[,5:ncol(selectedPhenoData)])
rownames(medianEuclideanDist_allTraits) <- colnames(selectedPhenoData[,5:ncol(selectedPhenoData)])

meanEuclideanDist <- as.data.frame(meanEuclideanDist_allTraits[1:14,])
meanEuclideanDist <- setDT(meanEuclideanDist, keep.rownames = T)
colnames(meanEuclideanDist) <- c("Trait", "Genotype", "Environment")
meanEuclideanDist$Trait <- factor(meanEuclideanDist$Trait, levels = unique(meanEuclideanDist$Trait))

meanEuclideanDist <- melt(meanEuclideanDist, id.vars = "Trait")

G_E_Total_euclideanDistance <- as.data.frame(medianEuclideanDist_allTraits)
G_E_Total_euclideanDistance$medianEuclideanDist_geno <- as.numeric(myround(G_E_Total_euclideanDistance$medianEuclideanDist_geno, 2))
G_E_Total_euclideanDistance$medianEuclideanDist_enviro <- as.numeric(myround(G_E_Total_euclideanDistance$medianEuclideanDist_enviro, 2))

G_E_Total_euclideanDistance$Total <- G_E_Total_euclideanDistance$medianEuclideanDist_geno + G_E_Total_euclideanDistance$medianEuclideanDist_enviro
colnames(G_E_Total_euclideanDistance) <- c("Genotypes", "Environments", "Total")

write.csv(G_E_Total_euclideanDistance, "total_euc_distances_from_AMMI.csv", row.names = F)

G_E_Total_euclideanDistance_2 <- setDT(G_E_Total_euclideanDistance, keep.rownames = T)
colnames(G_E_Total_euclideanDistance_2)[1] <- "Trait"

euclideanDistanceTotal <-
  ggplot(data = G_E_Total_euclideanDistance_2, mapping = aes(x = Trait, y = Total)) +
  geom_col(fill = "darkred") +
  geom_col(mapping = aes(x = Trait, y = Environments), fill = "steelblue") + 
  scale_x_discrete(limits = G_E_Total_euclideanDistance_2$Trait) +
  labs(x = "Trait", y = "Euclidean Distance") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16, angle = 50, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = euclideanDistanceTotal, filename = "medianEuclideanDistanceTotal_AMMIplot.tif", device = "tiff", 
       width = 10, height = 5, units = "in")


euclideanDistanceTotal_ordered <-
  ggplot(data = G_E_Total_euclideanDistance_2, mapping = aes(x = reorder(x = G_E_Total_euclideanDistance_2$Trait, X = G_E_Total_euclideanDistance_2$Total), y = G_E_Total_euclideanDistance_2$Total)) +
  geom_col(fill = "darkred") +
  geom_col(mapping = aes(x = reorder(x = Trait, X = Total, FUN = mean), y = Environments), fill = "steelblue") +
  labs(x = "Trait", y = "Euclidean Distance") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16, angle = 60, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 16),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(plot = euclideanDistanceTotal_ordered, filename = "medianEuclideanDistanceTotal_AMMIplot_ordered.tif", device = "tiff",
       width = 10, height = 5, units = "in")

#### Correlation of Euclidean distances and percent variance explained ####
varianceComponents_transposed <- t(varianceComponents)
colnames(varianceComponents_transposed) <- as.character(varianceComponents_transposed[1,])
varianceComponents_transposed <- varianceComponents_transposed[-1,]
varianceComponents_transposed <- as.data.frame(varianceComponents_transposed)

varianceComponents_transposed[,1] <- as.numeric.factor(varianceComponents_transposed[,1])
varianceComponents_transposed[,2] <- as.numeric.factor(varianceComponents_transposed[,2])
varianceComponents_transposed[,3] <- as.numeric.factor(varianceComponents_transposed[,3])
varianceComponents_transposed[,4] <- as.numeric.factor(varianceComponents_transposed[,4])
varianceComponents_transposed[,5] <- as.numeric.factor(varianceComponents_transposed[,5])

corr_percVarGE_EucTotal <- cor(varianceComponents_transposed$`Genotype x
Enviro`, G_E_Total_euclideanDistance$Total, method = "pearson")

corr_percVarG_EucTotal <- cor(varianceComponents_transposed$Genotype, G_E_Total_euclideanDistance$Total, method = "pearson")
corr_percVarE_EucTotal <- cor(varianceComponents_transposed$Enviro, G_E_Total_euclideanDistance$Total, method = "pearson")

corr_percVarG_EucG <- cor(varianceComponents_transposed$Genotype, G_E_Total_euclideanDistance$Genotypes, method = "pearson")
corr_percVarE_EucE <- cor(varianceComponents_transposed$Enviro, G_E_Total_euclideanDistance$Environments, method = "pearson")

