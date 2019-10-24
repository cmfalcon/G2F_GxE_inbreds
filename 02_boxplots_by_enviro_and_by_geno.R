##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: Boxplots by enviro and by geno,
##                 pheno means by geno (for Table
##                 _ of manuscript)
## Date: 2018-09-18
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(ggplot2)
library(broman)
library(ggpubr)

median_ <- function(...) median(..., na.rm = T)
mean_ <- function(...) mean(..., na.rm = T)
#### Read in, format, and organize data ####
selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)

#### Boxplots of trait data by environment ####

colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length",
                                 "Ear width", "Kernels per row", "Kernel row number", "Kernel weight", 
                                 "Kernel area", "Kernel length",  "Kernel width", "Kernel thickness")

for (i in 5:ncol(selectedPhenoData)) {
  # i = 5
  dat <- subset(selectedPhenoData, !is.na(selectedPhenoData[,i]))
  
  boxplotByEnviro <- ggplot(selectedPhenoData, mapping = aes(x = reorder(factor(Enviro), selectedPhenoData[,i], median_), selectedPhenoData[,i])) +
    geom_boxplot(aes(alpha = .9)) +
    scale_color_hue(l = 40, c = 100) +
    labs(x = "", y = colnames(dat)[i]) +
    theme_bw() +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12, angle = 70, hjust = 1),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.line = element_line(color = "black"),
          legend.position = 'none')
  ggsave(plot = boxplotByEnviro, filename = paste0("./boxplotsByEnviro/", colnames(dat[i]), ".tif"), device = "tiff",
         height = 3, width = 6, units = "in")
  
  assign(value = boxplotByEnviro, x = paste0("boxplotByEnviro_", colnames(selectedPhenoData)[i]))
}

mean_rank_enviros <- matrix(data = NA, nrow = 36, ncol = 0)
Enviro <- as.character(unique(selectedPhenoData$Enviro))
mean_rank_enviros <- cbind(mean_rank_enviros, Enviro)

for (i in 5:ncol(selectedPhenoData)) {
  # i = 7
  dat <- subset(selectedPhenoData, !is.na(selectedPhenoData[,i]))
  means <- aggregate(selectedPhenoData[,i] ~  Enviro, selectedPhenoData, mean_)
  means$rank <- transform(means, rank = ave(x = means$`selectedPhenoData[, i]`,
                                                           FUN = function(x) rank(-x, ties.method = "min")))
  
  mean_rank_enviros <- merge(mean_rank_enviros, means$rank, by = "Enviro", all.x = T)
  colnames(mean_rank_enviros)[i + (i - 7)] <- colnames(selectedPhenoData)[i]
}

enviro_boxplots_combined <- ggarrange(`boxplotByEnviro_Anthesis (GDU)`, `boxplotByEnviro_Silking (GDU)`, `boxplotByEnviro_Plant height`,
                                          `boxplotByEnviro_Ear height`, `boxplotByEnviro_Plot grain weight`, `boxplotByEnviro_Ear length`,
                                          `boxplotByEnviro_Ear width`, `boxplotByEnviro_Kernels per row`, `boxplotByEnviro_Kernel row number`,
                                          `boxplotByEnviro_Kernel weight`, `boxplotByEnviro_Kernel area`, `boxplotByEnviro_Kernel length`,
                                          `boxplotByEnviro_Kernel width`, `boxplotByEnviro_Kernel thickness`,
                                          # labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"),
                                          ncol = 4, nrow = 4)
ggsave(enviro_boxplots_combined, filename = "./enviro_boxplots_combined.png", device = "png",
       height = 15, width = 20, units = "in")

#### Calculate means for each genotype for each trait ####
phenotypicMeans <- aggregate(selectedPhenoData[,5:18], by = list(selectedPhenoData$Genotype), mean, na.rm = T)

colMins <- as.data.frame(apply(X = phenotypicMeans[,2:ncol(phenotypicMeans)], MARGIN = 2, min, na.rm = T))
colMaxs <-  as.data.frame(apply(phenotypicMeans[,2:ncol(phenotypicMeans)], 2, max, na.rm = T))
colRanges <-  as.data.frame(apply(phenotypicMeans[,2:ncol(phenotypicMeans)], 2, range, na.rm = T))
colMeans <-  as.data.frame(apply(phenotypicMeans[,2:ncol(phenotypicMeans)], 2, mean, na.rm = T))
traits <- as.data.frame(colnames(phenotypicMeans[2:ncol(phenotypicMeans)]))

phenotypicMeansSummary <- cbind(traits,colMins,colMeans,colMaxs)

colnames(phenotypicMeansSummary) <- c("Trait", "Min", "Mean", "Max")

rownames(phenotypicMeansSummary) <- seq(1:nrow(phenotypicMeansSummary))

phenotypicMeansSummary$Min <- as.numeric(myround(phenotypicMeansSummary$Min,2))
phenotypicMeansSummary$Mean <- as.numeric(myround(phenotypicMeansSummary$Mean,2))
phenotypicMeansSummary$Max <- as.numeric(myround(phenotypicMeansSummary$Max,2))

rm(colMaxs, colMeans, colMins, colRanges)

## Boxplots of trait data by genotype ####
for (i in 5:ncol(selectedPhenoData)) {
  # i = 10
  boxplotByGenotype <- ggplot(data = selectedPhenoData, mapping = aes(x = reorder(x = Genotype, X = selectedPhenoData[,i], FUN = var, na.rm = T), y = selectedPhenoData[,i])) +
    geom_boxplot() +
    labs(x = "", y = colnames(selectedPhenoData)[i]) +
    theme_bw() +
    theme(text = element_text(family = "Courier"),
          axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12, angle = 70, hjust = 1),
          axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
          panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          panel.border = element_blank(), axis.line = element_line(color = "black"),
          legend.position = 'none')
  ggsave(plot = boxplotByGenotype, filename = paste0("./boxplotsByGenotype/", colnames(dat[i]), ".png"), device = "png",
         height = 5, width = 6, units = "in")
  assign(value = boxplotByGenotype, x = paste0("boxplotByGenotype_", colnames(selectedPhenoData)[i]))
  
}

mean_rank_genotypes <- matrix(data = NA, nrow = 31, ncol = 0)

for (i in 5:ncol(selectedPhenoData)) {
  # i = 5
  dat <- subset(selectedPhenoData, !is.na(selectedPhenoData[,i]))
  means <- aggregate(selectedPhenoData[,i] ~  Genotype, selectedPhenoData, mean_)
  means$rank <- transform(means, rank = ave(x = means$`selectedPhenoData[, i]`,
                                            FUN = function(x) rank(-x, ties.method = "min")))
  
  mean_rank_genotypes <- cbind(mean_rank_genotypes, means$rank)
  colnames(mean_rank_genotypes)[3 * (i - 4)] <- colnames(selectedPhenoData)[i]
}

geno_boxplots_combined <- ggarrange(`boxplotByGenotype_Anthesis (GDU)`, `boxplotByGenotype_Silking (GDU)`, `boxplotByGenotype_Plant height`,
          `boxplotByGenotype_Ear height`, `boxplotByGenotype_Plot grain weight`, `boxplotByGenotype_Ear length`,
          `boxplotByGenotype_Ear width`, `boxplotByGenotype_Kernels per row`, `boxplotByGenotype_Kernel row number`,
          `boxplotByGenotype_Kernel weight`, `boxplotByGenotype_Kernel area`, `boxplotByGenotype_Kernel length`,
          `boxplotByGenotype_Kernel width`, `boxplotByGenotype_Kernel thickness`,
          # labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n"),
          ncol = 4, nrow = 4)
ggsave(geno_boxplots_combined, filename = "./geno_boxplots_combined.png", device = "png",
         height = 15, width = 20, units = "in")



