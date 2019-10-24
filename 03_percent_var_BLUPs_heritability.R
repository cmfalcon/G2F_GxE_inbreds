##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: Percent variation explained,
##                 BLUPs, and heritability
## Date: 2018-09-19
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(data.table)
library(lme4)
library(broman)
library(ggplot2)
library(lmerTest)

specify_decimal <- function(x, k) format(round(x, k), nsmall = k)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
get_stars = function(p) {
  stars = findInterval(p, c(0, 0.001, 0.01, 0.05, 0.1))
  codes = c("***" , "**","*", ".", "ns")
  codes[stars]
}

#### Read in, format, and organize data ####
selectedPhenoData <- read.csv("./selectedPhenoData.csv", header = T)
colnames(selectedPhenoData) <- c("Plot", "Enviro", "Rep", "Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                 "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                 "Kernel length", "Kernel width", "Kernel thickness")


#### Calculate percent variation explained, BLUPs, and heritability ####

### Get variance components
varComponents <- matrix(0, nrow = 5, ncol = 0)
heritability <- matrix(data = NA, nrow = 0, ncol = 2)
allBlups <- as.data.frame(unique(selectedPhenoData$Genotype))
setDT(allBlups, keep.rownames = T)[]
allBlups <- allBlups[,-1]
colnames(allBlups) <- "Genotype"

significance <- matrix(,nrow=6, ncol=0)
rownames(significance) <- (c("none", "Enviro:Rep", "Rep", "Enviro", "Genotype", "Genotype:Enviro"))

for (i in 1:14) {
  # i = 1
  
  trait = selectedPhenoData[,i + 4]
  
  df <- selectedPhenoData[!is.na(trait),]
  df <- droplevels(df)
  
  ## 1| denotes random effect
  tmpge <- lmer(trait ~ (1|Rep/Enviro) + (1|Enviro) + (1|Genotype) + (1|Genotype:Enviro), data = selectedPhenoData, na.action = na.omit)
  summ <- summary(tmpge)
  
  ranova_tmpge <- ranova(tmpge)
  
  genVar = summ$varcor$Genotype
  genEnviroVar = summ$varcor$`Genotype:Enviro`
  enviroVar = summ$varcor$Enviro
  repVar = summ$varcor$Rep
  repEnviroVar = summ$varcor$`Enviro:Rep`
  resVar = (summ$sigma)^2
  
  totalVar = sum(genVar,genEnviroVar,enviroVar,repVar,repEnviroVar,resVar)
  
  nEnv <- length(unique(df$Enviro))
  
  herit = genVar/(genVar + (genEnviroVar/nEnv) + (resVar/(2*nEnv)))
  herit2 <- c(colnames(selectedPhenoData)[i + 4], herit)
  heritability <- rbind(heritability, herit2)
  
  genVarProportion = genVar/totalVar * 100
  genEnviroVarProportion = genEnviroVar/totalVar * 100
  enviroVarProportion = enviroVar/totalVar  * 100
  resVarProportion = resVar/totalVar * 100
  repEnviroVarProportion = repEnviroVar/totalVar * 100
  
  sig <- get_stars(ranova_tmpge$`Pr(>Chisq)`)
  significance <- cbind(significance, sig)
  
  varComponents2 <- as.data.frame(c(enviroVarProportion, genVarProportion, genEnviroVarProportion, repEnviroVarProportion, resVarProportion))
  varComponents <- cbind(varComponents, varComponents2)
  
  ### This data frame contains the BLUPs:
  lmBlups = ranef(tmpge)$'Genotype'
  
  add1 = mean(selectedPhenoData[,i + 4], na.rm = TRUE)
  
  blupAdd = lmBlups[,1] + add1
  lmBlups = data.frame(lmBlups,blupAdd)
  setDT(lmBlups, keep.rownames = T)
  lmBlups <- lmBlups[,c(1,3)]
  colnames(lmBlups) = c("Genotype", paste(colnames(selectedPhenoData)[i + 4],"Blup",sep = "_"))
  
  allBlups = merge(allBlups, lmBlups, by = "Genotype", all = T)
  
  rm(tmpge,lmBlups,add1,trait)
}

colnames(allBlups) <- c("Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                        "Kernel length", "Kernel width", "Kernel thickness")
colnames(heritability) <- c("Trait", "Heritability")
heritability <- as.data.frame(heritability)
heritability$Heritability <- as.numeric.factor(heritability$Heritability)
heritability[,2] <- myround(as.numeric(heritability[,2]), 3)

#### Write .csv file of BLUPs and heritabilities
write.csv(allBlups,"gxe20142015_allTraits_BLUPs.csv",row.names = F,quote = F)
write.csv(heritability, "gxe20142015_allTraits_Heritability.csv", row.names = F, quote = F)

### Gather together variance components and significance for making plots
varComponents1 <- as.data.frame(c("Enviro", "Genotype", "Genotype x\nEnviro", "Rep/Enviro", "Residual"))
varianceComponents <- cbind(varComponents1, varComponents)
colnames(varianceComponents) <- c("Variance Component", 
                                  "Anthesis", "Silking", "Plant\nheight", "Ear\nheight",
                                  "Plot\nweight", "Ear\nlength", "Ear\nwidth",
                                  "Kernels\nper\nrow", "Kernel\nrow\nnumber", "Kernel\nweight", "Kernel\narea",
                                  "Kernel\nlength", "Kernel\nwidth", "Kernel\nthickness")
write.csv(varianceComponents, "variance_components.csv", row.names = F)

meltedVarComp <- melt(varianceComponents, id = "Variance Component")
meltedVarComp$value <- as.numeric(meltedVarComp$value)

colnames(significance) <- c("Anthesis", "Silking", "Plant\nheight", "Ear\nheight",
                            "Plot\nweight", "Ear\nlength", "Ear\nwidth",
                            "Kernels\nper\nrow", "Kernel\nrow\nnumber", "Kernel\nweight", "Kernel\narea",
                            "Kernel\nlength", "Kernel\nwidth", "Kernel\nthickness")
# meltedSig <- melt(significance, id = )

### Make bubble plots
plot <-
  ggplot(meltedVarComp, aes(x = variable, y = meltedVarComp[,1])) +
  geom_point(aes(size = value, colour = meltedVarComp[,1])) + 
  geom_text(vjust = 2.8, size = 6, label = specify_decimal(x = meltedVarComp$value,k = 2)) +
  # geom_text(vjust = 2, size = 6, label = meltedSig$value) +
  labs(x = "Trait", y = "") +
  scale_size(range = c(1,40)) +
  scale_color_hue(l = 30, c = 90) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')
ggsave(plot = plot, filename = "bubblePlot-varComps-gxe20142015Inbreds_withReps.tif", device = "tiff",
       height = 7, width = 20, units = "in")

varianceComponents_ordered <- varianceComponents[,c(1, 10, 2:4, 7:8, 5, 14, 6, 12, 11, 9, 15, 13)]
meltedVarComp_ordered <- melt(varianceComponents_ordered, id = "Variance Component")
meltedVarComp_ordered$value <- as.numeric(meltedVarComp_ordered$value)

plot_ordered <-
  ggplot(meltedVarComp_ordered, aes(x = variable, y = meltedVarComp_ordered[,1])) +
  geom_point(aes(size = value, colour = meltedVarComp_ordered[,1])) + 
  geom_text(vjust = 2.8, size = 6, label = specify_decimal(x = meltedVarComp_ordered$value,k = 2)) +
  # geom_text(vjust = 2, size = 6, label = meltedSig$value) +
  labs(x = "Trait", y = "") +
  scale_size(range = c(1,40)) +
  scale_color_hue(l = 30, c = 90) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 24), axis.text.x = element_text(size = 18), 
        axis.title.y = element_text(size = 24), axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')
ggsave(plot = plot_ordered, filename = "bubblePlot-varComps-gxe20142015Inbreds_withReps_ordered.tif", device = "tiff",
       height = 7, width = 20, units = "in")

colnames(varianceComponents) <- c("Variance Component", 
                                  "Anthesis", "Silking", 
                                  "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                                  "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                                  "Kernel length", "Kernel width", "Kernel thickness")
meltedVarComp <- melt(varianceComponents, id = "Variance Component")
meltedVarComp$value <- as.numeric(meltedVarComp$value)

#### Make panel bar plot
plot2 <-
  ggplot(meltedVarComp, aes(variable, value)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  facet_grid(~ meltedVarComp[,"Variance Component"]) +
  scale_x_discrete(limits = rev(levels(meltedVarComp$variable))) +
  labs(x = "Trait", y = "Percent phenotypic variance explained") +
  theme_bw() +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14), 
        axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 16), strip.background = element_rect(color = "white", fill = "white"),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), axis.line = element_line(color = "black"),
        legend.position = 'none')
ggsave(plot = plot2, filename = "panelBarPlot-varComps-gxe20142015Inbreds_withReps.tif", device = "tiff",
       width = 10, height = 3.5, units = "in")

rm(df, enviroVar, enviroVarProportion, genEnviroVar, genEnviroVarProportion, genVar, 
   genVarProportion, herit, repEnviroVar, repEnviroVarProportion, repVar, 
   varComponents, varComponents1, varComponents2, blupAdd, herit2, i, nEnv, plot, resVar,
   resVarProportion, summ, totalVar, traits)

#### BLUPs Summary ####
allBlups <- as.data.frame(read.csv("gxe20142015_allTraits_BLUPs.csv", header = T))
colnames(allBlups) <- c("Genotype", 
                        "Anthesis", "Silking", 
                        "Plant height", "Ear height", "Plot grain weight", "Ear length", "Ear width",
                        "Kernels per row", "Kernel row number", "Kernel weight", "Kernel area",
                        "Kernel length", "Kernel width", "Kernel thickness")

colMins <- as.data.frame(apply(X = allBlups[,2:ncol(allBlups)], MARGIN = 2, min, na.rm = T))
colMaxs <-  as.data.frame(apply(allBlups[,2:ncol(allBlups)], 2, max, na.rm = T))
colRanges <-  as.data.frame(apply(allBlups[,2:ncol(allBlups)], 2, range, na.rm = T))
colMeans <-  as.data.frame(apply(allBlups[,2:ncol(allBlups)], 2, mean, na.rm = T))
traits <- as.data.frame(colnames(allBlups[2:ncol(allBlups)]))

blupsSummary <- cbind(traits,colMins,colMeans,colMaxs)

colnames(blupsSummary) <- c("Trait", "Min", "Mean", "Max")

rownames(blupsSummary) <- seq(1:nrow(blupsSummary))

blupsSummary$Min <- as.numeric(myround(blupsSummary$Min,2))
blupsSummary$Mean <- as.numeric(myround(blupsSummary$Mean,2))
blupsSummary$Max <- as.numeric(myround(blupsSummary$Max,2))

rm(colMaxs, colMeans, colMins, colRanges)


#### With location and year separated ####
selectedPhenoData_loc_year <- selectedPhenoData
selectedPhenoData_loc_year$Loc <- as.factor(gsub("_.*", "", selectedPhenoData_loc_year$Enviro))
selectedPhenoData_loc_year$Year <- as.factor(gsub(".*_", "", selectedPhenoData_loc_year$Enviro))

### Get variance components
varComponents_loc_year <- matrix(0, nrow = 8, ncol = 0)

significance <- matrix(,nrow=7, ncol=0)
rownames(significance) <- (c("none", "Rep:Loc:Year", "Year", "Loc", "Genotype", "Genotype:Year", "Genoype:Loc"))

for (i in 1:14) {
  # i = 1
  
  trait = selectedPhenoData[,i + 4]
  
  df <- selectedPhenoData[!is.na(trait),]
  df <- droplevels(df)
  
  ## 1| denotes random effect
  tmpge <- lmer(trait ~ (1|Rep/Loc:Year) + (1|Loc) + (1|Year) + (1|Genotype) + (1|Genotype:Loc) + (1|Genotype:Year), data = selectedPhenoData_loc_year, na.action = na.omit)
  summ <- summary(tmpge)
  
  ranova_tmpge <- ranova(tmpge)
  
  genVar = summ$varcor$Genotype
  genLocVar = summ$varcor$`Genotype:Loc`
  genYearVar = summ$varcor$`Genotype:Year`
  locVar = summ$varcor$Loc
  yearVar = summ$varcor$Year
  repLocYearVar = summ$varcor$`Loc:Year:Rep`
  repVar = summ$varcor$Rep
  resVar = (summ$sigma)^2
  
  totalVar = sum(genVar,genLocVar,genYearVar,locVar,yearVar,repLocYearVar,repVar,resVar)
  
  genVarProportion = genVar/totalVar * 100
  genLocVarProportion = genLocVar/totalVar * 100
  genYearVarProportion = genYearVar/totalVar * 100
  locVarProportion = locVar/totalVar * 100
  yearVarProportion = yearVar/totalVar * 100
  repLocYearVarProportion = repLocYearVar/totalVar * 100
  repVarProportion = repVar/totalVar * 100
  resVarProportion = resVar/totalVar * 100
  
  sig <- get_stars(ranova_tmpge$`Pr(>Chisq)`)
  significance <- cbind(significance, sig)

  varComponents_loc_year2 <- as.data.frame(c(locVarProportion, yearVarProportion, genVarProportion, genLocVarProportion, genYearVarProportion, repVarProportion, repLocYearVarProportion, resVarProportion))
  varComponents_loc_year <- cbind(varComponents_loc_year, varComponents_loc_year2)
}

  varComponents_loc_year1 <- as.data.frame(c("Loc", "Year", "Genotype", "Genotype x Loc", "Genotype x Year", "Rep", "Rep/Loc*Year", "Residual"))
  varianceComponents_loc_year <- cbind(varComponents_loc_year1, varComponents_loc_year)
  colnames(varianceComponents_loc_year) <- c("Variance Component", 
                                    "Anthesis", "Silking", "Plant\nheight", "Ear\nheight",
                                    "Plot\nweight", "Ear\nlength", "Ear\nwidth",
                                    "Kernels\nper\nrow", "Kernel\nrow\nnumber", "Kernel\nweight", "Kernel\narea",
                                    "Kernel\nlength", "Kernel\nwidth", "Kernel\nthickness")  
  
  