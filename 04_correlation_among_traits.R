##################################################
## Project: GxE 2014-2015 Inbreds
## Script purpose: Correlation among traits
## Date: 2018-09-19
## Author: Celeste Falcon
##################################################

#### Load packages and functions ####
library(here)
library(corrplot)
library(ggcorrplot)
library(ggplot2)

#### Read in, format, and organize data ####
allBlups <- as.data.frame(read.csv("gxe20142015_allTraits_BLUPs.csv", header = T))

colnames(allBlups) <- c("Genotype", "Anthesis (GDU)", "Silking (GDU)", 
                                 "Plant height", "Ear height", "Plot grain weight", "Ear length",
                                 "Ear width", "Kernels per row", "Kernel row number", "Kernel weight", 
                                 "Kernel area", "Kernel length",  "Kernel width", "Kernel thickness")


#### Correlation among traits ####
data <- allBlups[,-1]
cormat <- cor(data, use = "p", method = "pearson")

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
p.mat <- cor.mtest(allBlups[,-1])
head(p.mat[, 1:5])

tiff("gxe20142015Inbreds-correlationPlot.tif")
corrplot(cormat, method = "color", 
         # addCoef.col = "gray40",
         tl.col = "black",  
         type = "upper", diag = F, p.mat = p.mat, 
         mar = c(0,3,0,0), insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "white")
dev.off()

ggcorrplot <- ggcorrplot(cor(allBlups[,-1]), p.mat = p.mat, type = "upper",
                         colors = c("darkred", "white", "steelblue"),
                         ggtheme = ggplot2::theme_bw, legend.title = "Correlation\nCoefficient") + 
  labs(x = NULL, y = NULL) +
  theme(text = element_text(family = "Courier"),
        axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 12, angle = 45, hjust = 1), 
        axis.title.y = element_text(size = 18), axis.text.y = element_text(size = 12),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        panel.border = element_blank(), axis.line = element_line(color = "black"))
ggsave(ggcorrplot, filename = "ggCorrPlot.tif", device = "tiff", width = 6, height = 4.5, units = "in")


cormat2 <- cormat

diag(cormat2) <- NA

inds <- arrayInd(which.min(cormat2), dim(cormat2))
rnames = rownames(cormat2)[inds[,1]]
cnames = colnames(cormat2)[inds[,2]]

inds2 <- arrayInd(which.max(cormat2), dim(cormat2))
rnames2 = rownames(cormat2)[inds2[,1]]
cnames2 = colnames(cormat2)[inds2[,2]]

posCorr <- which(cormat2 >= 0.4 & cormat2 < 1.0, arr.ind = T)
posCorr[,"row"] <- rownames(cormat2)[as.numeric(posCorr[,"row"])]
posCorr[,"col"] <- colnames(cormat2)[as.numeric(posCorr[,"col"])]

negCorr <- which(cormat <= -0.4, arr.ind = T)
negCorr[,"row"] <- rownames(cormat2)[as.numeric(negCorr[,"row"])]
negCorr[,"col"] <- colnames(cormat2)[as.numeric(negCorr[,"col"])]


