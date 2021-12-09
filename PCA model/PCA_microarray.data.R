
#Install BiocManager package

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maEndToEnd", version = "devel")


install.packages("devtools")
library(devtools)

devtools::install_github("r-lib/remotes")
library(remotes)
packageVersion("remotes") # has to be 1.1.1.9000 or later

remotes::install_github("b-klaus/maEndToEnd", ref="master")

#Install workflow package
suppressPackageStartupMessages({library("maEndToEnd")})

#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
#library(devtools)


raw_data_dir <- tempdir()

if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}

anno_AE <- getAE("E-MTAB-76", path = raw_data_dir, type = "raw")


sdrf_location <- file.path(raw_data_dir, "E-MTAB-76.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)


raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir, 
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))

#Subselect the columns of interest for the analysis
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[,c("Source.Name",
                                                        "Factor.Value..strain.",
                                                        "Factor.Value..treatment.",
                                                        "Factor.Value..time.")]
#Boxplot of the raw data
boxplot(raw_data, target = "core", 
        main = "Boxplot of log2-intensitites for the raw data")


#Import pre-processed data
data <- read.csv(file.choose())
dim(data)

#Assign first column to be rownames
rownames(data) <- data[,1]
data <- data[,-1]


#Perform PCA of transposed data
pca <- prcomp(t(data), center = TRUE, scale. = TRUE)
summary(pca)

#PCA scree plot
plot(pca, type="l", main= "Scree plot")

#Save loadings
loadings <- pca$rotation
write.csv(loadings, file = "loadings.csv")

#Save scores
scores <- pca$x
write.csv(scores, file= "scores.csv")

#PCA scores plot 

percentVariance <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVariance[2] / percentVariance[1])
dataframe <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                    Time = pData(raw_data)$Factor.Value..time.,
                    Treatment = pData(raw_data) $Factor.Value..treatment.)

ggplot(dataframe, aes(PC1, PC2)) +
            geom_point(aes(shape = Treatment, colour = factor(Time))) +
        ggtitle("Scores plot of PC1 and PC2") +
        xlab(paste0("PC1, VarExp: ", percentVariance[1], "%")) +
        ylab(paste0("PC2, VarExp: ", percentVariance[2], "%")) +
        theme(plot.title = element_text(hjust = 0.5))+
        coord_fixed(ratio = sd_ratio) +
        scale_shape_manual(values = c(0, 15, 4)) +
        scale_color_manual(values = c( "1"= "#0000ad", "2"= "#0000be", "3"= "#0000ce", "4"= "#0000df", "5"= "#0000f0", "6"= "#0a04ff", "7"= "#3318ff", 
                                              "8"= "#4825ff", "10"= "#572fff", "12"= "#6539ff", "14"= "#7142ff", "16"= "#7d4aff", "18"= "#8852ff", "20"= "#925aff",
                                              "22"= "#9c62ff", "23"= "#a66aff", "24"= "#af71ff", "28"= "#b879ff", "32"= "#c280ff", "36"= "#cb88ff", "40"= "#d38fff",
                                              "41"= "#dc97ff", "46"= "#dc97ff", "48"="#eea6ff", "51"="#f6adff", "52"="#feb5ff", "56"="#ffc0ff")) +
     labs(color='Time')



#new PCA model for B4053

percentVariance <- round(100*pca$sdev^2/sum(pca$sdev^2), 1)
sd_ratio <- sqrt(percentVariance[2] / percentVariance[1])
dataframe <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                     Time = pData$Factor.Value..time.,
                     Treatment = pData$Factor.Value..treatment.)
ggplot(dataframe, aes(PC1, PC3)) +
  geom_point(size=3, aes(shape = Treatment, colour = factor(Time))) +
  ggtitle("Scores plot of PC1 and PC2") +
  xlab(paste0("PC1, VarExp: ", percentVariance[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVariance[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(15, 4)) +
  scale_color_manual(values = c( "1"= "#000099", "2"= "#004C99", "3"= "#0066CC", "4"= "#0080FF", "5"= "#3399FF", "6"= "#66B2FF", "7"= "#99CCFF", 
                                 "8"= "#99CCFF")) +
  labs(color='Time') + theme(axis.text=element_text(size=12),
                             axis.title=element_text(size=12), plot.title=element_text(size=14))

