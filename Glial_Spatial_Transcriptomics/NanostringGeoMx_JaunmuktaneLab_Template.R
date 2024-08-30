
#NANOSTRING_DATA_LOAD##################################################################################

#Load Libraries

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(gridExtra)
library(standR)
library(tidyverse)
library(SpatialExperiment)
library(svDialogs)
library(htmlwidgets)

#Set Work Directory#
setwd("...")
datadir <- "..."

list.files(datadir)


#Load DCC files
DCCFiles <- dir(file.path(datadir, "GFAP_P2RY12_DCC"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
#Load pkc files (nanostring website)
PKCFiles <- unzip(zipfile = dir(file.path(datadir, "pkcs"), pattern = ".zip$",
                                full.names = TRUE, recursive = TRUE))
#Load annotation file/ Lab Worksheet
BULKAnnotationFile <- "....xlsx"

#Compile data into a GeoMxSet

Data <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = BULKAnnotationFile,
  phenoDataSheet = "Sheet1",
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = c("aoi", "roi"),
  experimentDataColNames = c("panel")
)

#Load GO-Term annotation Database
library(msigdb)
library(GSEABase)

msigdb_hs <- getMsigdb(version = '7.2')
msigdb_hs <- appendKEGG(msigdb_hs)


#RESULTSDIRECTORY_SET##################################################################################
#SET SUBDIRECTORY FOR RESULTS#

RESULTSDIR <- paste0(datadir,"/RESULTS_APOE/")
dir.create(RESULTSDIR)
setwd(RESULTSDIR)





#CELLDIR_SET#####################################################################################################################
#SET SUBDIRECTORY FOR CellType#
CELLTYPEDLG <-  dlgInput("ENTER CELL TYPE (GFAP/P2RY12)")
if (CELLTYPEDLG$res == "GFAP"){ CELLTYPEVAL <- "GFAP"}
if (CELLTYPEDLG$res == "P2RY12"){ CELLTYPEVAL <- "P2RY12"}

celldemoData <- subset(Data, select = phenoData(Data)[["CellType"]] %in% c(CELLTYPEVAL))

CELLDIR <- paste0(RESULTSDIR,"/", CELLTYPEVAL)
dir.create(CELLDIR)
setwd(CELLDIR)




#REGIONDIR_SET#####################################################################################################################
#SET SUBDIRECTORY FOR REGION#


GM <- c("D", "M", "S")

REGIONDLG <-  dlgInput("ENTER REGION (D/M/S/W/GM)")
if (REGIONDLG$res == "D"){ REGIONVAL <- "D"}
if (REGIONDLG$res == "M"){ REGIONVAL <- "M"}
if (REGIONDLG$res == "S"){ REGIONVAL <- "S"}
if (REGIONDLG$res == "W"){ REGIONVAL <- "W"}
if (REGIONDLG$res == "GM"){ REGIONVAL <- GM}


demoData <- subset(celldemoData, select = phenoData(celldemoData)[["region"]] %in% c(REGIONVAL))

if (REGIONDLG$res == "GM"){ REGIONVAL <- "GM"}

REGIONDIR <- paste0(CELLDIR,"/REGION_", REGIONVAL)
dir.create(REGIONDIR)
setwd(REGIONDIR)


library(knitr)
pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

#NANOSTRING_QC###############


#NANOSTRING_QC_1.Sample Overview
library(dplyr)
library(ggforce)
library(networkD3)
sankeyCols <- c("source", "target", "value")
link1 <- dplyr::count(pData(demoData), 'slide name', class) 
link2 <- dplyr::count(pData(demoData), class,pathology) 
link3 <- dplyr::count(pData(demoData),  pathology, region)
colnames(link1) <- sankeyCols
colnames(link2) <- sankeyCols
colnames(link3) <- sankeyCols
links <- rbind(link1,link2,link3)
nodes <- unique(data.frame(name=c(links$source, links$target)))
# sankeyNetwork is 0 based, not 1 based
links$source <- as.integer(match(links$source,nodes$name)-1)
links$target <- as.integer(match(links$target,nodes$name)-1)
sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30)
S_NET <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                       Target = "target", Value = "value", NodeID = "name",
                       units = "TWh", fontSize = 12, nodeWidth = 30)
saveWidget(S_NET, 'NETWORK_SANKEY.html')

#NANOSTRING_QC_2. Shift counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)


#NANOSTRING_QC_3. Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more detail

QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 30, # Minimum sequencing saturation (30%)
       minArea = 100, 
       minNegativeCount = 1)          # Minimum segment area (100)

demoData <- setSegmentQCFlags(demoData, 
                              qcCutoffs = QC_params)        

#NANOSTRING_QC_4. Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")

library(ggplot2)

col_by <- "LBD"

#NANOSTRING_QC_5. Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

Trimmed <- QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
Trimmed
png (file="trimmed.png", width=1000, height=1000)
Trimmed
dev.off()

Stitched <- QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
Stitched
png (file="Stitched.png", width=1000, height=1000)
Stitched
dev.off()

Aligned <- QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
Aligned
png (file="Aligned.png", width=1000, height=1000)
Aligned
dev.off()

Saturated <- QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
Saturated
png (file="Saturated.png", width=1000, height=1000)
Saturated
dev.off()

Area <- QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
Area
png (file="Area.png", width=1000, height=1000)
Area
dev.off()

QC_histogram(sData(demoData), "nuclei", col_by, 20)

#NANOSTRING_QC_6. calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

#NANOSTRING_QC_7. explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
NegGeoMeanQC <- for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}
png (file="NegGeoMeanQC.png", width=1000, height=1000)
plt
dev.off()


#NANOSTRING_QC_8. detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:

kable(table(NTC_Count = sData(demoData)$NTC),
col.names = c("NTC Count", "# of Segments"))
kable(QC_Summary, caption = "QC Summary Table for each Segment")

#NANOSTRING_QC_10. QC SUMMARY
library(gridExtra)
library(grid)

g <- tableGrob(QC_Summary)
grid.newpage()
grid.draw(g)
png (file="SUMMARYQC.png", width=1000, height=1000)
grid.draw(g)
dev.off()

demoData <- demoData[, QCResults$QCStatus == "PASS"]

#NANOSTRING_QC_11. Subsetting our dataset has removed samples which did not pass QC
dim(demoData)

demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

#NANOSTRING_QC_12. Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
qc_df
ProbeQCPassed <- 
  subset(demoData, 
         fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
demoData <- ProbeQCPassed 

length(unique(featureData(demoData)[["TargetName"]]))

target_demoData <- aggregateCounts(demoData)

#NANOSTRING_QC_13. Define Limit Of Quantification (LOQ) SD threshold and minimum value
cutoff <- 1
minLOQ <- 1

#NANOSTRING_QC_14. Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ

#15.Filetering
#After determining the limit of quantification (LOQ) per segment, 
#we recommend filtering out either segments and/or genes with abnormally low signal. 
#Filtering is an important step to focus on the true biological data of interest.

#We determine the number of genes detected in each segment across the dataset

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_demoData)$Module == module
  Mat_i <- t(esApply(target_demoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_demoData)$TargetName, ]


#NANOSTRING_QC_16.Segment Gene Detection
#Filter out segments with exceptionally low signal.

#NANOSTRING_QC_17. Save detection rate information to pheno data
pData(target_demoData)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_demoData)$GeneDetectionRate <-
  pData(target_demoData)$GenesDetected / nrow(target_demoData)

#NANOSTRING_QC_18. Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_demoData)$DetectionThreshold <- 
  cut(pData(target_demoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

#NANOSTRING_QC_19. stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_demoData),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = class)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "region, #",
       fill = "Segment Type")
LOQPLOT <- ggplot(pData(target_demoData),
                  aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = class)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "region, #",
       fill = "Segment Type")
png (file="LOQPLOT.png", width=1000, height=1000)
LOQPLOT
dev.off()

kable(table(pData(target_demoData)$DetectionThreshold,
            pData(target_demoData)$LBD))

#NANOSTRING_QC_20.Gene Detection Rate
target_demoData <-
  target_demoData[, pData(target_demoData)$GeneDetectionRate >= .02]

dim(target_demoData)



#NANOSTRING_QC_21. Gene Filtering
library(scales) 

#NANOSTRING_QC_22. Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_demoData)]
fData(target_demoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_demoData)$DetectionRate <-
  fData(target_demoData)$DetectedSegments / nrow(pData(target_demoData))

#NANOSTRING_QC_23. Plot detection rate:

plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_demoData)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_demoData))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

GENEDETECTPLOT <- ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

png (file="GENE_DETECT_PLOT.png", width=1000, height=1000)
GENEDETECTPLOT
dev.off()









#FILTVAL#########################################################################################
# Set Gene Filter Value
FILTDLG <-  dlgInput("ENTER FILTERING VALUE % (10/20/30/50)")


if (FILTDLG$res == "10"){ FILTVAL <- 0.1}
if (FILTDLG$res == "20"){ FILTVAL <- 0.2}
if (FILTDLG$res == "30"){ FILTVAL <- 0.3}
if (FILTDLG$res == "50"){ FILTVAL <- 0.5}

FILTDIR <- paste0(REGIONDIR, "/FILTER", FILTVAL*100)
dir.create(FILTDIR)
setwd(FILTDIR)
# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use

negativeProbefData <- subset(fData(target_demoData), CodeClass == "Negative")

neg_probes <- unique(negativeProbefData$TargetName)

FILT_demoData <- 
  target_demoData[fData(target_demoData)$DetectionRate >= FILTVAL |
                    fData(target_demoData)$TargetName %in% neg_probes, ]

dim(FILT_demoData)

Dim <- as.data.frame(dim(FILT_demoData))
grid.newpage()
grid.table(Dim)
png (file="DIMENSIONS.png", width=1000, height=1000)
grid.table(Dim)
dev.off()

pData(FILT_demoData)$LOQ <- NULL




#SPE_RUV##################################################################################
#SPE_RUV_1. Convert data to SpatialExperiment
library(SpatialExperiment)

spe <- as.SpatialExperiment(FILT_demoData, normData = "exprs", forceRaw = TRUE, )

library(standR)
library(SpatialExperiment)
assayNames(spe)
assay(spe, "counts") <- assay(spe, "GeoMx")
assay(spe, "counts")[1:5,1:5]
assay(spe, "GeoMx") <- NULL
assay(spe, "logcounts") <- log(assay(spe, "counts"),2)
assay(spe, "logcounts")[1:5,1:5]
assayNames(spe)
colData(spe)[1:5,1:5]

rowData(spe)[1:5,1:5]

metadata(spe)$NegProbes <- negativeGeoMeans
metadata(spe)$NegProbes
metadata(spe)$QCMetrics

#SPE_RUV_2. Raw Data Plots

library(ggplot2)
library(ggalluvial)

SampleInfo <- plotSampleInfo(spe, column2plot = c("class", "region", "pathology", "scan_name"))
SampleInfo
png (file="SampleInfo.png", width=1000, height=1000)
SampleInfo
dev.off()


spe <- addPerROIQC(spe, min_count=3)
dim(spe)

metadata(spe) |> names()               
plotGeneQC(spe, ordannots = "class", col = APOE, point_size = 2, top_n=12)

plotROIQC <- plotROIQC(spe,  x_axis = "area", x_lab = "AreaSize", y_axis = "lib_size", y_lab = "Library size", col = APOE)
plotROIQC
png (file="plotROIQC.png", width=1000, height=1000)
plotROIQC
dev.off()


plotRLExpr(spe)
plotRLExpr(spe, ordannots = "segment", assay = 2, color = APOE)
RLE_spe<- plotRLExpr(spe, ordannots = "segment", assay = 2, color = APOE)
png (file="RLE_speraw.png", width=1000, height=1000)
RLE_spe
dev.off()

plotRLExpr(spe, ordannots = "segment", assay = 2, color = segment)

set.seed(100)
spe <- scater::runPCA(spe)

drawPCA(spe, assay = 2, color = APOE)
RAWPCA <- drawPCA(spe, assay = 2, color = APOE)
png (file="plotScreePCA_speraw.png", width=1000, height=1000)
RAWPCA
dev.off()

pca_results <- reducedDim(spe, "PCA")
drawPCA(spe, precomputed = pca_results, col = APOE)
plotScreePCA(spe, precomputed = pca_results)
png (file="plotScreePCA_speraw.png", width=1000, height=1000)
plotScreePCA(spe, precomputed = pca_results)
dev.off()

plotPairPCA(spe, col = APOE, precomputed = pca_results, n_dimension = 5)
png (file="plotPairPCA_speraw.png", width=1000, height=1000)
plotPairPCA(spe, col = APOE, precomputed = pca_results, n_dimension = 5)
dev.off()

plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = APOE)
png (file="plotPCAbiplot_speraw.png", width=1000, height=1000)
plotPCAbiplot(spe, n_loadings = 10, precomputed = pca_results, col = APOE)
dev.off()

standR::plotMDS(spe, color = APOE)


set.seed(100)
spe <- scater::runUMAP(spe)
spe <- scater::runUMAP(spe, dimred = "PCA")
plotDR(spe, dimred = "UMAP", col = `slide name`)
png (file="UMAP_speraw.png", width=1000, height=1000)
plotDR(spe, dimred = "UMAP", col = `slide name`)
dev.off()

spe <- scater::runTSNE(spe)
plotDR(spe, dimred = "TSNE", col = `slide name`)
png (file="TSNE_speraw.png", width=1000, height=1000)
plotDR(spe, dimred = "TSNE", col = `slide name`)
dev.off()

##SPE_RUV_2. RUV NORMALISATION

spe <- findNCGs(spe, batch_name = "slide name", top_n = 300)

metadata(spe) |> names()


for(i in seq(10)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "APOE", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = APOE, title = paste0("DIM4, k = ", i)))
  
  png (file=paste0("DIM4k = ", i, ".png"), width=1000, height=1000)
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 4, color = APOE, title = paste0("DIM4, k = ", i)))
  dev.off()
  
}

for(i in seq(10)){
  spe_ruv <- geomxBatchCorrection(spe, factors = "APOE", 
                                  NCGs = metadata(spe)$NCGs, k = i)
  
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 2, color = APOE, title = paste0("DIM2, k = ", i)))
  
  png (file=paste0("DIM2k = ", i, ".png"), width=1000, height=1000)
  print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 2, color = APOE, title = paste0("DIM2, k = ", i)))
  dev.off()
  
}






#SELECT K-Value from RUV Normalisation##############

KVALDLG <-  dlgInput("ENTER K_VALUE (Values 1-10)")
if (KVALDLG$res == "1"){ KVAL <- 1}
if (KVALDLG$res == "2"){ KVAL <- 2}
if (KVALDLG$res == "3"){ KVAL <- 3}
if (KVALDLG$res == "4"){ KVAL <- 4}
if (KVALDLG$res == "5"){ KVAL <- 5}
if (KVALDLG$res == "6"){ KVAL <- 6}
if (KVALDLG$res == "7"){ KVAL <- 7}
if (KVALDLG$res == "8"){ KVAL <- 8}
if (KVALDLG$res == "9"){ KVAL <- 9}
if (KVALDLG$res == "10"){ KVAL <- 10}

#SET SUBDIRECTORY FOR K-VALUE#
KDIR <- paste0(FILTDIR,"/DGE_KVAL", KVAL)

dir.create(KDIR)
setwd(KDIR)




#DE_EXPRESSION - EDGER/VOOM#############################################
#DE_EXPRESSION_1. Set K Value for spe-ruv
spe_ruv <- geomxBatchCorrection(spe, factors = "LBD", 
                                NCGs = metadata(spe)$NCGs, k = KVAL)

#DE_EXPRESSION_2. spe-ruv normalised plots

print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 2, color = LBD, title = paste0("DIM2, k = ", KVAL)))
png (file="KVAL_PCA.png", width=1000, height=1000)
print(plotPairPCA(spe_ruv, assay = 2, n_dimension = 2, color = LBD, title = paste0("DIM2, k = ", KVAL)))
dev.off()
set.seed(100)

spe_ruv <- scater::runPCA(spe_ruv)

pca_results_ruv <- reducedDim(spe_ruv, "PCA")

plotScreePCA(spe_ruv, precomputed = pca_results_ruv)
png (file="plotScreePCA_spe_ruv.png", width=1000, height=1000)
plotScreePCA(spe_ruv, precomputed = pca_results_ruv)
dev.off()


plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = LBD, title = paste0("DIM4, k = ", KVAL), n_dimension = 4)
png (file="plotPairPCA_spe_ruv.png", width=1000, height=1000)
plotPairPCA(spe_ruv, precomputed = pca_results_ruv, color = LBD, title = paste0("DIM4, k = ", KVAL), n_dimension = 4)
dev.off()


spe_ruv <- scater::runUMAP(spe_ruv, dimred = "PCA")
plotDR(spe_ruv, dimred = "UMAP", col = LBD)
png (file="UMAP_spe_ruv_LBD.png", width=1000, height=1000)
plotDR(spe_ruv, dimred = "UMAP", col = LBD)
dev.off()

plotDR(spe_ruv, dimred = "UMAP", col = segment)
png (file="UMAP_spe_ruv_slide.png", width=1000, height=1000)
plotDR(spe_ruv, dimred = "UMAP", col = segment)
dev.off()


spe_ruv <- scater::runTSNE(spe_ruv, dimred = "PCA")
plotDR(spe_ruv, dimred = "TSNE", col = LBD)
png (file="tSNE_spe_ruv_LBD.png", width=1000, height=1000)
plotDR(spe_ruv, dimred = "TSNE", col = LBD)
dev.off()

plotDR(spe_ruv, dimred = "TSNE", col = segment)
png (file="tSNE_spe_ruv_LBD.png", width=1000, height=1000)
plotDR(spe_ruv, dimred = "TSNE", col = segment)
dev.off()

colData(spe_ruv)[,seq(ncol(colData(spe_ruv))-3, ncol(colData(spe_ruv)))] |>
  head()

#DE_EXPRESSION_3. Voom-StandR Differential Expression

library(edgeR)
library(limma)

dge <- SE2DGEList(spe_ruv)

if  (KVAL == 1) {design <- model.matrix(~ 0 + LBD + ruv_W1, data = colData(spe_ruv))}
if  (KVAL == 2) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2, data = colData(spe_ruv))}
if  (KVAL == 3) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W3, data = colData(spe_ruv))}
if  (KVAL == 4) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4, data = colData(spe_ruv))}
if  (KVAL == 5) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4 + ruv_W5, data = colData(spe_ruv))}
if  (KVAL == 6) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4 + ruv_W5 + ruv_W6, data = colData(spe_ruv))}
if  (KVAL == 7) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4 + ruv_W5 + ruv_W6 + ruv_W7, data = colData(spe_ruv))}
if  (KVAL == 8) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4 + ruv_W5 + ruv_W6 + ruv_W7 + ruv_W8, data = colData(spe_ruv))}
if  (KVAL == 9) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4 + ruv_W5 + ruv_W6 + ruv_W7 + ruv_W8 + ruv_W9, data = colData(spe_ruv))}
if  (KVAL == 10) {design <- model.matrix(~ 0 + LBD + ruv_W1+ ruv_W2+ ruv_W1+ ruv_W2+ ruv_W3 + ruv_W4 + ruv_W5 + ruv_W6 + ruv_W7 + ruv_W8 + ruv_W9 + ruv_W10, data = colData(spe_ruv))}

#DE_EXPRESSION_4. Define Design
design 

colnames(design)
colnames(design) <- gsub("^LBD","",colnames(design))
colnames(design) <- gsub(" ","_",colnames(design))
colnames(design)

#DE_EXPRESSION_5. Define Contrasts
contr.matrix <- makeContrasts( 
  CTRLvLBD = CTRL- LBD,
  levels = colnames(design))

#DE_EXPRESSION_6. Filter Genes by expression
keep <- filterByExpr(dge, design)
table(keep)
rownames(dge)[!keep]
dge_all <- dge[keep, ]

#DE_EXPRESSION_7.Estimate Dispersion
dge_all <- estimateDisp(dge_all, design = design, robust = TRUE)

#DE_EXPRESSION_8.Plot Biological coefficient of variation (BCV). 
plotBCV(dge_all)

bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)

highbcv <- bcv_df$BCV > 0.4
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4, size=2)

png (file="bcv.png", width=1000, height=1000)
plotBCV(dge_all)
bcv_df <- data.frame(
  'BCV' = sqrt(dge_all$tagwise.dispersion),
  'AveLogCPM' = dge_all$AveLogCPM,
  'gene_id' = rownames(dge_all)
)

highbcv <- bcv_df$BCV > 0.4
highbcv_df <- bcv_df[highbcv, ]
points(highbcv_df$AveLogCPM, highbcv_df$BCV, col = "red")
text(highbcv_df$AveLogCPM, highbcv_df$BCV, labels = highbcv_df$gene_id, pos = 4, size=2)
dev.off()

#DE_EXPRESSION_9.Plot Voom Mean Variance Trend.
v <- voom(dge_all, design, plot = TRUE)
png (file="voom.png", width=1000, height=1000)
v <- voom(dge_all, design, plot = TRUE) 
dev.off()


#DE_EXPRESSION_10.Perform linear Model fit and differential expression using StandR.
fit <- lmFit(v)

fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)

efit <- eBayes(fit_contrast, robust = TRUE)

top.table <- topTable(efit, sort.by = "P", n = Inf)
head(top.table, 20)

results_efit<- decideTests(efit, p.value = 0.05)

summary_efit <- summary(results_efit)
summary_efit

#DE_EXPRESSION_11.DE Results.
Results_Summary <- as.data.frame.matrix(summary_efit)
grid.newpage()
grid.table(Results_Summary)
png (file="Results_Summary.png", width=1000, height=1000)
grid.table(Results_Summary)
dev.off()


res_efit<-as.data.frame(results_efit)

#DE_EXPRESSION_12.DE Result export and Plots.
library(ggrepel)
library(tidyverse)



de_results <- topTable(efit, coef = 1, sort.by = "P", n = Inf)

de_genes_toptable <- topTable(efit, coef = 1, sort.by = "P", n = Inf, p.value = 0.05)

library(xlsx)
print_results <- de_results
row.names(print_results) <- NULL
write.xlsx(print_results, "DEResults.xlsx", sheetName="Sheet1")


png (file="MAplot.png", width=1000, height=1000)
de_results %>% 
  mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                     ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
  ggplot(aes(AveExpr, logFC, col = DE)) + 
  geom_point(shape = 1, size = 1) + 
  geom_text_repel(data = de_results %>% 
                    mutate(DE = ifelse(logFC > 0 & adj.P.Val <0.05, "UP", 
                                       ifelse(logFC <0 & adj.P.Val<0.05, "DOWN", "NOT DE"))) %>%
                    rownames_to_column(), aes(label = rowname)) +
  theme_bw() +
  xlab("Average log-expression") +
  ylab("Log-fold-change") +
  ggtitle("CTRL vs LBD (limma-voom)") +
  scale_color_manual(values = c("blue","gray","red")) +
  theme(text = element_text(size=15))
dev.off()


library(EnhancedVolcano)

SIG_VOLCANO <- EnhancedVolcano(de_results,
                               lab = rownames(de_results),
                               x = 'logFC',
                               y = 'adj.P.Val',
                               pCutoff = 0.05,
                               FCcutoff = 0.5,
                               pointSize = 0.5,
                               xlim = c(-3, 3),
                               cutoffLineType = 'twodash',
                               cutoffLineWidth = 0.5,
                               hline = c(0.01,
                                         0.01 * 0.001,
                                         0.001 * 0.0001),
                               
                               hlineCol = c('black', 'blue', 'blue'),
                               hlineType = c('twodash'),
                               hlineWidth = c(0.5),
                               labSize = 5.0,
                               title = "EdgeR results - Differential expression",
                               subtitle = "Log2FC cutoff, 0.5; FDR, 0.05",
                               legendLabels=c('Not sig.','Log2FC','p-value',
                                              'p-value & Log2FC'),
                               legendPosition = 'bottom',
                               legendLabSize = 10,
                               legendIconSize = 5.0,
                               drawConnectors = FALSE,
                               widthConnectors = 1,
                               max.overlaps = 200,
                               colAlpha = 4/5,
                               colConnectors = 'black',
                               selectLab = c("GFAP","GBA","LRRK2", "S100B", "SNCA", "PRKN", "PARK2", "PARK7", "APP", "MAPT", "PINK1", "NFKB", "AGT", "SOX9", "AQP4")) 
SIG_VOLCANO
png (file="SigVolcanoplot.png", width=1000, height=1000)
SIG_VOLCANO
dev.off()




FULL_VOLCANO <- EnhancedVolcano(de_results,
                                lab = rownames(de_results),
                                x = 'logFC',
                                y = 'adj.P.Val',
                                pCutoff = 0.05,
                                FCcutoff = 0.5,
                                pointSize = 0.5,
                                xlim = c(-3, 3),
                                cutoffLineType = 'twodash',
                                cutoffLineWidth = 0.5,
                                hline = c(0.01,
                                          0.01 * 0.001,
                                          0.001 * 0.0001),
                                
                                hlineCol = c('black', 'blue', 'blue'),
                                hlineType = c('twodash'),
                                hlineWidth = c(0.5),
                                labSize = 4,
                                title = "EdgeR results - Differential expression",
                                subtitle = "Log2FC cutoff, 0.5; FDR, 0.05",
                                legendLabels=c('Not sig.','Log2FC','p-value',
                                               'p-value & Log2FC'),
                                legendPosition = 'bottom',
                                legendLabSize = 10,
                                legendIconSize = 5.0,
                                drawConnectors = FALSE,
                                widthConnectors = .1,
                                max.overlaps = 30,
                                colAlpha = 4/5,
                                colConnectors = 'black',
                                selectLab = rownames(de_genes_toptable))
FULL_VOLCANO
png (file="FullVolcanoplot.png", width=1000, height=1000)
FULL_VOLCANO
dev.off()


library(DT)

updn_cols <- c(RColorBrewer::brewer.pal(6, 'Greens')[2], RColorBrewer::brewer.pal(6, 'Purples')[2])

de_results %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  DT::datatable(caption = paste0('CTRLvsLBD K=',KVAL)) %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)


DT_result <- de_results %>% 
  dplyr::select(c("logFC", "AveExpr", "P.Value", "adj.P.Val")) %>%
  DT::datatable(caption = paste0('CTRLvsLBD K=',KVAL )) %>%
  DT::formatStyle('logFC',
                  valueColumns = 'logFC',
                  backgroundColor = DT::styleInterval(0, rev(updn_cols))) %>%
  DT::formatSignif(1:4, digits = 4)


DT::saveWidget(DT_result, 'DT_RESULT.html')

#GO-TERM_ANALYSIS#############################################
library(msigdb)
library(GSEABase)

#msigdb_hs <- getMsigdb(version = '7.2')
#msigdb_hs <- appendKEGG(msigdb_hs)

sc <- listSubCollections(msigdb_hs)

gsc <- c(subsetCollection(msigdb_hs, c('h')),
         subsetCollection(msigdb_hs, 'c2', sc[grepl("^CP:",sc)]),
         subsetCollection(msigdb_hs, 'c5', sc[grepl("^GO:",sc)])) %>%
  GeneSetCollection()

fry_indices <- ids2indices(lapply(gsc, geneIds), rownames(v), remove.empty = FALSE)
names(fry_indices) <- sapply(gsc, setName)

gsc_category <- sapply(gsc, function(x) bcCategory(collectionType(x)))
gsc_category <- gsc_category[sapply(fry_indices, length) > 5]

gsc_subcategory <- sapply(gsc, function(x) bcSubCategory(collectionType(x)))
gsc_subcategory <- gsc_subcategory[sapply(fry_indices, length) > 5]

fry_indices <- fry_indices[sapply(fry_indices, length) > 5]

names(gsc_category) = names(gsc_subcategory) = names(fry_indices)

fry_indices_cat <- split(fry_indices, gsc_category[names(fry_indices)])
fry_res_out <- lapply(fry_indices_cat, function (x) {
  limma::fry(v, index = x, design = design, contrast = contr.matrix[,1], robust = TRUE)
})

post_fry_format <- function(fry_output, gsc_category, gsc_subcategory){
  names(fry_output) <- NULL
  fry_output <- do.call(rbind, fry_output)
  fry_output$GenesetName <- rownames(fry_output)
  fry_output$GenesetCat <- gsc_category[rownames(fry_output)]
  fry_output$GenesetSubCat <- gsc_subcategory[rownames(fry_output)]
  return(fry_output)
}

fry_res_sig <- post_fry_format(fry_res_out, gsc_category, gsc_subcategory) %>%
  as.data.frame() %>%
  filter(FDR < 0.05) 

fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Up",) %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "red") +
  theme_bw() +
  coord_flip() +
  ggtitle("Up-regulated")


fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Down") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw() +
  coord_flip() +
  ggtitle("Down-regulated")


png (file="GOTERM_UP.png", width=1000, height=1000)
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Up") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "red") +
  theme_bw() +
  coord_flip() +
  ggtitle("Up-regulated")
dev.off()

png (file="GOTERM_DOWN.png", width=1000, height=1000)
fry_res_sig %>%
  arrange(FDR) %>%
  filter(Direction == "Down") %>%
  .[seq(20),] %>%
  mutate(GenesetName = factor(GenesetName, levels = .$GenesetName)) %>%
  ggplot(aes(GenesetName, -log(FDR))) +
  geom_bar(stat = "identity", fill = "blue") +
  theme_bw() +
  coord_flip() +
  ggtitle("Down-regulated")
dev.off()
