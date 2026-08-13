GSE275938_count <- read.csv(file = 'GSE275938_compiled_counts.csv', row.names = 1)
GSE275938_count <- t(GSE275938_count)
GSE275938_metadata <- read.csv(file = '/home/fanjie/Yao/GSE275938/GSE275938_cell_metadata.csv',row.names = 1)


GSE275938_sce <- CreateSeuratObject(counts = GSE275938_count,
                                    meta.data = GSE275938_metadata)

GSE275938_sce.list <- SplitObject(GSE275938_sce, split.by="dataset")

GSE275938_sce.list <- lapply(GSE275938_sce.list, function(x) {x <- SCTransform(x,verbose=FALSE)})

GSE275938_sce.list.features <- SelectIntegrationFeatures(object.list=GSE275938_sce.list, nfeatures=3000)
GSE275938_sce.list <- PrepSCTIntegration(object.list=GSE275938_sce.list, anchor.features=GSE275938_sce.list.features, verbose=FALSE)
GSE275938_sce.int.anchors <- FindIntegrationAnchors(GSE275938_sce.list, normalization.method="SCT", anchor.features=GSE275938_sce.list.features, verbose=FALSE)
GSE275938_sce.combined <- IntegrateData(anchorset=GSE275938_sce.int.anchors, normalization.method="SCT", verbose=FALSE)

DefaultAssay(GSE275938_sce.combined) <- "integrated"

GSE275938_sce.combined <- subset(GSE275938_sce.combined, dataset != 'Acute preterm injury 1')
GSE275938_sce.combined <- RunPCA(GSE275938_sce.combined, features = VariableFeatures(object = GSE275938_sce.combined))
GSE275938_sce.combined <- FindNeighbors(GSE275938_sce.combined, dims = 1:50)
GSE275938_sce.combined <- FindClusters(GSE275938_sce.combined, resolution = 2)
GSE275938_sce.combined <- RunUMAP(GSE275938_sce.combined, dims=1:50)



Idents(GSE275938_sce.combined) <- 'celltype'
DimPlot(GSE275938_sce.combined, label = F, cols = c(brewer.pal(9,"Set1"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"),brewer.pal(12,"Set3"),brewer.pal(8,"Dark2"),
                                                  brewer.pal(8,"Pastel2"),brewer.pal(8,"Set2"),brewer.pal(8,"Paired")))

### Figure 1

library(scplotter)
library(scCustomize)
library(viridis)
library(dplyr)
library(ggtext)

# Gene sets
# Mitophagy
mitophagy <- list(mitophagy = c("ATG12", "ATG5", "CSNK2A1", "CSNK2A2", "CSNK2B", "FUNDC1", "MAP1LC3A", "MAP1LC3B", "MFN1", "MFN2", "MTERF3", "PGAM5", "PINK1", "PRKN", "RPS27A", "SQSTM1", "SRC", "TOMM20", "TOMM22", "TOMM40", "TOMM5", "TOMM6", "TOMM7", "TOMM70", "UBA52", "UBB", "UBC", "ULK1", "VDAC1", "BNIP3", "BNIP3L", "OPTN", "NDP52", "ATG7", "BECN1", "LAMP1", "LAMP2", "DNM1L", "OPA1", "PPARGC1A", "TFAM", "NRF1"))

# cGAS-STING
cGAS_STING <- list(cGAS_STING = c('CGAS', 'TMEM173', 'TBK1', 'IRF3', 'IKBKG', 'CHUK', 'IKBKB', 'NFKB1', 'RELA', 'TREX1', 'ENPP1', 'NLRP4', 'NLRC3', 'TRIM21', 'TRIM56', 'IFI16', 'DDX41', 'ZDHHC1', 'PRKDC', 'XRCC6', 'STAT6', 'DTX4', 'XRCC5', 'MRE11'))

#cGAS_STING <- list(cGAS_STING = c("CGAS", "TMEM173", "TBK1", "IRF3", "IKBKG", "CHUK", "IKBKB", "NFKB1", "RELA", "TREX1", "ENPP1", "NLRP4", "NLRC3", "TRIM21", "TRIM56", "IFI16", "DDX41", "ZDHHC1", "MB21D1", "TMEM173", "ZBP1", "IKBKE", "IRF7", "NFKBIA", "STAT1", "STAT2", "IRF9", "JAK1", "TYK2", "ISG15", "IFI44", "IFI44L", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM2", "IFITM3", "HERC5", "USP18", "CXCL10", "CXCL11", "CCL5", "TNF", "IL6"))

# TLR9
TLR9 <- list(TLR9 = c(c("TLR9", "MYD88", "IRAK4", "IRAK1", "TRAF6", "MAP3K7", "TAB1", "TAB2", "TAB3", "CHUK", "IKBKB", "IKBKG", "NFKBIA", "NFKB1", "RELA", "MAP2K", "MAPK", "FOS", "JUN", "IRF7", "TNF", "IL6", "IL12", "IFNA", "IFNB", "IRF5", "IFNB1", "IFNA1", "IL1B", "CCL5", "CXCL10")))

# inflammasome
#inflammasome <- list(inflammasome = c("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD", "NEK7", "CARD8", "IL1RN", "NLRP1", "NLRC4", "AIM2", "MEFV", "CASP4", "CASP5", "IL33"))

inflammasome <- list(inflammasome = c("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD", "NEK7","NLRP1", "NLRC4", "AIM2", "NLRP6", "NLRP12","CARD8", "IL1RN", "SIRT1", "SIRT2", "TXNIP","P2RX7", "PANX1", "TLR4", "TLR2","IL1A", "IL33", "CASP4", "CASP5", "GSDME","NFKB1", "RELA", "NFKBIA"))

