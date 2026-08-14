GSE275938_count <- read.csv(file = 'GSE275938_compiled_counts.csv', row.names = 1)
GSE275938_count <- t(GSE275938_count)
GSE275938_metadata <- read.csv(file = 'GSE275938_cell_metadata.csv', row.names = 1)


GSE275938_sce <- CreateSeuratObject(counts = GSE275938_count,
                                    meta.data = GSE275938_metadata)

GSE275938_sce.list <- SplitObject(GSE275938_sce, split.by = "dataset")

GSE275938_sce.list <- lapply(GSE275938_sce.list, function(x) {x <- SCTransform(x, verbose = FALSE)})

GSE275938_sce.list.features <- SelectIntegrationFeatures(object.list = GSE275938_sce.list, nfeatures = 3000)
GSE275938_sce.list <- PrepSCTIntegration(object.list = GSE275938_sce.list, anchor.features = GSE275938_sce.list.features, verbose = FALSE)
GSE275938_sce.int.anchors <- FindIntegrationAnchors(GSE275938_sce.list, normalization.method = "SCT", anchor.features = GSE275938_sce.list.features, verbose = FALSE)
GSE275938_sce.combined <- IntegrateData(anchorset = GSE275938_sce.int.anchors, normalization.method = "SCT", verbose = FALSE)

DefaultAssay(GSE275938_sce.combined) <- "integrated"

GSE275938_sce.combined <- subset(GSE275938_sce.combined, dataset != 'Acute preterm injury 1')
GSE275938_sce.combined <- RunPCA(GSE275938_sce.combined, features = VariableFeatures(object = GSE275938_sce.combined))
GSE275938_sce.combined <- FindNeighbors(GSE275938_sce.combined, dims = 1:50)
GSE275938_sce.combined <- FindClusters(GSE275938_sce.combined, resolution = 2)
GSE275938_sce.combined <- RunUMAP(GSE275938_sce.combined, dims = 1:50)



library(scCustomize)
library(viridis)
library(ggpubr)

# Gene sets
# Mitophagy
mitophagy <- list(mitophagy = c("ATG12", "ATG5", "CSNK2A1", "CSNK2A2", "CSNK2B", "FUNDC1", "MAP1LC3A", "MAP1LC3B", "MFN1", "MFN2", "MTERF3", "PGAM5", "PINK1", "PRKN", "RPS27A", "SQSTM1", "SRC", "TOMM20", "TOMM22", "TOMM40", "TOMM5", "TOMM6", "TOMM7", "TOMM70", "UBA52", "UBB", "UBC", "ULK1", "VDAC1", "BNIP3", "BNIP3L", "OPTN", "NDP52", "ATG7", "BECN1", "LAMP1", "LAMP2", "DNM1L", "OPA1", "PPARGC1A", "TFAM", "NRF1"))

# cGAS-STING
cGAS_STING <- list(cGAS_STING = c('CGAS', 'TMEM173', 'TBK1', 'IRF3', 'IKBKG', 'CHUK', 'IKBKB', 'NFKB1', 'RELA', 'TREX1', 'ENPP1', 'NLRP4', 'NLRC3', 'TRIM21', 'TRIM56', 'IFI16', 'DDX41', 'ZDHHC1', 'PRKDC', 'XRCC6', 'STAT6', 'DTX4', 'XRCC5', 'MRE11'))

# TLR9
TLR9 <- list(TLR9 = c("TLR9", "MYD88", "IRAK4", "IRAK1", "TRAF6", "MAP3K7", "TAB1", "TAB2", "TAB3", "CHUK", "IKBKB", "IKBKG", "NFKBIA", "NFKB1", "RELA", "MAP2K", "MAPK", "FOS", "JUN", "IRF7", "TNF", "IL6", "IL12", "IFNA", "IFNB", "IRF5", "IFNB1", "IFNA1", "IL1B", "CCL5", "CXCL10"))

# inflammasome
inflammasome <- list(inflammasome = c("NLRP3", "PYCARD", "CASP1", "IL1B", "IL18", "GSDMD", "NEK7","NLRP1", "NLRC4", "AIM2", "NLRP6", "NLRP12","CARD8", "IL1RN", "SIRT1", "SIRT2", "TXNIP","P2RX7", "PANX1", "TLR4", "TLR2","IL1A", "IL33", "CASP4", "CASP5", "GSDME","NFKB1", "RELA", "NFKBIA"))

GSE275938_sce.combined <- AddModuleScore(object = GSE275938_sce.combined, 
                                         features = mitophagy,
                                         name = 'mitophagy',
                                         assay = 'SCT',
                                         slot = 'data')
GSE275938_sce.combined <- AddModuleScore(object = GSE275938_sce.combined, 
                                         features = cGAS_STING,
                                         name = 'cGAS_STING',
                                         assay = 'SCT',
                                         slot = 'data')
GSE275938_sce.combined <- AddModuleScore(object = GSE275938_sce.combined, 
                                         features = TLR9,
                                         name = 'TLR9',
                                         assay = 'SCT',
                                         slot = 'data')
GSE275938_sce.combined <- AddModuleScore(object = GSE275938_sce.combined, 
                                         features = inflammasome,
                                         name = 'inflammasome',
                                         assay = 'SCT',
                                         slot = 'data')

### Figure 1 ###
GSE275938_sce.combined_celltype <- GSE275938_sce.combined@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% cbind(Celltype = GSE275938_sce.combined@meta.data$celltype)

GSE275938_sce.combined_celltype$Celltype <- factor(GSE275938_sce.combined_celltype$Celltype, 
                                levels = c("Basal", "Multiciliated", "Secretory MUC5B", "Secretory -3A1, -3A2", "RASC", "AT1", "AT2", # Epithelium
                                           "gCap", "aCap", "abCap", "Arterial EC", "Pulmonary venous EC", "Systemic venous EC", "Lymphatic", # Endothelium
                                           "Alveolar FB", "Adventitial FB", "Pericyte", "Activated FB", "Ductal MyoFB", "VSMC", "Alveolar MyoFB", # Stroma
                                           "Alveolar Macrophage", "Monocyte", "B Cell", "Plasma cell", "T Cell", "NKT Cell", "NK Cell", "cDC", "pDC", "Mast cell", "Basophil", "Neutrophil"
                                ))
GSE275938_celltype_col <- c(brewer.pal(12, "Set3"),brewer.pal(9, "Pastel1"),brewer.pal(8, "Accent"),brewer.pal(9, "Set1"),brewer.pal(8, "Dark2"),
                      brewer.pal(8, "Pastel2"),brewer.pal(8, "Set2"),brewer.pal(8, "Paired"))

GSE275938_cluster_centers <- GSE275938_sce.combined_celltype%>%
  group_by(Celltype)%>%
  dplyr::summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))%>%
  dplyr::mutate(cluster_num = paste0("c", "", seq_along(Celltype))) %>%
  as.data.frame()

GSE275938_cluster_centers$x <- 1
GSE275938_cluster_centers$y <- seq_len(nrow(GSE275938_cluster_centers))

cleg_color <- GSE275938_celltype_col
cleg <- ggplot(GSE275938_cluster_centers,aes(x = x, y = y))+
  geom_point(
    aes(color = Celltype),
    show.legend = FALSE,
    size = 5
  )+
  geom_text(aes(label = cluster_num))+
  geom_text(
    aes(label=Celltype),
    hjust = 0,
    nudge_x = 0.1,
    size = 4
  )+
  scale_color_manual(values = cleg_color)+
  scale_y_reverse()+
  xlim(0, 2)+
  theme_void()
cleg


ggplot(GSE275938_sce.combined_celltype,aes(x = UMAP_1, y = UMAP_2, color = Celltype))+
  geom_point(size = 0.3, alpha = 1)+
  geom_point(data = GSE275938_cluster_centers,aes(x = UMAP_1, y = UMAP_2),
             size = 8, color = "white", fill = 'black', alpha = 0.5, stroke = 2)+
  geom_text(data = GSE275938_cluster_centers,
            aes(x = UMAP_1, y = UMAP_2, label = cluster_num),
            size = 3.5,
            color = "black")+
  scale_color_manual(values = GSE275938_celltype_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background = element_rect(fill = "white"),
        plot.title = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none',
        legend.text = element_text(family = 'Arial', size=18, color = "black"),
        legend.key = element_rect(fill = 'transparent'),
        legend.background = element_rect(fill = "transparent"),
        legend.justification = c(1, 0),
        legend.key.size=unit(0.5, "cm"))+
  geom_segment(aes(x = min(GSE275938_sce.combined_celltype$UMAP_1) , y = min(GSE275938_sce.combined_celltype$UMAP_2) ,
                   xend = min(GSE275938_sce.combined_celltype$UMAP_1) +3, yend = min(GSE275938_sce.combined_celltype$UMAP_2) ),
               colour = "black", size = 0.5, arrow = arrow(length = unit(0.3, "cm")))+ 
  geom_segment(aes(x = min(GSE275938_sce.combined_celltype$UMAP_1)  , y = min(GSE275938_sce.combined_celltype$UMAP_2)  ,
                   xend = min(GSE275938_sce.combined_celltype$UMAP_1) , yend = min(GSE275938_sce.combined_celltype$UMAP_2) + 3),
               colour = "black", size = 0.5, arrow = arrow(length = unit(0.3, "cm"))) +
  annotate("text", x = min(GSE275938_sce.combined_celltype$UMAP_1) + 1.5, y = min(GSE275938_sce.combined_celltype$UMAP_2)-1, label = "UMAP_1",
           color = "black", size = 5, fontface = "plain", family = 'Arial') + 
  annotate("text", x = min(GSE275938_sce.combined_celltype$UMAP_1) -1, y = min(GSE275938_sce.combined_celltype$UMAP_2)+1.5, label = "UMAP_2",
           color = "black", size = 5, fontface = "plain" , family = 'Arial', angle = 90) 


pal <- viridis(n = 33, option = "A", direction = -1)
FeaturePlot_scCustom(seurat_object = GSE275938_sce.combined,
                     colors_use = pal,
                     features = 'Mitophagy pathway score')+
  theme(title = element_text(size = 20),
    panel.border=element_rect(fill = NA, color = "black",
                                  size = 1, linetype = "solid"))

FeaturePlot_scCustom(seurat_object = GSE275938_sce.combined,
                     colors_use = pal,
                     features = 'cGAS-STING pathway score')+
  theme(title = element_text(size = 20),
    panel.border=element_rect(fill = NA, color = "black",
                                  size = 1, linetype = "solid"))

FeaturePlot_scCustom(seurat_object = GSE275938_sce.combined,
                     colors_use = pal,
                     features = 'TLR9 pathway score')+
  theme(title = element_text(size = 20),
    panel.border=element_rect(fill = NA, color = "black",
                                  size = 1, linetype = "solid"))

FeaturePlot_scCustom(seurat_object = GSE275938_sce.combined,
                     colors_use = pal,
                     features = 'Inflammasome pathway score')+
  theme(title = element_text(size = 20),
    panel.border=element_rect(fill = NA, color = "black",
                                  size = 1, linetype = "solid"))

### Figure 2 ###
comparisons <- list( c("Term infant", "BPD"))

## 1. Miotphagy ##
# Global
p <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD')), 
             assay = 'SCT', feature = 'Mitophagy', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Global')+
  ylab('Estimated score')

p + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.25, max(p$data[tail('Mitophagy', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.42),
                     method = 'wilcox.test',
                     label = "p.signif")


# Epithelial cells
p1 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Epithelial'), assay = 'SCT', feature = 'Mitophagy', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Epithelial cells')+
  ylab('Estimated score')

p1 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p1$data[tail('Mitophagy', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.38),
                     method = 'wilcox.test',
                     label = "p.signif")


# Endothelial cells
p2 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Endothelial'), assay = 'SCT', feature = 'Mitophagy', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Endothelial cells')+
  ylab('Estimated score')

p2 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.24, max(p2$data[tail('Mitophagy', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.41),
                     method = 'wilcox.test',
                     label = "p.signif")


# Immune
p3 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Immune'), assay = 'SCT', feature = 'Mitophagy', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Immune cells')+
  ylab('Estimated score')

p3 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face ='bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5,face = 'bold',colour = 'black'),
        #panel.border = element_rect(color = 'black',size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.26, max(p3$data[tail('Mitophagy', 1)]+0.09))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.43),
                     method = 'wilcox.test',
                     label = "p.signif")


# Mesenchyme
p4 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Mesenchymal'), assay = 'SCT', feature = 'Mitophagy', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Mesenchymal cells')+
  ylab('Estimated score')

p4 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p4$data[tail('Mitophagy', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.41),
                     method = 'wilcox.test',
                     label = "p.signif")


## 2. cGAS-STING ##
# Global
p5 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD')), 
             assay = 'SCT', feature = 'cGAS-STING', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Global')+
  ylab('Estimated score')

p5 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"),0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p5$data[tail('cGAS-STING', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.38),
                     method = 'wilcox.test',
                     label = "p.signif")


# Epithelium
p6 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Epithelial'), assay = 'SCT', feature = 'cGAS-STING', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Epithelial cells')+
  ylab('Estimated score')

p6 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p6$data[tail('cGAS-STING', 1)]+0.06))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.19),
                     method = 'wilcox.test',
                     label = "p.signif")


# Endothelium
p7 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Endothelial'), assay = 'SCT', feature = 'cGAS-STING', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Endothelial cells')+
  ylab('Estimated score')

p7 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color='black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p7$data[tail('cGAS-STING', 1)]+0.06))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.32),
                     method = 'wilcox.test',
                     label = "p.signif")


# Immune
p8 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Immune'), assay = 'SCT', feature = 'cGAS-STING', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Immune cells')+
  ylab('Estimated score')

p8 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p8$data[tail('cGAS-STING', 1)]+0.09))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.38),
                     method = 'wilcox.test',
                     label = "p.signif")


# Mesenchyme
p9 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Mesenchymal'), assay = 'SCT', feature = 'cGAS-STING', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black",fill = "white")+NoLegend()+
  ggtitle('Mesenchymal cells')+
  ylab('Estimated score')

p9 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p9$data[tail('cGAS-STING', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.32),
                     method = 'wilcox.test',
                     label = "p.signif")


## 3. TLR9 ##
# Global
p10 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD')), 
             assay = 'SCT', feature = 'TLR9', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Global')+
  ylab('Estimated score')

p10 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.26, max(p10$data[tail('TLR9', 1)]+0.10))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.63),
                     method = 'wilcox.test',
                     label = "p.signif")


# Epithelium
p11 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Epithelial'), assay = 'SCT', feature = 'TLR9', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Epithelial cells')+
  ylab('Estimated score')

p11 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.20, max(p11$data[tail('TLR9', 1)]+0.07))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.45),
                     method = 'wilcox.test',
                     label = "p.signif")


# Endothelium
p12 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Endothelial'), assay = 'SCT', feature = 'TLR9', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Endothelial cells')+
  ylab('Estimated score')

p12 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1,fill = NA),
        panel.grid = element_blank())+
  ylim(-0.24, max(p12$data[tail('TLR9', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.47),
                     method = 'wilcox.test',
                     label = "p.signif")


# Immune
p13 <- VlnPlot(subset(GSE275938_sce.combined1, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Immune'), assay = 'SCT', feature = 'TLR9', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Immune cells')+
  ylab('Estimated score')

p13 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p13$data[tail('TLR9', 1)]+0.1))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.63),
                     method = 'wilcox.test',
                     label = "p.signif")


# Mesenchyme
p14 <- VlnPlot(subset(GSE275938_sce.combined1, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Mesenchymal'), assay = 'SCT', feature = 'TLR9', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Mesenchymal cells')+
  ylab('Estimated score')

p14 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.24, max(p14$data[tail('TLR9', 1)]+0.09))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.54),
                     method = 'wilcox.test',
                     label = "p.signif")


## 4. Inflammasome ##
# Global
p15 <- VlnPlot(subset(GSE275938_sce.combined1, group %in% c('Term infant', 'BPD')), 
             assay = 'SCT', feature = 'Inflammasome', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Global')+
  ylab('Estimated score')

p15 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1,fill = NA),
        panel.grid = element_blank())+
  ylim(-0.26, max(p15$data[tail('Inflammasome', 1)]+0.10))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.75),
                     method = 'wilcox.test',
                     label = "p.signif")


# Epithelium
p16 <- VlnPlot(subset(GSE275938_sce.combined1, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Epithelial'), assay = 'SCT', feature = 'Inflammasome', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black",fill = "white")+NoLegend()+
  ggtitle('Epithelial cells')+
  ylab('Estimated score')

p16 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.22, max(p16$data[tail('Inflammasome', 1)]+0.04))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.22),
                     method = 'wilcox.test',
                     label = "p.signif")


# Endothelium
p17 <- VlnPlot(subset(GSE275938_sce.combined1, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Endothelial'), assay = 'SCT', feature = 'Inflammasome', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Endothelial cells')+
  ylab('Estimated score')

p17 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.20, max(p17$data[tail('Inflammasome', 1)]+0.08))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.34),
                     method = 'wilcox.test',
                     label = "p.signif")


# Immune
p18 <- VlnPlot(subset(GSE275938_sce.combined1, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Immune'), assay = 'SCT', feature = 'Inflammasome', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Immune cells')+
  ylab('Estimated score')

p18 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.23, max(p18$data[tail('Inflammasome', 1)]+0.1))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.76),
                     method = 'wilcox.test',
                     label = "p.signif")


# Mesenchyme
p19 <- VlnPlot(subset(GSE275938_sce.combined, group %in% c('Term infant', 'BPD') &
                      celltype_lineage == 'Mesenchymal'), assay = 'SCT', feature = 'Inflammasome', group.by = 'group', pt.size = 0)+
  geom_boxplot(width = .15, col = "black", fill = "white")+NoLegend()+
  ggtitle('Mesenchymal cells')+
  ylab('Estimated score')

p19 + 
  scale_fill_manual(values = alpha(c("#f4a4c9", "#4955d0"), 0.5))+
  theme_classic()+
  theme(axis.text.y = element_text(size = 18, face = 'bold', colour = 'black'),
        axis.text.x = element_text(colour = 'black', 
                                   face = 'bold',
                                   size = 18, 
                                   angle = 45, 
                                   vjust = 1, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,
                                    colour = 'black',
                                    face = 'bold',
                                    angle = 90, 
                                    vjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 20, hjust = 0.5, face = 'bold', colour = 'black'),
        #panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.grid = element_blank())+
  ylim(-0.24, max(p$data[tail('Inflammasome', 1)]+0.07))+
  stat_compare_means(comparisons = comparisons,
                     tip.length = 0,
                     label.y = c(0.24),
                     method = 'wilcox.test',
                     label = "p.signif")
