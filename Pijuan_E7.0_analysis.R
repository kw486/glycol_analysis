suppressMessages(c(library(Seurat),
                   library(ggplot2),
                   library(dplyr),
                   library(tidyverse),
                   library(rdist),
                   library(tidyr),
                   library(pheatmap),
                   library(pals), 
                   library(readxl),
                   library(EnhancedVolcano),
                   library(writexl)
                   
))

#read in glycolysis genes from uniprot
load("glycolgenes.RData")
glucose.transporters <- c('Slc2a8', 'Slc2a9', 'Slc5a1','Slc2a2','Slc2a3','Slc2a10', 'Slc2a12','Slc2a1','Slc2a4','Slc2a5') 
lactate.transporters <- c('Slc5a12','Slc5a8','Slc16a1','Slc16a3','Slc16a4','Slc16a5')
unip.glycol <- c(unip.glycol,glucose.transporters,lactate.transporters)

#read in oxphos genes from uniprot
load("oxphosgenes.RData")

#load mouse data
load("pijuan_data.RData")

#subset to E7.0 and load celltype annotations
timestage <- "E7.0"
pijuan.subset<- subset(pijuan,subset= stage %in% timestage)
Idents(pijuan.subset) <- "celltype.pijuan"

#normalise
pijuan.subset <- NormalizeData(pijuan.subset, normalization.method = "LogNormalize")

#differential expression testing
cluster.markers <- FindMarkers(pijuan.subset, ident.1="Primitive Streak", ident.2=c("ExE ectoderm", "Surface ectoderm","Rostral neurectoderm"), min.pct = 0.01, logfc.threshold = 0.05) 
cluster.markers$gene <- rownames(cluster.markers)
cluster.markers <- cluster.markers %>% filter (!(p_val_adj == 1)) 
cluster.markers <- cluster.markers %>% filter (!(p_val_adj == 0)) 
sig.glycol <- cluster.markers %>% filter(gene %in% unip.glycol & p_val_adj <=0.05)

#volcano plot
cluster.markers$anno<-ifelse(cluster.markers$gene %in% c(unip.glycol,oxphosmarkers), 1, 0 )
cluster.markers <- cluster.markers %>%
  arrange(anno)

annogene <- cluster.markers$gene[cluster.markers$gene %in% c(unip.glycol,oxphosmarkers) & abs(cluster.markers$avg_log2FC)>=0.5 & cluster.markers$p_val_adj<0.05 ]

keyvals.col <- ifelse(cluster.markers$gene %in% unip.glycol, "deeppink",
                      ifelse(cluster.markers$gene %in% oxphosmarkers, "#4381C1", "grey"))

keyvals.col[is.na(names(keyvals.col))] <- "grey"
names(keyvals.col)[keyvals.col == "grey"] <- 'All Genes'
names(keyvals.col)[keyvals.col == "deeppink"] <- 'Glycolysis'
names(keyvals.col)[keyvals.col == "#4381C1"] <- 'Oxidative Phosphorylation'


p <- EnhancedVolcano(cluster.markers , 
                     lab = rownames(cluster.markers),
                     selectLab = annogene , 
                     colCustom = keyvals.col,
                     x = 'avg_log2FC',
                     y = 'p_val_adj',
                     title = paste(unique(pijuan.subset$stage), "Primitive Streak vs Ectoderm", sep=" | "),
                     caption="",
                     subtitle = "",
                     legendLabels = c('Glycolysis', 'Oxidative Phosphorylation'),
                     colAlpha = ifelse(cluster.markers$gene %in% c(oxphosmarkers,unip.glycol), 1,0.1),
                     pointSize = 4,
                     drawConnectors = T,
                     colConnectors = "grey",
                     widthConnectors = 0,
                     arrowheads = F,
                     labSize = 5,
                     pCutoff = 0.05,
                     FCcutoff = 0.5, 
) + xlim(-7,7) + ylim(0,300)

# export to pdf
pdf("Volcano_PrimEcto_E7.0_jitter.pdf", width = 16, height = 12)
print(p)
dev.off()

#pseudobulk expression
mat<- AggregateExpression(
  pijuan.subset,
  group.by = 'celltype.pijuan',
  assays = 'RNA',
  features = sig.glycol$gene) 
mat<-as.data.frame(mat[[1]]) 

#normalise
mat <- NormalizeData(mat)

#select cell types >300 cells
celltype.numbers <- as.data.frame(table(pijuan.subset$celltype.pijuan))
celltype.numbers %>% filter(Freq>300)
mat <- as.data.frame(mat) %>% dplyr::select(celltype.nos)

#heatmap
p <- pheatmap(
  mat,
  color=viridis::viridis(100), #scale to make between -1 0 1
  cluster_rows = T,
  cluster_cols = T, 
  treeheight_row = 0,
  scale = "row") 

#export
pdf("Heatmap_GlycolMarkers.pdf", width = 6, height = 8)
print(p)
dev.off()