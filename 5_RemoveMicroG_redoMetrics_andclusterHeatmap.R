######remove microglia cliuster 8 to get better cluster resolution
Justine_sub_no8<-subset(x = Justine_sub, idents = 8, invert=T)
Justine_sub_no8<- NormalizeData(Justine_sub_no8) %>% ScaleData()
dim(Justine_sub_no8)
#[1] 17236  7163
Justine_sub_no8 <- FindVariableFeatures(Justine_sub_no8, nfeatures = 8618)#50%
top10 <- head(VariableFeatures(Justine_sub_no8), 10)
plot1 <- VariableFeaturePlot(Justine_sub_no8, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

set.seed(42)
Justine_sub_no8 <- RunPCA(Justine_sub_no8, features = VariableFeatures(object = Justine_sub_no8))
Justine_sub_no8 <- FindNeighbors(Justine_sub_no8, dims = 1:10)
Justine_sub_no8 <- FindClusters(Justine_sub_no8, resolution = 0.5)
Justine_sub_no8 <- RunUMAP(Justine_sub_no8, reduction = "pca", dims = 1:10)
DimPlot(Justine_sub_no8, reduction = "umap")
UMAPPlot(object = Justine_sub_no8, group.by ="orig.ident",cols = c('Control' = '#B40F20', 'N' = '#E58601','O'='#46ACC8','S'='#E2D200'))
UMAPPlot(object = Justine_sub_no8,cols = c(wes_palette("Darjeeling2", 9, type = c("continuous"))), label=T)
UMAPPlot(object = Justine_sub_no8, group.by ="orig.ident",split.by="orig.ident", cols = c('Control' = '#B40F20', 'N' = '#E58601','O'='#46ACC8','S'='#E2D200'))

Idents(object = Justine_sub_no8) <- Justine_sub_no8@meta.data[["seurat_clusters"]]
Justine_sub.markers <- FindAllMarkers(Justine_sub_no8, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Justine_sub.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
write.table(top10,"clusters_top10.txt")

saveRDS(Justine_sub.markers, file = "Justine_subno8.markers.RDS")

pdf(file='clusterHmap_noMicroG.pdf', width=20, height=10,bg="white")
DoHeatmap(Justine_sub_no8, features = top10$gene,group.colors =c(wes_palette("Darjeeling2", 9, type = c("continuous"))),angle=0,size = 6,hjust=0.8,group.bar.height=0.02) + 
  scale_fill_viridis_b("Expression",option = "E") + theme(axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(vjust=0.5, hjust=1)) 
dev.off()
