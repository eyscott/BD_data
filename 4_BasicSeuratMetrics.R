Justine <- merge(seu_control, y = c(seu_S,seu_O,seu_N), add.cell.ids = c("Control","S","O","N"), project = "Justine")
Justine[["percent.mt"]] <- PercentageFeatureSet(Justine, pattern = "^mt-")

#Violin plots for QC
VlnPlot(Justine, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,cols =(wes_palette(n=4, name="Royal1")),
        pt.size=0.001, group.by ="orig.ident") 

Justine_sub <- subset(Justine, subset = nFeature_RNA > 200)
Justine_sub <- NormalizeData(Justine_sub) %>% ScaleData()
dim(Justine_sub)
#[1] 17236  7420
Justine_sub <- FindVariableFeatures(Justine_sub, nfeatures = 8618)#80%
top10 <- head(VariableFeatures(Justine_sub), 10)
plot1 <- VariableFeaturePlot(Justine_sub, log = FALSE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

ElbowPlot(Justine_sub, ndims = 20)

set.seed(42)
Justine_sub <- RunPCA(Justine_sub, features = VariableFeatures(object = Justine_sub))
Justine_sub <- FindNeighbors(Justine_sub, dims = 1:10)
Justine_sub <- FindClusters(Justine_sub, resolution = 0.5)
Justine_sub <- RunUMAP(Justine_sub, reduction = "pca", dims = 1:10)
DimPlot(Justine_sub, reduction = "umap")

UMAPPlot(object = Justine_sub, group.by ="orig.ident",cols = c('Control' = '#B40F20', 'N' = '#E58601','O'='#46ACC8','S'='#E2D200'))
UMAPPlot(object = Justine_sub,cols = c(wes_palette("Darjeeling2", 11, type = c("continuous"))), label=T)
UMAPPlot(object = Justine_sub, group.by ="orig.ident",split.by ="orig.ident",cols = c('Control' = '#B40F20', 'N' = '#E58601','O'='#46ACC8','S'='#E2D200'))
