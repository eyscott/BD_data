sds <- slingshot(Embeddings(Justine_sub_no8, "umap"), clusterLabels = Justine_sub_no8$seurat_clusters, 
                 start.clus = 1, stretch = 0)
source_col<-c('Control' = '#B40F20', 'N' = '#E58601','O'='#46ACC8','S'='#E2D200')
cols_O<-source_col[Justine_sub_no8$orig.ident]
cols = c(wes_palette("Darjeeling2", 9, type = c("continuous")))
cols_C<-cols[Justine_sub_no8$seurat_clusters]

plot(sds@elementMetadata@listData[["reducedDim"]], col = cols_O, pch = 16, cex = 0.5,xlab = "UMAP 1", ylab = "UMAP 2")
lines(SlingshotDataSet(sds), lwd = 2, type = 'lineages', col = 'black')

#separate lineages and plot directionality
options(repr.plot.width=15, repr.plot.height=6)
nc <- 4
pt <- slingPseudotime(sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(sds@elementMetadata@listData[["reducedDim"]], col = colors, pch = 16, main = i)
  lines(SlingshotDataSet(sds), lwd = 2, col = 'black', type = 'lineages')
}

##random forest for genes driving trajectory analysis
top_hvg <- HVFInfo(Justine_sub_no8) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(variance.standardized)) %>% 
  top_n(300, variance.standardized) %>% 
  pull(bc)
# Prepare data for random forest
dat_use <- t(GetAssayData(Justine_sub_no8, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(sds)[,2], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])
library(tidymodels)
dat_split <- initial_split(dat_use_df)
dat_train <- training(dat_split)
dat_val <- testing(dat_split)
structure(dat_val)
#[1] 1122  301 (1122 cells, 300 genes)
#floor((300)/3)==100=mtry
#dim(dat_val)

model <- rand_forest(mtry = 100, trees = 500, min_n = 8, mode = "regression") %>%
  set_engine("ranger", importance = "impurity", num.threads = 3) %>%
  fit(pseudotime ~ ., data = dat_train)

val_results <- dat_val %>% 
  mutate(estimate = predict(model, .[,-1]) %>% pull()) %>% 
  select(truth = pseudotime, estimate)

var_imp <- sort(model$fit$variable.importance, decreasing = TRUE)
top_genes <- names(var_imp)[1:12]
pal <- viridis(100, end = 0.95)

pdf(file='lineage_genes_NoMicroglia.pdf', width=12, height=8,bg="white")
par(mfrow = c(3, 4))
for (i in seq_along(top_genes)) {
  colors <- pal[cut(dat_use[,top_genes[i]], breaks = 100)]
  plot(sds@elementMetadata@listData[["reducedDim"]], col = colors, pch = 16, cex=0.5,main = top_genes[i])
  lines(SlingshotDataSet(sds), lwd = 2, col = 'black', type = 'lineages')
}
dev.off()
