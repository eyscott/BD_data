read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}
control_mat <- read_count_output("/Users/erica_1/Desktop/Faiz/Justine/Justine_2/Control/", name = "control_counts")
S_mat <- read_count_output("/Users/erica_1/Desktop/Faiz/Justine/Justine_2/S/", name = "S_counts")
O_mat <- read_count_output("/Users/erica_1/Desktop/Faiz/Justine/Justine_2/O/", name = "O_counts")
N_mat <- read_count_output("/Users/erica_1/Desktop/Faiz/Justine/Justine_2/N/", name = "N_counts")

# Filter the matrix
dim(control_mat)#[1]  36047 232194
umi.per.barcode <- colSums(control_mat)
x <- sort(umi.per.barcode, decreasing = TRUE)
plot(x, log="xy",type="l", xlab="Barcodes", ylab="UMI counts") +abline(a=8000, b=0, h=8000, v=0,lty=2, col="red") 
control_mat_test <- control_mat[ , which(colSums(control_mat) >= 8000) ]
dim(control_mat_test)
#[1] 36047  1889
#head(control_mat_test)
umi.per.barcode <- colSums(S_mat)
x <- sort(umi.per.barcode, decreasing = TRUE)
plot(x, log="xy",type="l", xlab="Barcodes", ylab="UMI counts")+abline(a=4000, b=0, h=4000, v=0,lty=2, col="red")
S_mat_test <- S_mat[ , which(colSums(S_mat) >= 4000) ]
dim(S_mat_test)
#[1] 36047  1812

umi.per.barcode <- colSums(O_mat)
x <- sort(umi.per.barcode, decreasing = TRUE)
plot(x, log="xy",type="l", xlab="Barcodes", ylab="UMI counts")+abline(a=4000, b=0, h=4000, v=0,lty=2, col="red")
O_mat_test <- O_mat[ , which(colSums(O_mat) >= 4000) ]
dim(O_mat_test)
#[1] 36047  1654


umi.per.barcode <- colSums(N_mat)
x <- sort(umi.per.barcode, decreasing = TRUE)
plot(x, log="xy",type="l", xlab="Barcodes", ylab="UMI counts")+abline(a=4000, b=0, h=4000, v=0,lty=2, col="red")
N_mat_test <- N_mat[ , which(colSums(N_mat) >= 4000) ]
dim(N_mat_test)

tr2g <- read_tsv("transcripts_to_genes.txt", col_names = c("transcript", "gene", "gene_symbol")) %>%
  select(-transcript) %>%
  distinct()

# Convert from Ensembl gene ID to gene symbol
rownames(control_mat_test) <- tr2g$gene_symbol[match(rownames(control_mat_test), tr2g$gene)]
rownames(S_mat_test) <- tr2g$gene_symbol[match(rownames(S_mat_test), tr2g$gene)]
rownames(O_mat_test) <- tr2g$gene_symbol[match(rownames(O_mat_test), tr2g$gene)]
rownames(N_mat_test) <- tr2g$gene_symbol[match(rownames(N_mat_test), tr2g$gene)]

seu_control <- CreateSeuratObject(control_mat_test, min.cells = 3, project = "Control")
seu_S <- CreateSeuratObject(S_mat_test, min.cells = 3, project = "S")
seu_O <- CreateSeuratObject(O_mat_test, min.cells = 3, project = "O")
seu_N <- CreateSeuratObject(N_mat_test, min.cells = 3, project = "N")
