# Astrocyte to oligodendrocyte reprogramming scRNA-seq data (BD Biosciences)
Use of scRNAseq for transcription factor-mediated astrocyte to oligodendrocyte reprogramming.
scRNA-seq data generated using BD Rhapsody technology and processed using Kallisto package

## General Workflow
The provided scripts go through alignment, demultiplexing and generation of genecount tables as well as analysis using seurat and slingshot trajectory analysis: \
1.Alignment of fastq files with Kallisto \
2.Generation of genecount tables with Kallisto (includes demultiplexing)\
3.Quality control and preparation of Seurat object \
4.Analysis of genetables using Seurat \
5.Analysis after removing cluster identified as microglia \
6.Trajectory analysis using Slingshot 

## Publication
In Progress:  
Justine Bajohr, Erica Y. Scott, Arman Olfat, Kevin Lee, Maria Fahim, Hiba T. Taha, Daniela Lozano, Casasbuenas, Ann Derham, Scott A. Yuzwa & Maryam Faiz. "Astrocyte to oligodendrocyte reprogramming is transcription factor specific". 2022. Submitted
