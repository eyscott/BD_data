# Analysis of BD Rhapsody derived scRNA-seq data
These data come from astrocytes in culture (purified via immunopanning) and track the conversion of these astrocytes into various forms of oligodendrocytes after being pushed by the transcription factors Sox10, Olig2, and Nkx6.2

## General Workflow
The provided scripts go through alignment, demultiplexing and generation of genecount tables as well as analysis using seurat and slingshot trajectory analysis: \
1.Alignment of fastq files with Kallisto \
2.Generation of genecount tables with Kallisto (includes demultiplexing)\
3.Quality control and preparation of Seurat object \
4.Analysis of genetables using Seurat \
5.Analysis after removing cluster identified as microglia \
6.Trajectory analysis using Slingshot \
7.Analysis for Sox10 timepoint data in R (Figure 3) \
8.Cell Oracle analysis (Figure 3) 

## Datasets
Datatable will be made public after publication

### Publication:
Bajohr, J., Scott, E.Y., Olfat, A., Sadria, M., Lee, K., Fahim, M., Taha, H.T., Lozano Casasbuenas, D., Derham, A., Yuzwa,,S.A, Bader, G.D., & M. Faiz (2024) Direct lineage conversion of postnatal mouse cortical astrocytes to oligodendrocyte lineage cells. eLife. 13:RP98632 https://doi.org/10.7554/eLife.98632.1
