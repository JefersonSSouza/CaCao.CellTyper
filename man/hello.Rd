\name{Identify.CellTypes}
\alias{Identify.CellTypes}
\title{Identify.CellTypes}
\usage{
------------------------------------------
Intallation
------------------------------------------
install.packages('remotes')
library('remotes')
remotes::install_github('JefersonSSouza/CaCao.CellTyper')


------------------------------------------
Run finAllMarkers from Seurat package
------------------------------------------
  allmarkers <- findAllMarkers(object.seurat,only.pos=T, min.pct = 0.25, logfc.threshold = 0.25)


------------------------------------------
Get top 50 variable gene
------------------------------------------
  allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top.sig.genes #top 50 genes by cluster


------------------------------------------
load CaCao.CellTyper and dplyr packages
------------------------------------------
  library(CaCao.CellTyper)
  library(dplyr)


------------------------------------------
View all tissues from cellMarkerDB
------------------------------------------
  CaCao.Tissues()


------------------------------------------
Identifying all cell types an save a PDF file to path
------------------------------------------
  cluster.percentage <- Identify.CellTypes(
    all.markers.sig=top.sig.genes, #top 50 variales genes
    tissue='Lung', # tissue
    specie = 'Human', # Specie 'Human' or 'Mouse'
    cancer = 'Normal cell', # cancer 'Normal cell' or 'Cancer cell'
    path_to_save = getwd(), # path to save PDF
    plot_name = 'BARPLOT_CELLTYPES') # PDF name to sav


------------------------------------------
Generate a annotation into Seurat object
------------------------------------------
  object.seurat <- seurat.CellTyper(object=object.seurat,
                 CaCao.cluster.percent = cluster.percentage,
                 rm.celltype = c('Alveolar macrophage'),  # remove cell types annotation attributed to the same percentage, in this case, will be removed 'Alveolar macrophage'.
                 rm.cluster = c('2'), #cluster with  duplicated cell types annotation attributed to the same percentage.
                 remove.duplicated.cluster = T)


------------------------------------------
Visualize the cell types
------------------------------------------
  object.seurat$cell.type.CaCao
  Idents(object.seurat)
  DimPlot(object.seurat, reduction = "umap",label = TRUE, pt.size = 0.5)+NoLegend()
}
\description{
Identifying cell types from CellMarker DB.
https://github.com/JefersonSSouza/CaCao.CellTyper
}

