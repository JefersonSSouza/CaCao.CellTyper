# CaCao.CellTyper
A Seurat cell typer package.

![](header.png)

## Installation

### remotes package is required

```sh
install.packages('remotes')
library('remotes')
remotes::install_github('JefersonSSouza/CaCao.CellTyper')
```


## Usage example
### Run finAllMarkers from Seurat package
```sh
allmarkers <- findAllMarkers(object.seurat,only.pos=T, min.pct = 0.25, logfc.threshold = 0.25)
```
### Get top 50 variable gene
```sh
allmarkers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top.sig.genes #top 50 genes by cluster
```
### load CaCao.CellTyper and  dplyr packages
```sh
library(CaCao.CellTyper)
library(dplyr)
```
## View all tissues from cellMarkerDB
```sh
CaCao.Tissues()
```
 Returns a tissue list - please check the file on [CaCao.Tissues.txt](https://github.com/JefersonSSouza/CaCao.CellTyper/blob/master/CaCao.Tissues.txt)

### identifying all cell types an save a PDF file to path 
```sh
cluster.percentage <- Identify.CellTypes(
  all.markers.sig=top.sig.genes, #top 50 variales genes
  tissue='Lung', # tissue
  specie = 'Human', # Specie 'Human' or 'Mouse'
  cancer = 'Normal cell', # cancer 'Normal cell' or 'Cancer cell'
  path_to_save = getwd(), # path to save PDF
  plot_name = 'BARPLOT_CELLTYPES') # PDF name to save
```

### Generate a annotation into Seurat object
```sh
object.seurat <- seurat.CellTyper(object=object.seurat,
                 CaCao.cluster.percent = cluster.percentage,
                 rm.celltype = c('Alveolar macrophage'=2),  # remove cell types annotation attributed to the same percentage, in this case, will be removed 'Alveolar macrophage'.
                 remove.duplicated.cluster = T)
```

### Visualize the cell types
```sh
object.seurat$cell.type.CaCao
Idents(object.seurat)
DimPlot(object.seurat, reduction = "umap",label = TRUE, pt.size = 0.5)+NoLegend()
```



