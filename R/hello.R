# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'




CaCao.Tissues<- function(){
  require("readxl")
  require("stringr")
  tf <- tempfile(fileext=".xlsx")
  download.file("http://yikedaxue.slwshop.cn/Cell_marker_Seq.xlsx", tf, quiet=TRUE)
  rstudioapi::navigateToFile(tf)
  marker_file <- read_excel(tf)
  colnames(marker_file) <- str_replace(colnames(marker_file),' ','.')
  tissues <- c(unique(marker_file$Tissue.class))
  return(tissues)
}

Identify.CellTypes <- function(all.markers.sig,specie='Human',tissue,cancer='Normal cell',path_to_save=getwd(),plot_name='bar_plot_celltypes'){
  require("readxl")
  require("ggplot2")
  require("dplyr")
  require("reshape2")
  require('RColorBrewer')
  require('patchwork')
  require('ggforce')
  require('ggrepel')
  tf <- tempfile(fileext=".xlsx")
  download.file("http://yikedaxue.slwshop.cn/Cell_marker_Seq.xlsx", tf, quiet=TRUE)
  rstudioapi::navigateToFile(tf)
  marker_file <- read_excel(tf)
  celltypes.df <- data.frame(cluster=c(all.markers.sig$cluster),marker=c(all.markers.sig$gene))
  colnames(marker_file) <- str_replace(colnames(marker_file),' ','.')
  df.marker <- filter(marker_file,Species==specie)
  df.marker <- filter(marker_file,Cell.type==cancer)
  df.marker <- filter(marker_file,Tissue.class==tissue)
  merged.celltypes <- merge(celltypes.df, df.marker, by.x='marker',by.y='Marker')
  merged.celltypes$Cell.name <- as.factor(merged.celltypes$Cell.name)
  df_bar <- data.frame(cluster=c(merged.celltypes$cluster),celltype=c(merged.celltypes$Cell.name))
  df_bar <- melt(df_bar)
  clusters <- as.factor(c(0:31))
  pp_list<-list()

  for(clust in clusters){
    df_bar.2 <- df_bar[which(df_bar[,'cluster']==clust),]
    pp<- ggplot(df_bar.2, aes(x=celltype, y=cluster,fill=celltype)) +
      geom_bar(stat = "identity")+theme(
        axis.text.x = element_text(angle = 45,size = 8,face='bold',  vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18,face = 'bold'),
        legend.position = 'none'
      )
    pp<- pp+scale_fill_manual(values=colorRampPalette(brewer.pal(11,"Set1"))(31))
    pp_list[[clust]]<- pp
  }
  pp_list[[clust]]
  pp <- wrap_plots(pp_list,ncol=2)+plot_layout(heights=2, ncol=2)
  ggsave(paste(path_to_save,'/',plot_name,'.pdf',sep=''),plot =pp  ,width = 15, height = 45)
  print(paste('save file in : ',path_to_save,'/',plot_name,'.pdf',sep=''))

}
