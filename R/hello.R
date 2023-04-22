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



#save(tissues, marker_file, file = "mydata.rda")
CaCao.Tissues<- function(){
  return(tissues)
}

Identify.CellTypes <- function(all.markers.sig,specie='Human',tissue,cancer='Normal cell',path_to_save=getwd(),plot_name='bar_plot_celltypes'){
  require("ggplot2")
  require("dplyr")
  require("reshape2")
  require('RColorBrewer')
  require('patchwork')
  require('ggforce')
  require('ggrepel')

  celltypes.df <- data.frame(cluster=c(all.markers.sig$cluster),marker=c(all.markers.sig$gene))
  df.marker <- filter(marker_file,Species==specie)
  df.marker <- filter(marker_file,Cell.type==cancer)
  df.marker <- filter(marker_file,Tissue.class==tissue)
  merged.celltypes <- merge(celltypes.df, df.marker, by.x='marker',by.y='Marker')
  merged.celltypes$Cell.name <- as.factor(merged.celltypes$Cell.name)
  df_bar <- data.frame(cluster=c(merged.celltypes$cluster),celltype=c(merged.celltypes$Cell.name))
  df_bar <- melt(df_bar)
  clusters <- as.factor(c(0:length(unique(all.markers.sig$cluster))))
  pp_list<-list()
  cluster.cell.types<-list()
   df.final<-data.frame()
  for(clust in clusters){
    df_bar.2 <- df_bar[which(df_bar[,'cluster']==clust),]
    count.df<- df_bar.2 %>%
      dplyr::count(celltype)%>%
      dplyr::mutate(perc=n/sum(n)*100,
                    cluster=clust)
      count.df <-count.df[which(count.df[,"perc"]==max(count.df[,'perc'])),]
      df.final<-bind_rows(df.final,count.df)
    pp<- ggplot(df_bar.2, aes(x=celltype, y=cluster,fill=celltype)) +
      geom_bar(stat = "identity")+theme(
        axis.text.x = element_text(angle = 45,size = 8,face='bold',  vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18,face = 'bold'),
        legend.position = 'none'
      )
    pp<- pp+scale_fill_manual(values=colorRampPalette(brewer.pal(11,"Set1"))(200))
    pp_list[[clust]]<- pp
  }
  pp_list[[clust]]
  pp <- wrap_plots(pp_list,ncol=2)+plot_layout(heights=2, ncol=2)
  ggsave(paste(path_to_save,'/',plot_name,'.pdf',sep=''),plot =pp  ,width = 15, height = 45)
  print(paste('save file in : ',path_to_save,'/',plot_name,'.pdf',sep=''))
 
return(df.final)
}

seurat.CellTyper <- function(object,CaCao.cluster.percent,rm.celltype=c(), rm.cluster=c(),remove.duplicated.cluster=F){
  if(remove.duplicated.cluster==T){
    for(rm.indx in 1:length(rm.celltype)){
    print(paste('removing : ',rm.celltype[rm.indx],sep=''))
    CaCao.cluster.percent <- CaCao.cluster.percent[-which((CaCao.cluster.percent[,'celltype']==rm.celltype[rm.indx])& (CaCao.cluster.percent[,'cluster']==rm.cluster[rm.indx])) ,]
    
  }
  }
                

for(i in CaCao.cluster.percent$cluster){
  print(i)
  current_cluster= c(paste(i,sep=''))
  df <- as.data.frame(CaCao.cluster.percent[which(CaCao.cluster.percent['cluster']==i),])
  df$celltype <- as.character(df$celltype)
  new_cluster = c(df[1,'celltype'])
  object@active.ident <- plyr::mapvalues(x = object@active.ident, from = current_cluster, to = new_cluster)
  object$cell.type.CaCao <- Idents(object)
}
return(object)
}