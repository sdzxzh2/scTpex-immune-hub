library(dplyr)
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 


#fig4a
#refer to fig3a

#fig4b
library(Startrac)
library(ggplot2)
library(tictoc)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(sscVis)
meta=data.combined_fib@meta.data
celltype_sub=unique(meta$fincell)
celltype=celltype[celltype%in%celltype_sub]
meta$fincell=factor(meta$fincell,levels = celltype)
R_oe <- calTissueDist(meta,
                      byPatient = F,
                      colname.cluster = "fincell",
                      colname.patient = "orig.ident",
                      colname.tissue = "group",
                      method = "chisq", 
                      min.rowSum = 0) 
R_oe


col_fun  <-  colorRamp2(c(min(R_oe,na.rm  =  TRUE),1,max(R_oe,na.rm  =TRUE)),  
                        c("#f6f8e6",  "#f9a33e",  "red"))

Roe=R_oe
write.csv(Roe,"Roe.csv")
roe_df <- as.data.frame(as.table(as.matrix(Roe)))
colnames(roe_df)<-c("CellType","Tissue","Roe")
roe_df <-roe_df%>%mutate(
  Enrichment = ifelse(Roe >= 1, "Enrichment", "Depletion"), 
  Roe = ifelse(Roe == 0, NA, Roe)  
)  %>% 
  filter(!is.na(Roe))  
roe_df$CellType=factor(roe_df$CellType,levels = rev(celltype))
p <- ggplot(roe_df, aes(x = Tissue, y = CellType)) +
  coord_flip() +
  geom_point(
    aes(size = Roe, color = Enrichment)
  ) +
  scale_size_continuous(
    name = "Ro/e",
    breaks = c(0.5, 1.0, 1.5),
    range = c(1,7) 
  ) +
  scale_color_manual(
    name = "Status",
    values = c("Enrichment" = "#E36C8C", "Depletion" = "#2E75B6"),  
    #labels = c("Enrichment", "Depletion"),
    guide = guide_legend(  
      override.aes = list(
        size = 4  
      )
    )
  ) +
  scale_y_discrete(limits = rev) + 
  labs(
    title = "Cell Type Enrichment",
    x = "Tissue Group",
    y = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.box = "horizontal" 
  )
ggsave("roe.pdf", plot = p, width =6, height = 3.5)


#fig4d
#refer to figS1d


#fig4f,e

sigScores <- as.matrix(t(getSignatureScores(vis)))
data=as.data.frame(t(sigScores))
id=data.combined_fib@meta.data$name
type=factor(data.combined_fib@meta.data$fincell,levels=celltype_sub)
group=as.vector(data.combined_fib@meta.data$group)
for (i in gene_use) {
  Signature_score=data[i]
  Signature_score=Signature_score[,1]
  pdata_melt=data.frame(id,group,type,Signature_score)
    c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = group )) +
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")+stat_compare_means(label = "p.signif")  
  ggsave(paste0("cell&group1_",i,".pdf"),plot = c, width=14, height=5)
  
  c<-ggplot(pdata_melt, aes(x = type, y = Signature_score, fill = group)) +
    geom_violin( aes(color=group),width=.8)+
    scale_fill_manual(values = sample_color)+
    scale_color_manual(values = sample_color)+
    theme_classic() +
    geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means(label = "p.signif")
  ggsave(paste0("cell&group2_",i,".pdf"),plot = c, width=20, height=5)
   c <- ggplot(pdata_melt, aes(x = group, y = Signature_score, fill = group)) +
    geom_violin( aes(color=group))+
    scale_fill_manual(values = sample_color)+
    scale_color_manual(values = sample_color)+
    theme_classic() +
    geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means()
  p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
  ggsave(paste0("group2_",i,".pdf"),plot = p, width=3.5, height=5)
  c <- ggplot(pdata_melt,
              aes(x=type, y=Signature_score, 
                  fill = type, 
              )) + 
    geom_boxplot(notch = F, alpha = 0.95, 
                 outlier.shape = 16,
                 outlier.colour = "black", 
                 outlier.size = 0.65) +
    scale_fill_manual(values = celltype_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")  
  ggsave(paste0("cell1_",i,".pdf"),plot = c, width=8, height=5)
  
  c <- ggplot(pdata_melt, aes(x = type, y = Signature_score, fill = type)) +
    geom_violin( aes(color=type))+
    scale_fill_manual(values = celltype_col)+
    scale_color_manual(values = celltype_col)+
    theme_classic() +
    geom_boxplot( aes(color=type),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means()
  ggsave(paste0("cell2_",i,".pdf"),plot = c, width=9, height=4)
    
  
}


male <- data.combined_1
femlale <- data.combined_2
male_DGEs <- list()
cells1 <- unique(male$fincell)
table(mlale$fincell,mlale$group)
for (i in 1:length(unique(male$fincell))) {
  data = subset(male, fincell==cells1[i])
  df <- FindMarkers(data,
                    group.by="group",
                    ident.1="B_ca",
                    ident.2="A_ca",
                    logfc.threshold = 0,
                    min.cells.group = 1,
                    min.pct = 0)
  male_DGEs[[i]] <- df
}

names(male_DGEs) <- cells1
femlale_DGEs <- list()
cells2 <- unique(femlale$fincell)
table(femlale$fincell,femlale$group)
for (i in 1:length(unique(femlale$fincell))) {
  data = subset(femlale, fincell==cells2[i])
  df <- FindMarkers(data,
                    group.by="group",
                    ident.1="B_lymph",
                    ident.2="A_lymph",
                      min.cells.group = 1,
                    logfc.threshold = 0,
                    min.pct = 0)
  femlale_DGEs[[i]] <- df
}

names(femlale_DGEs) <- cells2
for(i in seq_along(male_DGEs)){
  male_DGEs[[i]]$sex <- "M"
  male_DGEs[[i]]$cluster <- names(male_DGEs)[i]
  male_DGEs[[i]]$gene <- rownames(male_DGEs[[i]])
}
for(i in seq_along(femlale_DGEs)){
  femlale_DGEs[[i]]$sex <- "F"
  femlale_DGEs[[i]]$cluster <- names(femlale_DGEs)[i]
  femlale_DGEs[[i]]$gene <- rownames(femlale_DGEs[[i]])
}

myorder <- celltype_sub
male_DGEs <- male_DGEs[myorder]
femlale_DGEs <- femlale_DGEs[myorder]
for(i in seq_along(male_DGEs)){
  male_DGEs[[i]] <- male_DGEs[[i]][which(male_DGEs[[i]]$p_val < 1),]
  male_DGEs[[i]] <- male_DGEs[[i]][order(male_DGEs[[i]]$avg_log2FC, decreasing = TRUE), ]
}


for(i in seq_along(femlale_DGEs)){
    femlale_DGEs[[i]] <- femlale_DGEs[[i]][which(femlale_DGEs[[i]]$p_val < 1),]
  femlale_DGEs[[i]] <- femlale_DGEs[[i]][order(femlale_DGEs[[i]]$avg_log2FC, decreasing = TRUE), ]
  
}

top_male <- list()
for(i in seq_along(male_DGEs)){
  top_male[[i]] <- rbind(head(male_DGEs[[i]], 5), tail(male_DGEs[[i]],15))
}
top_female <- list()
for(i in seq_along(femlale_DGEs)){
  top_female[[i]] <- rbind(head(femlale_DGEs[[i]], 5), tail(femlale_DGEs[[i]],15))  
}

 male_DEGs_file <- do.call(rbind, male_DGEs)
 femlale_DGEs_file <- do.call(rbind, femlale_DGEs)
 write.csv(male_DEGs_file, file = "male_DEGs_file.csv")
 write.csv(femlale_DGEs_file, file = "femlale_DGEs_file.csv")
top_female <- do.call(rbind, top_female)
top_female <- top_female[,c(2,6:8)]
top_female$group <- paste0(top_female$sex,"_",top_female$cluster)
top_female <- top_female[,-c(2:3)]
top_female<-spread(top_female, group, avg_log2FC)
rownames(top_female) <- top_female$gene
top_female <- top_female[,-1]
library(ComplexHeatmap)
library(circlize)

cell=gsub("^..","",colnames(A))

annotation_col = data.frame(
  group = c(rep("B_ca.vs.A_ca",length(unique(male$fincell))),rep("B_lymph.vs.A_lymph",length(unique(male$fincell)))),
  celltype = c(cell))

row.names(annotation_col) <- colnames(A)
groupcolor <- c("#2f5688","#CC0000")

celltypecolor <- celltype_col

ann_colors <- list(group=groupcolor, celltype=celltypecolor) 
col_fun=colorRamp2(c(-2,0,2),c("#2f5688", "white","#CC0000"))

pdf("heat_two_compare_1.pdf",5,4.5)
ComplexHeatmap::pheatmap(as.matrix(A), 
                         scale = "none",
                         show_rownames = F,
                         show_colnames = F,
                         #col = colorRampPalette(c("#2f5688", "white","#CC0000"))(10),
                         col =col_fun,
                         column_title_gp  = gpar(fontsize = 10),
                         fontsize_row   =  gpar(fontsize = 8),
                         row_title = NULL,
                         cluster_rows = T，
                         cluster_cols = F,
                         column_split = annotation_col$group,
                         heatmap_legend_param = list(
                           title='Log2FC'),
                         annotation_col  = annotation_col,
                         annotation_colors = ann_colors,
                         annotation_legend = T,
                         annotation_names_col = F,
                         border_color = NA
                         )
dev.off()

#fig4g
library(DESeq2)
library(Seurat)
library(IHW)
library(monocle)
library(tidyverse)
library(magrittr)
library(Seurat) 
library(pheatmap)
library(RColorBrewer)
library(aplot)
Sys.setenv(LANGUAGE = "en") 
options(stringsAsFactors = FALSE) 

p1 <- plot_cell_trajectory(monocle_cds, color_by = "group", cell_size = 2)

p2 <- plot_cell_trajectory(monocle_cds, color_by = "fincell", cell_size = 2)

cowplot::plot_grid(p2,p1, ncol = 2)
ggsave(paste0('psedu2','.pdf'), width=8, height=5)
dev.off()

plot_cell_trajectory(monocle_cds, color_by = "clone", cell_size = 0.9,alpha=0.9)+
  scale_colour_manual(values = c("#F05662","grey90"))

#fig4i
library(Seurat)
library(scRepertoire)

cell.info=data.combined@meta.data
cell.info$clone_group=paste0(cell.info$clone,"_",cell.info$group)
cell.info$UMAP_1=data.combined@reductions$umap@cell.embeddings[,1]
cell.info$UMAP_2=data.combined@reductions$umap@cell.embeddings[,2]
cell.info$fincell=as.vector(cell.info$fincell)
col2=c("grey90")
names(col2)="Others"
celltype_col_raw=celltype_col
celltype_col=c(celltype_col,col2)
umapPlot.tissue <- function(data, tissue.use, topn) {
  clusters.use <- table(subset(data, Tissue == tissue.use)$ClusterID)
  if(topn > length(clusters.use)) {
    topn <- length(clusters.use)
  }
  clusters.use <- names(head(sort(clusters.use, decreasing = T), n=topn))
  data$new.cluster <- ifelse(data$Tissue == tissue.use & data$ClusterID %in% clusters.use, data$CellType, "Others")
  data$new.cluster <- factor(data$new.cluster, levels = c(setdiff(names(table(data$new.cluster)), "Others"), "Others"))
  data$pt.size <- ifelse(data$new.cluster == "Others", 0.4, 0.4)
  data1=data[data$new.cluster=="Others",]
  data2=data[data$new.cluster!="Others",]
  ggplot() + 
    geom_point(data = data1,
               aes(UMAP_1, UMAP_2, color=new.cluster),
               size=data1$pt.size,
               alpha=0.4) + 
    geom_point(data =data2,
               aes(UMAP_1, UMAP_2, color=new.cluster),
               size=data2$pt.size,
               alpha=0.5)+
    #scale_color_manual(values = c(pal_d3("category10")(topn), "grey")) + 
    scale_color_manual(values =celltype_col) + 
    guides(colour = guide_legend(keyheight = .7, 
                                 keywidth = .1, 
                                 override.aes = list(size=3))) + 
    theme_bw(base_size = 12) + 
    ggtitle(tissue.use) + 
    theme(#legend.justification = c(0,0), 
      legend.title = element_blank(),
      legend.key = element_rect(fill = alpha("white", 0)),
      legend.background = element_rect(fill=alpha('white', 0)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color="black"),
      plot.title = element_text(hjust = .5, face = "bold")
    )
  
}

cell.info$clone_sam=cell.info$clone_group
cell.info$Tissue=cell.info$clone_group
cell.info$ClusterID=cell.info$fincell
cell.info$CellType=cell.info$fincell
no_clone_group=paste0("no_clone","_",unique(cell.info$group))
clone_group=paste0("clone","_",unique(cell.info$group))
for (i in c(no_clone_group,clone_group)) {
  p<-umapPlot.tissue(cell.info, i, 25)
  ggsave(paste0(i,"_type1.pdf"),plot = p,width = 6.5,height = 5)
  
  
}

celltype_col=rep("#3A7BBF",length(celltype))
names(celltype_col)=celltype

celltype_col=c(celltype_col,col2)

for (i in c(no_clone_group)) {
  p<-umapPlot.tissue(cell.info, i, 25)
  ggsave(paste0(i,"_type2_A_no_clone.pdf"),plot = p,width = 6.5,height = 5)
  
  
}

celltype_col=rep("#E6315E",length(celltype))
names(celltype_col)=celltype
celltype_col=c(celltype_col,col2)

for (i in c(clone_group)) {
  p<-umapPlot.tissue(cell.info, i, 25)
  ggsave(paste0(i,"_type2_B_clone.pdf"),plot = p,width = 6.5,height = 5)
  
  
}
