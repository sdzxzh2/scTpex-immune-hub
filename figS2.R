
library(vistributions)
library(copykat)
library(infercnv)
library(ggplot2)
library(ggpubr)
library(ggsci)

#figS2a,b
infercnv_obj = readRDS("run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$cCON
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- as.vector(unlist(test_loc))
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)
gn <- rownames(expr)
geneFile <- read.table("./all_sam_geneFile.txt", header=FALSE, row.names=NULL)
rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]
set.seed(20210418)
kmeans.result <- kmeans(t(expr), 7)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB=rownames(kmeans_df)
kmeans_df=kmeans_df%>%inner_join(anno.df,by="CB") #合并
kmeans_df_s=arrange(kmeans_df,kmeans_class) #排序
rownames(kmeans_df_s)=kmeans_df_s$CB
kmeans_df_s$CB=NULL
kmeans_df_s$kmeans_class=as.factor(kmeans_df_s$kmeans_class) #将kmeans_class转换为因子，作为热图的一个注释，最终在热图上就会按照1:7排序
head(kmeans_df_s)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
#top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:19,labels_gp = gpar(cex = 1.5)))
color_v=RColorBrewer::brewer.pal(8, "Dark2")[1:7] 
names(color_v)=as.character(1:7)
write.csv(kmeans_df_s,"kmeans_df_s_normal.csv")
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("normal"="#E29DCE"),
                                                     kmeans_class=color_v))
pdf("try_fin_normal.pdf",width = 15,height = 5)
ht = Heatmap(t(expr)[rownames(kmeans_df_s),], 
             col = colorRamp2(c(0.5,1,1.5), c("#377EB8","#F0F0F0","#E41A1C")), 
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"),
             
             heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.5,1,1.5),legend_height = unit(3, "cm")),
             
             top_annotation = top_anno,left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
draw(ht, heatmap_legend_side = "right")
dev.off()

#figS2c
data.combined=data.combined
df=data.combined@meta.data
A <- prop.table(table(df$fincell, df$group), margin = 2)
A <- as.data.frame(A)
colnames(A) <- c("celltype", "time", "Freq")
ggplot(A,aes(x = time,y =Freq,
             group=celltype))+
  stat_summary(geom = 'line',fun='mean',cex=1,col='white')+
  geom_area(data = A,aes(fill=celltype))+
  scale_fill_manual(values=celltype_col)+
  labs(x=NULL,y=NULL)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.text = element_text(color = "black",size = 10))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#figS2d
library(stringr)
library(Seurat)
library(miloR)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(qs)
library(BiocParallel)
register(SnowParam(workers =28, progressbar = TRUE)) 
scRNA_pre <- as.SingleCellExperiment(data.combined_fib)

scRNA_milo <- Milo(scRNA_pre)
reducedDim(scRNA_milo,"UMAP") <- reducedDim(scRNA_pre,"UMAP")
scRNA_milo
scRNA_milo <- buildGraph(scRNA_milo, k = 30, d = 15,reduced.dim = "PCA")
scRNA_milo <- makeNhoods(scRNA_milo, prop = 0.2, k = 30, d = 15, 
                         refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(scRNA_milo)
scRNA_milo <- countCells(scRNA_milo, 
                         meta.data = as.data.frame(colData(scRNA_milo)),
                         sample="orig.ident") 
traj_design <- data.frame(colData(scRNA_milo))[,c("orig.ident", "group")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
traj_design <- traj_design[colnames(nhoodCounts(scRNA_milo)), , drop=FALSE]
traj_design$group <- factor(traj_design$group)
traj_design
scRNA_milo <- calcNhoodDistance(scRNA_milo, d=15, reduced.dim = "PCA")

da_results <- testNhoods(scRNA_milo,design = ~ group, design.df = traj_design)
da_results %>%
  arrange(SpatialFDR) %>%
  head() 
ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) 

scRNA_milo <- buildNhoodGraph(scRNA_milo)

umap_pl <- plotReducedDim(scRNA_milo, dimred = "UMAP", 
                          colour_by="group", text_by = "fincell", 
                          text_size = 3, point_size=0.5) + guides(fill="none")
umap_pl
nh_graph_pl <- plotNhoodGraphDA(scRNA_milo, da_results,
                                layout="UMAP",alpha=0.9)+scale_fill_gradient2(low = "#1CC5FE",mid = "white",high = "#FB7D80")  #alpha默认0.1
nh_graph_pl
ggsave("nh_graph_pl2.pdf",  width =6, height = 5)
da_results <- annotateNhoods(scRNA_milo, 
                             da_results, 
                             coldata_col = "fincell")
head(da_results)
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)
str(da_results)
table(da_results$celltype)
range(da_results$SpatialFDR)
plotDAbeeswarm(da_results, group.by = "fincell",alpha = 0.9)+scale_color_gradient2(low = "#1CC5FE",mid = "white",high = "#FB7D80")
ggsave("miloR2.pdf",  width =6, height = 7)

#fig2b,g
some.sample.stat=as.data.frame(some.sample.stat)
rownames(some.sample.stat)=alltype
colnames(some.sample.stat)=sam
library(RColorBrewer)
library(scales)
library(pheatmap)
alltype_order=c()
for (i in 1:22) {
  for (j in c("p","q")) {
    for (k in c("gain","loss")) {
      alltype_order=append(alltype_order,paste("chr",i,j,"_",k,sep = ""))
    }
  }
}
alltype_order=alltype_order[alltype_order%in%rownames(some.sample.stat)]
some.sample.stat=some.sample.stat[alltype_order,]
pdf(paste0("heat_cnv1.pdf"),3.8,12)
pheatmap(some.sample.stat,
         cluster_rows = F,
         cluster_cols = F,
         color = colorRampPalette(brewer.pal(9,"PuRd"))(100),
         show_rownames=T,
         border_color = "grey"      
)
dev.off()
pdf(paste0("heat_cnv2.pdf"),4,12)
pheatmap(some.sample.stat,
         cluster_rows = T,
         cluster_cols = F,
         color = colorRampPalette(brewer.pal(9,"PuRd"))(100),
         show_rownames=T,
         border_color = "grey"
         
)
dev.off()

infercnv_obj = readRDS("./run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
expr<-expr[,test_loc]
expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
meta=data.combined_fib@meta.data
meta$CB=rownames(meta)
CNV_score=CNV_score%>%inner_join(meta,by="CB")
CNV_score%>%ggplot(aes(group,CNV_score))+geom_violin(aes(fill=group),color="NA")+
   scale_fill_manual(values = c("#1CC5FE","#6FC7CF","#FBA27D","#FB7D80"))+ylim(0,0.01)+
  theme_bw()+ 
  stat_compare_means(comparisons =my_comparisons, simplify =FALSE)
ggsave("CNV_group.pdf",width = 4.5,height = 5)

figS2e
#refer to fig3i

#figS2f
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)

scenicOptions <- initializeScenic(
  org="mgi", 
  dbDir="./cisTarget_databases", # RcisTarget databases location
  dbs="hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",# file name of motif database
  datasetTitle="SCENIC on Mouse Cell Atlas", # choose a name for your analysis
  nCores=1
)

genesKept <- geneFiltering(
  exprMat=dge, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)
dge <- dge[genesKept, ]
dim(dge)
if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")
source("utils/add_cellAnnotation.R")
saveLoom <- function(exprMat, output){
  cellInfo2 <- data.frame(
    row.names = rownames(cell.info),
    cellType =cell.info$fincell
  )
  loom <- build_loom(output, dgem=exprMat)
  loom <- add_cellAnnotation(loom, cellInfo2)
  close_loom(loom)
}
saveLoom(dge, "output/s1_avg20_rep1.loom")
loom <- build_loom("output/s1_exprMat.loom", dgem=dge)
loom <- add_cellAnnotation(loom, cell.info)
close_loom(loom)

rasMat <- readRDS("./output/s5_avg20_rep1.rasMat.rds")
s3_avg20_rep1.regulons <- read.delim("./output/s3_avg20_rep1.regulons.txt", header=FALSE, row.names=1)
rasMat2=t(rasMat)
rasMat_sub=rasMat2[,colnames(rasMat2)%in%rownames(data.combined@meta.data)]
row.names(rasMat_sub)=row.names(s3_avg20_rep1.regulons) 
fib.cell=colnames(rasMat_sub)
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined=sub_CT
data.combined$fincell=Idents(data.combined)
sigScores <-as.matrix(rasMat_sub)
Idents(data.combined) <- "group"
data.combined@assays$RNA@counts=sigScores
data.combined@assays$RNA@data=sigScores
data.combined@assays$RNA@scale.data=sigScores
data.combined.markers <- FindAllMarkers(data.combined, 
                                        only.pos = TRUE, 
                                        min.pct = 0.25,
                                        assay ="RNA" ,
                                        logfc.threshold = 0)

saveRDS(data.combined.markers, "degs.rds")
top10 <- data.combined.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
RotatedAxis2 <- function (...) 
{
  rotated.theme <- theme(axis.text.x = element_text(angle = 90, 
                                                    hjust = 1), validate = TRUE, ...)
  return(rotated.theme)
}
pdf("p.pdf",8,3.5)
DotPlot(data.combined, 
        assay ="RNA",
        features = unique(top10$gene))+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  RotatedAxis2()+scale_x_discrete("")+scale_y_discrete("")
dev.off()

#fig2h
#refer to fig4g

#fig2i
library(viridis)
plot_umap_density<-function(data,type,base_size=12,base_family=""){
  ggplot(data=data,aes(x=UMAP_1,y=UMAP_2))+
    theme_black(base_size,base_family)+
    xlim(c(min(data$UMAP_1)-2,max(data$UMAP_1)+2))+
    ylim(c(min(data$UMAP_2)-2,max(data$UMAP_2)+2))+
    stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE,alpha=1)+
    geom_point(color='grey90',size=0.05,alpha=0.4)+
    scale_fill_viridis(option="magma",alpha=1)+
    labs(title=paste0("DensityPlotfor",type))+
    theme(legend.position="none")
}
theme_black<-function(base_size=12,base_family=""){
  theme_grey(base_size=base_size,base_family=base_family)%+replace%
    theme(
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      legend.background=element_rect(color=NA,fill="black"),
      legend.key=element_rect(color="white",fill="black"),
      legend.key.size=unit(1.2,"lines"),
      legend.key.height=NULL,
      legend.key.width=NULL,
      legend.text=element_text(size=base_size*0.8,color="white"),
      legend.title=element_text(size=base_size*0.8,face="bold",hjust=0,color="white"),
      legend.position="none",
      legend.text.align=NULL,
      legend.title.align=NULL,
      legend.direction="vertical",
      legend.box=NULL,
      panel.background=element_rect(fill="black",color=NA),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.spacing=unit(0,"lines"),
      strip.background=element_rect(fill="grey30",color="grey10"),
      strip.text.x=element_text(size=base_size*0.8,color="white"),
      strip.text.y=element_text(size=base_size*0.8,color="white",angle=-90),
      plot.background=element_rect(color="black",fill="black"),
      plot.title=element_text(size=base_size*1.2,color="white"),
      plot.margin=unit(rep(0,4),"lines")
    )
}

p<-plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1.5,)
meta_monocle=p$data
meta_monocle$UMAP_1=meta_monocle$data_dim_1
meta_monocle$UMAP_2=meta_monocle$data_dim_2
unique_types<-unique(meta_monocle$group)
meta=meta_monocle
for (i in unique_types) {
  tissue_data<-meta[meta$group==i,]
  p<-plot_umap_density(tissue_data,i)
  ggsave(filename=paste0(i,"_density.pdf"),
         width=3,height=3,plot=p)
}

for (i in unique_types) {
  tissue_data<-meta[meta$group==i,]
  p<-ggplot(tissue_data,aes(x=UMAP_1,y=UMAP_2))+
    stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE)+
    # scale_fill_distiller(palette=4,direction=1)+
    scale_fill_gradient(low = "white", high = "darkblue") + 
    xlim(c(min(meta$UMAP_1)-2,max(meta$UMAP_1)+2))+
    ylim(c(min(meta$UMAP_2)-2,max(meta$UMAP_2)+2))+
    theme(
      legend.position='right'
    )
  ggsave(filename=paste0(i,"_density2.pdf"),
         width=5,height=4,plot=p)
}


