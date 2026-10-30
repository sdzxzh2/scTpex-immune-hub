
library(Seurat,lib.loc="~/rpackage2")
packageVersion("Seurat")
library(SeuratObject)
library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)
library(harmony)
library(arrow)
library(scran)
library(data.table)
library(harmony)
library(easySCFr)
library(png)

#data prepare
expr_raw=dir(path ="./expr_fin/")
tmp= data_raw <- data.table::fread(file = paste0("./expr_fin/",expr_raw[1]), #sep = sep,
                                   data.table = FALSE, 
                                   verbose = FALSE)
name=colnames(tmp)
lis=list()
for (i in expr_raw) {
  
  
  data_raw <- data.table::fread(file = paste0("./expr_fin/",i), #sep = sep,
                                data.table = FALSE, 
                                verbose = FALSE)
  lis[[i]]=sum(colnames(data_raw)%in%name)
  
}

lambda <- c(0.2)
k_geom <- c(30)
seu_list=list()
spe_list=list()
compute_agf <- TRUE
SEED=1000
use_agf<-TRUE
npcs <- 20
name1=paste0("PCA_M1_lam",lambda)
name2=paste0("UMAP_M1_lam",lambda)
file=dir(path ="./expr_fin/")
allsam=strsplit(file,"\\.")
allsam=lapply(allsam,FUN = function(x){x[1]})
allsam=unlist(allsam)
allsam=gsub("fin_cell_mean_", "", allsam)

for (sam in allsam) {
  
    sample_raw <- LoadAkoya( paste0("./expr_fin/fin_cell_mean_",sam,".csv"),
                           type="processor",
                           fov=sam)
  colnames(sample_raw)=paste0(sam,"_",1:ncol(sample_raw))
  sample_raw@assays$Spatial=sample_raw@assays$Akoya
  sample_raw$sample_id=sam
   image <- readPNG(source = paste0("./geojson/",sam,".geojson.png"))
  load("./scale.factors.Rdata")
  fov=sample_raw@images[[1]]
  scale.factors$hires=dim(image)[1]/max(fov$centroids@coords[,1])
  scale.factors$lowres=dim(image)[1]/max(fov$centroids@coords[,1])
  
  visium.v2 <- new(Class = "VisiumV2", boundaries = fov@boundaries, 
                   molecules = fov@molecules, assay = fov@assay, key = fov@key, 
                   image = image, scale.factors = scale.factors)
  sample_raw@images[[sam]]<-visium.v2
  sample_raw@assays[["Spatial"]]@key<-"spatial_"
  saveH5(sample_raw, paste0(sam,"_step1_python.h5"),assay = "Spatial")
  qs::qsave(sample_raw, paste0(sam,"_step1_seurat.qs"))
  sample1=sample_raw
  gcm <- GetAssay(sample1)
  gcm=gcm$counts
  locs <- as.matrix(GetTissueCoordinates(sample1)[,1:2])
  sample1 <- SpatialExperiment(assay = list(counts = gcm), spatialCoords = locs)
  sample1$sample_id=sam
  qs::qsave(sample1, paste0(sam,"_step1_spe.qs"))
  sample_clean=sample_raw
  
  seu_list[[sam]]=sample_clean
  spe_list[[sam]]=sample1
  
}

SpatialDimPlot(sample_raw,pt.size.factor = 0.1)
spe=qs::qread(paste0(allsam[1],"_step1_spe.qs"))
for (sam in allsam[2:length(allsam)]) {
  tmp=qs::qread(paste0(sam,"_step1_spe.qs"))
  spe=cbind(spe,tmp)
}

qs::qsave(spe, paste0("all_sam_spe_mihc.qs"))
locs <- spatialCoords(spe)
locs <- cbind(locs, sample_id = factor(spe$sample_id))
locs_dt <- data.table(locs)
colnames(locs_dt) <- c("sdimx", "sdimy", "group")
locs_dt[, sdimx := sdimx + group * global_max]
locs <- as.matrix(locs_dt[, 1:2])
rownames(locs) <- colnames(spe)
spatialCoords(spe) <- locs

ST<-merge(seu_list[[1]], y=seu_list[2:length(seu_list)]) 
ST <- FindVariableFeatures(ST, nfeatures = 2000)
VariableFeatures(ST)<-rownames(ST)
aname <- "normcounts"
assay(spe, aname) <- normalizeCounts(spe, log = FALSE)
spe <- spe[VariableFeatures(ST),]

qs::qsave(spe, paste0("ST_raw_mihc.qs"))
spe <- computeBanksy(spe, assay_name = aname, compute_agf = compute_agf, k_geom = k_geom)
spe <- runBanksyPCA(spe, use_agf = use_agf, lambda = lambda, npcs = npcs, seed = SEED)

set.seed(SEED)
harmony_embedding <- HarmonyMatrix(
  data_mat = reducedDim(spe, name1),
  meta_data = colData(spe),
  vars_use = c("sample_id"),
  do_pca = FALSE,
  max.iter.harmony = 20,
  verbose = FALSE
)


reducedDim(spe, "Harmony_BANKSY") <- harmony_embedding
spatial_plots <- lapply(unique(spe$sample_id), function(snm) {
  x <- spe[, spe$sample_id == snm]
  ari <- aricode::ARI(x$clust_Harmony_BANKSY_k30_res0.7, colData(x)[, cnm])
  
  df <- cbind.data.frame(clust=colData(x)[[cnm]], spatialCoords(x))
  ggplot(df, aes(x=sdimy, y=sdimx, col=clust)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = pals::kelly()[-1]) +
    theme_classic() + 
    theme(
      legend.position = "none",
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank()) +
    labs(title = sprintf("Sample %s - ARI: %s", snm, round(ari, 3))) +
    coord_equal()
})


for (i in 1:length(spatial_plots)) {
  
  sam=unique(spe$sample_id)[i]
  p=spatial_plots[[i]]
  factor=30/(max(p$data$sdimx)-min(p$data$sdimx))
  y=(max(p$data$sdimy)-min(p$data$sdimy))*factor
  
  pdf(paste0("Cluster_",sam,"_type1.pdf"),30,30)
  plot(p)
  dev.off()
  
}

Idents(ST)<-"Harmony_BANKSY"
names(col)=unique(ST$Harmony_BANKSY)


allsam=unique(ST$sample_id)
for (sam in allsam) {
  
  tmp=ST@meta.data
  
  fib.cell <- rownames(tmp[tmp$sample_id==sam,])
  sub_CT<- subset(ST, name %in% fib.cell)
  sub_CT$fincell=Idents(sub_CT)
  p<-SpatialDimPlot(sub_CT,group.by='Harmony_BANKSY',ncol =1,label=F,cols = col,pt.size.factor = 0.11,image.alpha = 0)
  tmp=p[[1]]
  factor=10/(max(tmp$data$y)-min(tmp$data$y))
  y=(max(tmp$data$x)-min(tmp$data$x))*factor
  
  ggsave(plot = p,paste0("cluster_",sam,"_type2.pdf"),width = 10,height = y,limitsize = FALSE)
  
}

ST@assays[["Spatial"]]@key<-"spatial_"
qs::qsave(ST, "ST_mihc_allsam.qs")
qs::qsave(spe, "spe_mihc_allsam.qs")

ST[["Akoya"]]<-JoinLayers(ST[['Akoya']])
DefaultAssay(ST)<-"Akoya"
ST <-NormalizeData(ST)
ST <- FindVariableFeatures(ST)
VariableFeatures(ST)<-rownames(ST)
ST <- ScaleData(ST)

ST[["Spatial"]]<-JoinLayers(ST[['Spatial']])
DefaultAssay(ST)<-"Spatial"
ST <-NormalizeData(ST)
ST <- FindVariableFeatures(ST)
VariableFeatures(ST)<-rownames(ST)
ST <- ScaleData(ST)
qs::qsave(ST, "ST_mihc_allsam_JoinLayers.qs")

library(dplyr) 
library(presto)
presto <- wilcoxauc(ST,seurat_assay = "Spatial",group_by = "Harmony_BANKSY") 
presto_top5 <- presto %>%
  filter(logFC > 0 & pct_in > 25 & padj < 0.05) %>%
  group_by(group) %>%
  arrange(desc(logFC), .by_group=T) %>%
  top_n(n = 5, wt = logFC) 
presto_top5
write.csv(presto_top5,"degs.csv")

gene=rownames(ST)
library(tidyverse)

cellratio1 <- cellratio%>%
  group_by(celltype, group) %>%
  summarise(avg.freq = mean(Freq),
            sd.freq = sd(Freq),
            se.freq = sd(Freq)/sqrt(n()))
ggplot(data = cellratio1, mapping = aes(x = celltype, y = avg.freq,
                                        color =group, fill = group)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(data = cellratio1,
                aes(ymin=avg.freq-se.freq, ymax=avg.freq+se.freq),
                width=0.3,
                size = 0.5,
                color = "black",
                position=position_dodge(0.9)) +
  geom_jitter(data = cellratio,
              mapping = aes(x = celltype, y = Freq),
              shape=21,alpha=0.8,
              color = "black",
              position = position_dodge(width = 0.9),
              size = 1) 
ggsave('box_all2.pdf', width = 13, height = 4.6)

ggplot(data = cellratio1, mapping = aes(x = celltype, y = avg.freq,
                                        color =group, fill = group)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(data = cellratio1,
                aes(ymin=avg.freq-se.freq, ymax=avg.freq+se.freq),#加误差线
                width=0.3,
                size = 0.5,
                color = "black",
                position=position_dodge(0.9)) 
ggsave('box_all4.pdf', width = 10, height = 4.6)

library(ggpubr)

for (i in unique(cellratio$celltype)) {
  
  cellratio_sub=cellratio[cellratio$celltype==i,]
  cellratio1_sub <- cellratio_sub%>%
    group_by(celltype, group) %>%
    summarise(avg.freq = mean(Freq),
              sd.freq = sd(Freq),
              se.freq = sd(Freq)/sqrt(n()))
  ebtop<-function(x){
    return(mean(x)+sd(x)/sqrt(length(x)))
  }
  ebbottom<-function(x){
    return(mean(x)-sd(x)/sqrt(length(x)))
  }
  
  c <- ggplot(data=cellratio_sub,aes(x=group,y=Freq,fill=group))+
    stat_summary(geom = "bar",
                 fun = mean,
    )+
    stat_summary(geom = "errorbar",
                 fun.min = ebbottom,
                 fun.max = ebtop,
                 width=0.2)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10))  + 
    geom_jitter(aes(fill=group),width = 0.2,size=1,shape=21,color="black") + 
    stat_compare_means(comparisons =my_comparisons, simplify =FALSE)
  ggsave(plot = c,paste0("box1_",i,".pdf"), width=3.8, height=4.3)
  
  
  p <- ggplot(cellratio_sub,
              aes(x=group, y=Freq, 
                  color=group
              )) + 
    geom_boxplot(position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge()) +
    theme_bw()+  
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10))  + 
    stat_compare_means(comparisons =my_comparisons, simplify =FALSE)
  ggsave(plot = p,paste0("box2_",i,".pdf"), width=3.8, height=4)
}
ann <- read.delim("./ann.txt", row.names=NULL)

new.cluster.ids<-ann$name
new.cluster.ids
names(new.cluster.ids)<-levels(ST)
ST<-RenameIdents(ST,new.cluster.ids)
ST@meta.data$fincell<-ST@active.ident
ST@meta.data$celltype<-ST@meta.data$fincell
ST@meta.data$name<-rownames(ST@meta.data)
qs::qsave(ST, "ST_ann_fin.qs")
library(easySCFr)
allsam=unique(ST$sample_id)
for (sam in allsam) {
  tmp=ST@meta.data
  fib.cell <- rownames(tmp[tmp$sample_id%in%sam,])
  sub_CT<- subset(ST, name %in% fib.cell)
  sub_CT$fincell=Idents(sub_CT)
  
  
  saveH5(sub_CT, paste0("mIHC_ann_",sam,".h5"),assay = "Spatial")
  
}
ST=qs::qread( "ST_ann_fin.qs")

#fig5a
for (sam in allsam) {
  tmp=ST@meta.data

  fib.cell <- rownames(tmp[tmp$sample_id%in%sam,])
  sub_CT<- subset(ST, name %in% fib.cell)
  sub_CT$fincell=Idents(sub_CT)
  p<-SpatialDimPlot(sub_CT,label=F,cols = cols,pt.size.factor = 0.2,image.alpha = 0)
  tmp=p[[1]]
  factor=20/(max(tmp$data$y)-min(tmp$data$y))
  y=(max(tmp$data$x)-min(tmp$data$x))*factor
  ggsave(plot = p,paste0("ann_",sam,"_type2.pdf"),width = 20,height = y,limitsize = FALSE)
  
}
#fig5b,c
library(stringr)
library(tidyverse)
library(ggplot2)

celltype_main_sort=c("Lymphocyte","Myeloid","Stromal","Malignant_total")
ST$celltype_main=cell
ST$celltype_main=factor(ST$celltype_main,levels = celltype_main_sort)
qs::qsave(ST, "ST_ann_fin3_with_main_type.qs")
data_fin=ST@meta.data
data_for_ratio2=data_fin[,c("sam_group","celltype_main")]
count2=prop.table(table(data_for_ratio2$celltype_main, data_for_ratio2$sam_group), margin = 2)
count2=as.data.frame(count2)
colnames(count2)<-c("cell","sam_group","Freq")

count2$sam_group2=count2$sam_group
count2 = count2 %>% 
  separate(sam_group2, into = c('sam','group'), sep = ':')
count_main=count2
count_sub=count2[count2$cell=="Malignant_total",]
count_sub=count_sub[order(count_sub$Freq,decreasing = T),]
count_sub=count_sub[order(count_sub$group),]
sam_sort=as.vector(count_sub$sam_group)
count2$sam_group=factor(count2$sam_group,levels = sam_sort)
count2$cell=factor(count2$cell,levels =celltype_main_sort)

p3<-ggplot(data = count2, aes(x = sam_group, y = Freq, fill = cell)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=c("grey80","grey60","grey40","grey10")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  rotate_x_text()+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_line(color =  NA))
p3

data_fin=ST@meta.data
data_for_ratio2=data_fin[,c("sam_group","celltype")]
count2=prop.table(table(data_for_ratio2$celltype, data_for_ratio2$sam_group), margin = 2)
count2=as.data.frame(count2)
colnames(count2)<-c("cell","sam_group","Freq")
count2$sam_group2=count2$sam_group
count2 = count2 %>% 
  separate(sam_group2, into = c('sam','group'), sep = ':')
count2$sam_group2
count_sub=count2[count2$cell=="Malignant",]
count_sub=count_sub[order(count_sub$Freq,decreasing = T),]
count_sub=count_sub[order(count_sub$group),]
count2$sam_group=factor(count2$sam_group,levels = sam_sort)
count2$cell=factor(count2$cell,levels =celltype)
count_all_type=count2

p1<-ggplot(data = count2, aes(x = sam_group, y = Freq, fill = cell)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=celltype_col) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  rotate_x_text()+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))
p1
data_for_ratio1=data_fin[,c("sam_group","celltype")]
count=as.data.frame(table(data_for_ratio1$sam_group))
count$sam_group=count$Var1
count = count %>% 
  separate(Var1, into = c('sam','group'), sep = ':')

count$sam_group=factor(count$sam_group,levels = sam_sort)

p2<-ggplot(data=count)+
  geom_bar(aes(x=sam_group, y=Freq, fill=group), stat='identity')+
  scale_fill_manual(values=sample_color) +theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="count")+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_line(color =  NA))

p2 

library(patchwork)
p <- p2 + p3 +p1 + plot_layout(ncol = 1, heights = c(1.2, 2,2))
p

ggsave("ratio_stack.pdf",width = 8,height = 8)
library(ggplot2)
count_all=count_main
for (i in as.vector(unique(count_all$cell))) {
  

  count_sub=count_all[count_all$cell==i,]
  p4 <- ggplot(count_sub, aes(x = group, y = Freq, color = group)) +
    geom_jitter(size = 7, width = 0.15, show.legend = FALSE,alpha=0.6) + 
    #geom_point(position = 'jitter', size = 1.5, show.legend = FALSE) +
    scale_color_manual(values = sample_color) +
    theme(panel.grid = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(color = 'black')) +
    labs(x = 'Group', y = 'Cell Freq') +
    stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, 
                 geom = 'crossbar', width = 0.3, size = 0.25) +
    stat_summary(fun.data = function(x) median_hilow(x, 0.5), 
                 geom = 'errorbar', width = 0.25, size = 0.4)+ 
    stat_compare_means(label = "p.signif",comparisons = my_comparisons,method = "t.test")
   ggsave(filename = paste0("ratio_",i,".pdf"),plot = p4,width = 3.8,height = 4.3)
  
}


#fig5d
library(pROC)
library(ROCR)

roc_input=rbind(count_main,count_all_type)
roc_input$group2=roc_input$group


roc_input = roc_input %>% 
  separate(group2, into = c('response','treat'), sep = '_')
auc=c()
for (i in unique(roc_input$cell)) {
  roc_input_sub=roc_input[roc_input$cell%in%i,]
  tmp=roc(roc_input_sub$response,roc_input_sub$Freq)
  tmp=tmp$auc
  names(tmp)=i
  auc=c(auc,tmp)
  
}

plot_dat=as.data.frame(auc)
plot_dat$cell=rownames(plot_dat)
plot_dat=plot_dat[order(plot_dat$auc),]
plot_dat$cell=factor(plot_dat$cell,levels = rownames(plot_dat))

ggplot(plot_dat, aes(x = auc, y = cell,fill = auc,color=auc)) +
  geom_col( width = 0.7) +
  scale_fill_gradient(low = "grey90", high = "grey9")+
  scale_color_gradient(low = "grey90", high = "grey9")+
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  labs(
    x = "AUC",
    y = NULL    ) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0("auc_all_sam_response.pdf"),width = 4.5,height = 6)

auc=c()
roc_input_pre=roc_input[roc_input$treat%in%c("before"),]

for (i in unique(roc_input$cell)) {
  
  roc_input_sub=roc_input_pre[roc_input_pre$cell%in%i,]
  tmp=roc(roc_input_sub$response,roc_input_sub$Freq)
  tmp=tmp$auc
  names(tmp)=i
  auc=c(auc,tmp)
  
}

plot_dat=as.data.frame(auc)
plot_dat$cell=rownames(plot_dat)
plot_dat=plot_dat[order(plot_dat$auc),]
plot_dat$cell=factor(plot_dat$cell,levels = rownames(plot_dat))

ggplot(plot_dat, aes(x = auc, y = cell,fill = auc,color=auc)) +
  geom_col( width = 0.7) +
  scale_fill_gradient(low = "grey90", high = "grey9")+
  scale_color_gradient(low = "grey90", high = "grey9")+
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  labs(
    x = "AUC",
    y = NULL    ) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 8),  
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0("auc_before_sam_response.pdf"),width = 4.5,height = 6)

auc=c()
roc_input_post=roc_input[roc_input$treat%in%c("post"),]

for (i in unique(roc_input$cell)) {
  
  roc_input_sub=roc_input_post[roc_input_post$cell%in%i,]
  tmp=roc(roc_input_sub$response,roc_input_sub$Freq)
  tmp=tmp$auc
  names(tmp)=i
  auc=c(auc,tmp)
  
}

plot_dat=as.data.frame(auc)
plot_dat$cell=rownames(plot_dat)
plot_dat=plot_dat[order(plot_dat$auc),]
plot_dat$cell=factor(plot_dat$cell,levels = rownames(plot_dat))

ggplot(plot_dat, aes(x = auc, y = cell,fill = auc,color=auc)) +
  geom_col( width = 0.7) +
  scale_fill_gradient(low = "grey90", high = "grey9")+
  scale_color_gradient(low = "grey90", high = "grey9")+
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  labs(
    x = "AUC",
    y = NULL    ) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),    
 panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0("auc_post_sam_response.pdf"),width = 4.5,height = 6)

auc=c()

for (i in unique(roc_input$cell)) {
  
  roc_input_sub=roc_input[roc_input$cell%in%i,]
  tmp=roc(roc_input_sub$treat,roc_input_sub$Freq)
  tmp=tmp$auc
  names(tmp)=i
  auc=c(auc,tmp)
  
}

plot_dat=as.data.frame(auc)
plot_dat$cell=rownames(plot_dat)
plot_dat=plot_dat[order(plot_dat$auc),]
plot_dat$cell=factor(plot_dat$cell,levels = rownames(plot_dat))

ggplot(plot_dat, aes(x = auc, y = cell,fill = auc,color=auc)) +
  geom_col( width = 0.7) +
  scale_fill_gradient(low = "grey90", high = "grey9")+
  scale_color_gradient(low = "grey90", high = "grey9")+
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  labs(
    x = "AUC",
    y = NULL    ) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(), 
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0("auc_all_sam_treat.pdf"),width = 4.5,height = 6)

auc=c()
roc_input_A=roc_input[roc_input$response%in%c("A"),]

for (i in unique(roc_input$cell)) {
  
  roc_input_sub=roc_input_A[roc_input_A$cell%in%i,]
  tmp=roc(roc_input_sub$treat,roc_input_sub$Freq)
  tmp=tmp$auc
  names(tmp)=i
  auc=c(auc,tmp)
  
}

plot_dat=as.data.frame(auc)
plot_dat$cell=rownames(plot_dat)
plot_dat=plot_dat[order(plot_dat$auc),]
plot_dat$cell=factor(plot_dat$cell,levels = rownames(plot_dat))

ggplot(plot_dat, aes(x = auc, y = cell,fill = auc,color=auc)) +
  geom_col( width = 0.7) +
  scale_fill_gradient(low = "grey90", high = "grey9")+
  scale_color_gradient(low = "grey90", high = "grey9")+
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  labs(
    x = "AUC",
    y = NULL    ) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 8),  
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0("auc_A_sam_treat.pdf"),width = 4.5,height = 6)

auc=c()
roc_input_B=roc_input[roc_input$response%in%c("B"),]

for (i in unique(roc_input$cell)) {
  
  roc_input_sub=roc_input_B[roc_input_B$cell%in%i,]
  tmp=roc(roc_input_sub$treat,roc_input_sub$Freq)
  tmp=tmp$auc
  names(tmp)=i
  auc=c(auc,tmp)
  
}

plot_dat=as.data.frame(auc)
plot_dat$cell=rownames(plot_dat)
plot_dat=plot_dat[order(plot_dat$auc),]
plot_dat$cell=factor(plot_dat$cell,levels = rownames(plot_dat))

ggplot(plot_dat, aes(x = auc, y = cell,fill = auc,color=auc)) +
  geom_col( width = 0.7) +
  scale_fill_gradient(low = "grey90", high = "grey9")+
  scale_color_gradient(low = "grey90", high = "grey9")+
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  labs(
    x = "AUC",
    y = NULL    ) +
  theme_bw()+
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    panel.grid.major.y = element_blank(),  
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0("auc_B_sam_treat.pdf"),width = 4.5,height = 6)

library(ggplot2)
library(ggmapcn)
library(sf)
library(ggspatial) 

allsam=unique(ST$sample_id)


for (sam in allsam[1:length(allsam)]) {
  
sf=paste0(sam,".geojson")
map= read_sf(sf)
if (map$objectType[1]=="annotation") {
  map=map[2:nrow(map),]
}
tmp=ST@meta.data
fib.cell <- rownames(tmp[tmp$sample_id%in%sam,])
tmp=tmp[fib.cell,]
p1<- ggplot(st_set_crs(map, 3006), aes(fill =fincell,color=fincell)) + 
  geom_sf(size=0.1) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


ggsave(plot = p1,paste0("ann_",sam,"_type1.pdf"), width = 120, height = 120,limitsize = F)

p2<- ggplot(st_set_crs(map, 3006), aes(fill =fincell,color=fincell)) + 
  geom_sf(size=0.1) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black"))


ggsave(plot = p2,paste0("ann_",sam,"_type2.pdf"), width = 120, height = 120,limitsize = F)
p2<- ggplot(st_set_crs(map, 3006), aes(fill =fincell,color=fincell)) + 
  geom_sf(size=0.1) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "grey80"))
ggsave(plot = p2,paste0("ann_",sam,"_type3.pdf"), width = 120, height = 120,limitsize = F)


}
#fig5e,f,g

library(stringr)
library(tidyverse)
library(pheatmap)
m <- as.data.frame(gibbs.res$Beta)
pheatmap::pheatmap(t(m))
niche=as.vector(gibbs.res$cell_topics)
niche=paste0("CN",niche)
ST$niche=niche
data.combined=ST
cellratio <- prop.table(table(data.combined$niche, data.combined$sam_group), margin = 2)%>%
  as.data.frame()%>%
  `colnames<-`(c('celltype','orig.ident','Freq'))
group=strsplit(as.vector(cellratio$orig.ident),":")
group=lapply(group,FUN = function(x){x[2]})
group=unlist(group)
cellratio$group <-group
cellratio1 <- cellratio%>%
  group_by(celltype, group) %>%
  summarise(avg.freq = mean(Freq),
            sd.freq = sd(Freq),
            se.freq = sd(Freq)/sqrt(n()))
ggplot(data = cellratio1, mapping = aes(x = celltype, y = avg.freq,
                                        color =group, fill = group)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=sample_color)+
  scale_color_manual(values=sample_color)+
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(data = cellratio1,
                aes(ymin=avg.freq-se.freq, ymax=avg.freq+se.freq),
                width=0.3,
                size = 0.5,
                color = "black",
                position=position_dodge(0.9)) +
  geom_jitter(data = cellratio,
              mapping = aes(x = celltype, y = Freq),
              shape=21,alpha=0.8,
              color = "black",
              position = position_dodge(width = 0.9),
              size = 1) 
ggsave('box_all2.pdf', width = 13, height = 4.6)

ggplot(data = cellratio1, mapping = aes(x = celltype, y = avg.freq,
                                        color =group, fill = group)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=sample_color)+
  scale_color_manual(values=sample_color)+
  geom_bar(stat = "identity", position=position_dodge()) +
  geom_errorbar(data = cellratio1,
                aes(ymin=avg.freq-se.freq, ymax=avg.freq+se.freq),#加误差线
                width=0.3,
                size = 0.5,
                color = "black",
                position=position_dodge(0.9)) 
ggsave('box_all4.pdf', width = 10, height = 4.6)

library(ggpubr)

for (i in unique(cellratio$celltype)) {
  
  cellratio_sub=cellratio[cellratio$celltype==i,]
  cellratio1_sub <- cellratio_sub%>%
    group_by(celltype, group) %>%
    summarise(avg.freq = mean(Freq),
              sd.freq = sd(Freq),
              se.freq = sd(Freq)/sqrt(n()))
  ebtop<-function(x){
    return(mean(x)+sd(x)/sqrt(length(x)))
  }
  ebbottom<-function(x){
    return(mean(x)-sd(x)/sqrt(length(x)))
  }
  
  c <- ggplot(data=cellratio_sub,aes(x=group,y=Freq,fill=group))+
    stat_summary(geom = "bar",
                 fun = mean,
    )+
    stat_summary(geom = "errorbar",
                 fun.min = ebbottom,
                 fun.max = ebtop,
                 width=0.2)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10))  + 
    scale_fill_manual(values=sample_color)+
    scale_color_manual(values=sample_color)+
    geom_jitter(aes(fill=group),width = 0.2,size=1,shape=21,color="black") + 
    stat_compare_means(comparisons =my_comparisons, simplify =FALSE)
  ggsave(plot = c,paste0("box1_",i,".pdf"), width=3.8, height=4.3)
  p <- ggplot(cellratio_sub,
              aes(x=group, y=Freq, 
                  color=group
              )) + 
    geom_boxplot(position = position_dodge(0.8)) +
    geom_point(position = position_jitterdodge()) +
    scale_color_manual(values = sample_color) +
    theme_bw()+#+stat_compare_means() 
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10))  + 
    stat_compare_means(comparisons =my_comparisons, simplify =FALSE)
  ggsave(plot = p,paste0("box2_",i,".pdf"), width=3.8, height=4)
}

cellratio <- prop.table(table(data.combined$niche, data.combined$sam_group), margin = 2)%>%
  as.data.frame()%>%
  `colnames<-`(c('celltype','orig.ident','Freq'))
niche_res=gibbs.res
cell_niches=niche_res$cell_topics
cell_niches=as.vector(cell_niches)
ST$cell_niche=cell_niches
ST$niche=factor(ST$niche,levels = paste0("CN",sort(unique(cell_niches))))
qs::qsave(ST, "ST_with_niche.qs")

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)  
library(RColorBrewer)

df_bar=as.data.frame(prop.table(table(ST$fincell)))
colnames(df_bar)=c("celltype","ratio")
df_bar=df_bar[order(df_bar$ratio,decreasing = T),]
df_bar$celltype=factor(df_bar$celltype,levels =df_bar$celltype )

p_top<-ggplot(df_bar, aes(x = celltype, y = ratio, fill = celltype)) +
  geom_bar(stat = "identity", width = 1) +  
  scale_fill_manual(values = celltype_col) +
  #scale_y_continuous(limits = c(0, 0.15), breaks = c(0, 0.15), labels = c("0.00", "0.15")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "none",  
    panel.grid = element_blank()
  ) +
  labs(title = "Cell ratio") 
p_top

table=prop.table(table(ST$niche,ST$fincell),margin = 2)
pheatmap(table)

cluster <- pheatmap(table, cluster_cols=F, silent=TRUE)$tree_row

write.csv(table,"table.csv")
table <- read.csv("./table.csv", row.names=NULL)
ratio_table=as.data.frame(table)
ratio_table=reshape2::melt(ratio_table)
ratio_table=ratio_table[,c(1,2,4)]
colnames(ratio_table)=c("niche","celltype","ratio")



niche_order=rev(paste0("CN",cluster$order))

color_bins <- seq(0, 1, length.out=10)
color_palette <- rev(brewer.pal(9, "Spectral")) 
ratio_table$niche=factor(ratio_table$niche,levels = niche_order)
ratio_table$celltype=factor(ratio_table$celltype,levels =df_bar$celltype )

p_point <- ggplot(ratio_table[,], aes(x=celltype, y=niche, size=ratio, color=ratio)) +
  geom_tile(fill = "white", color = "black", size = 0.8) +
  geom_point() +
  scale_color_gradientn(colors=color_palette, limits=c(0, 1), breaks=color_bins) +
  scale_size_continuous(range=c(1, 10), limits=c(0, 1)) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=10),
    axis.text.y = element_text(size=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(x="Cell Type", y=element_blank(), size="Cell Ratio", color="Cell Ratio")
p_point

library(ggdendro)
dend_data <- dendro_data(cluster)

p_dend <- ggplot(segment(dend_data),rotate = T) +
  geom_segment(aes(x=x, y=-y, xend=xend, yend=-yend)) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank())

p_dend

df_bar2=as.data.frame(prop.table(table(ST$niche)))
colnames(df_bar2)=c("niche","ratio")
df_bar2$niche=factor(df_bar2$niche,levels = niche_order )

p_right<-ggplot(df_bar2, aes(x = niche, y = ratio, fill = ratio)) +
  geom_bar(stat = "identity", width = 1) +  
  scale_fill_gradientn(colors = color_palette, limits = c(0, max(df_bar2$ratio))) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",  
    panel.grid = element_blank()
  ) +
  labs(title = "niche ratio") + 
  coord_flip() 
p_right

library(patchwork)

p_top 
ggsave("p_niche_celltype_top.pdf",width = 8.5,height = 2.5)

p_dend + p_point + p_right+ plot_layout(widths=c(0.8, 5,1))
ggsave("p_niche_celltype_point.pdf",width = 12,height = 7)

df_bar=as.data.frame(prop.table(table(ST$niche)))
colnames(df_bar)=c("niche","ratio")
df_bar=df_bar[order(df_bar$ratio,decreasing = T),]
df_bar$niche=factor(df_bar$niche,levels =df_bar$niche )

p_top<-ggplot(df_bar, aes(x = niche, y = ratio, fill = ratio)) +
  geom_bar(stat = "identity", width = 1) +  
  scale_fill_gradientn(colors = color_palette, limits = c(0, max(df_bar2$ratio))) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_blank(),
    legend.position = "none",  
    panel.grid = element_blank()
  ) +
  labs(title = "niche ratio") 
p_top

df_bar2=as.data.frame(prop.table(table(ST$sample_id,ST$group)))

colnames(df_bar2)=c("sample_id","group","ratio")

df_bar2=df_bar2[df_bar2$ratio>0,]

df_bar2=df_bar2[order(df_bar2$ratio,decreasing = F),]
df_bar2=df_bar2[order(df_bar2$group,decreasing = T),]

df_bar2$sample_id=factor(df_bar2$sample_id,levels = rev(df_bar2$sample_id) )


p_right<-ggplot(df_bar2, aes(x = sample_id, y = ratio, fill = ratio)) +
  geom_bar(stat = "identity", width = 1) + 
  scale_fill_gradientn(colors = color_palette, limits = c(0, max(df_bar2$ratio))) +
   theme(
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",  
    panel.grid = element_blank()
  ) +
  labs(title = "cell ratio") +  
  coord_flip()  
p_right

table=prop.table(table(ST$sample_id,ST$niche),margin = 1)
pheatmap(table)
cluster <- pheatmap(table, cluster_cols=F, silent=TRUE)$tree_row
ratio_table=as.data.frame(table)
colnames(ratio_table)=c("sample_id","niche","ratio")
color_bins <- seq(0, 1, length.out=10)
color_palette <- rev(brewer.pal(9, "Spectral"))  
red_palette <- colorRampPalette(c("#F5CBA7", "#F1948A", "#E74C3C", "#C0392B"))(10)

ratio_table$sample_id=factor(ratio_table$sample_id,levels = rev(df_bar2$sample_id))
ratio_table$niche=factor(ratio_table$niche,levels =df_bar$niche )
p_point <- ggplot(ratio_table[,], aes(x=niche, y=sample_id, size=ratio, color=ratio)) +
  geom_tile(fill = "white", color = "black", size = 0.8) +
  geom_point() +
  scale_color_gradientn(
    colors = red_palette,
    limits = c(0, max(ratio_table$ratio)),
    breaks = c(0.0, 0.1, 0.2, 0.3, 0.4),
    name = "Ratio"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=10),
    axis.text.y = element_blank(),
    #legend.position="right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(x="Niche", y=element_blank(), size="Niche Ratio", color="Niche Ratio")

p_point

p_left <- ggplot(df_bar2, aes(x = 1, y = sample_id, fill = group)) +
  geom_tile(width = 0.8, height = 1) +
  scale_fill_manual(values =sample_color) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

library(patchwork)

p_top 
ggsave("p_niche_sample_top.pdf",width = 5.7,height = 1.7)

p_left + p_point + p_right+ plot_layout(widths=c(0.1, 5,1.5))
ggsave("p_niche_sample_point.pdf",width = 8,height = 9)
meta=ST@meta.data


for (i in as.vector(unique(meta$niche))) {
  meta_sub=meta[meta$niche==i,]
  df_bar_sub=as.data.frame(prop.table(table(meta_sub$fincell)))
  colnames(df_bar_sub)=c("celltype","ratio")
  df_bar_sub$celltype=factor(df_bar_sub$celltype,levels = df_bar_sub$celltype )

  
  
  p_right<- ggplot(df_bar_sub, aes(x = celltype, y = ratio, fill = celltype)) +
    geom_bar(stat = "identity", width = 1) +  
    scale_fill_manual(values = celltype_col) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  
      axis.title.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",  
      panel.grid = element_blank()
    ) +
    labs(title = "cell ratio") 
  df_point_sub=df_bar_sub
  df_point_sub$celltype=factor(df_point_sub$celltype,levels = df_bar_sub$celltype)
  red_palette <- colorRampPalette(c("#F5CBA7", "#F1948A", "#E74C3C", "#C0392B"))(10)
  p_point <- ggplot(df_point_sub[,], aes(x=celltype, y=1, size=ratio, color=ratio)) +
    geom_tile(fill = "white", color = "black", size = 0.8) +
    geom_point() +
    scale_color_gradientn(
      colors = red_palette,
      limits = c(0, max(df_point_sub$ratio)),
      breaks = c(0,max(df_point_sub$ratio)),
      name = "Ratio"
    ) +    theme_minimal() +
    theme(legend.position = "bottom")+
    theme(
      axis.text.x = element_text(angle=45, hjust=1, size=10),
      axis.text.y = element_blank(),
      #legend.position="right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    labs(x="celltype", y=element_blank(), size="celltype Ratio", color="celltype Ratio")
   p_right / p_point+ plot_layout(height=c(2, 0.8))
  ggsave(paste0("niche_sub_",i,".pdf"),width = 8,height = 4)
}
library(ggplot2)
library(ggmapcn)
library(sf)
library(ggspatial) 
allsam=unique(ST$sample_id)
for (sam in allsam[11:length(allsam)]) {
  
  sf=paste0(sam,".geojson")
  map= read_sf(sf)
  if (map$objectType[1]=="annotation") {
    map=map[2:nrow(map),]
  }
  tmp=ST@meta.data
  fib.cell <- rownames(tmp[tmp$sample_id%in%sam,])
  tmp=tmp[fib.cell,]
  map$niche<-as.vector(tmp$niche)
  p1<- ggplot(st_set_crs(map, 3006), aes(fill =niche,color=niche)) + 
    geom_sf(size=0.8) +
    scale_fill_manual(values = niche_col)+
    scale_color_manual(values = niche_col)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  ggsave(plot = p1,paste0("ann_",sam,"_type1.pdf"), width = 50, height = 50,limitsize = F)
  p2<- ggplot(st_set_crs(map, 3006), aes(fill =niche,color=niche)) + 
    geom_sf(size=0.8) +
    scale_fill_manual(values = niche_col)+
    scale_color_manual(values = niche_col)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black"))
  ggsave(plot = p2,paste0("ann_",sam,"_type2.pdf"), width = 50, height = 50,limitsize = F)
  p2<- ggplot(st_set_crs(map, 3006), aes(fill =niche,color=niche)) + 
    geom_sf(size=0.8) +
    scale_fill_manual(values = niche_col)+
    scale_color_manual(values = niche_col)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey80"))
  ggsave(plot = p2,paste0("ann_",sam,"_type3.pdf"), width = 50, height = 50,limitsize = F)
  
  
}
#fig5i
library(Seurat)
library(ggvoronoi)
library(ggplot2)
library(dplyr)
library(ggforce)
celltype_niche=c(celltype,"out_of_niche")
cols=c(celltype_col,"white")
names(cols)=celltype_niche
celltype_col_niche=cols
allsam=unique(ST$sample_id)
for (sam in allsam) {
  tmp=ST@meta.data
  data_fin_sub=tmp[tmp$sample_id==sam,]
  data_fin_sub$fincell=as.vector(data_fin_sub$fincell)
  data_fin_sub$niche=as.vector(data_fin_sub$niche)
  for (niche_use in unique(data_fin_sub$niche)) {
    data_fin_sub_sub=data_fin_sub
    data_fin_sub_sub$celltype="out_of_niche"
    data_fin_sub_sub[data_fin_sub_sub$niche==niche_use,]$celltype<-as.vector(data_fin_sub_sub[data_fin_sub_sub$niche==niche_use,]$fincell)
    
    data_fin3=data_fin_sub_sub

    map<-  ggplot(data_fin3, aes(x, y, group = -1L)) +
      geom_voronoi_tile(aes(fill = celltype),color="grey90", normalize = F,expand = unit(-1.5, 'pt'), radius = unit(0, 'pt'), max.radius = 10)+
      scale_fill_manual(values=celltype_col_niche)+
      theme()+
      theme(panel.background = element_rect(fill = 'white'),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) 
    
    
    factor=250/(max(data_fin3$y)-min(data_fin3$y))
    y=(max(data_fin3$x)-min(data_fin3$x))*factor
    ggsave(paste0("voronoi_",sam,"_",niche_use,".pdf"),plot =map, width = y, height = 250,limitsize = F) 
  }
}
