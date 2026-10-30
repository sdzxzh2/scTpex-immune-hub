
library(SingleCellExperiment)
library(reshape2)
library(tidyverse)
library(dplyr)
library(data.table) 
library(ggplot2)
library(ComplexHeatmap)
library(colorRamps)
library(circlize)
library(RColorBrewer)
library(ggbeeswarm)
library(destiny)
library(scater)
library(dittoSeq)
library(gridExtra)
library(ggpmisc)
library(cowplot)
library(viridis)
library(ggpubr)
library(rstatix)
library(sf)
library(concaveman)
library(RANN)
library(pheatmap)
library(SummarizedExperiment)
library(imcRtools)
library(Seurat)


library(ggplot2)
library(ggmapcn)
library(sf)
library(ggspatial) 


#fig6a,b
allsam=unique(ST$sample_id)
celltype=names(celltype_col)
score_lis=list()
for (sam in allsam) {
  
  excluded_score_sub=paste0("excluded_score_",sam,".csv")
  excluded_score_sub=read.csv(excluded_score_sub, row.names=NULL)
  rownames(excluded_score_sub)=colnames(excluded_score_sub)
  new_mat=excluded_score_sub[celltype,]
  rownames(new_mat)=celltype
  new_mat=as.data.frame(t(new_mat))
  new_mat=new_mat[celltype,]
  rownames(new_mat)=celltype
  new_mat=as.data.frame(t(new_mat))
  new_mat[is.na(new_mat)]<-0
  new_mat[new_mat==Inf]<-min(new_mat)
  score_lis[[sam]]=new_mat

}

datR <- a_before
for(i in 1:nrow(datR)){
  datR[i,1:i] <- b_before[i,1:i]
}
datR[1:3,1:3]

datR[datR>10]<-10
datR[datR<(-10)]<- -10
datR=datR/10
p1<-Heatmap(datR, rect_gp = gpar(type = "none"), 
        show_heatmap_legend = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height,
                    gp = gpar(col = "grey", fill = NA))
          if(i == j) {
            grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
          }else if(i > j) {
            grid.circle(x = x, y = y, r =  0.9*abs(datR[i, j])/2 * min(unit.c(width, height)), 
                        gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
          } else {
            grid.circle(x = x, y = y, r =  0.9*abs(datR[i, j])/2 * min(unit.c(width, height)), 
                        gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
          }
        },
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = T, show_column_names = T, 
        row_names_side = "right", 
        row_names_rot = 45,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8)
)

p1

lgdRight <- Legend(col_fun = colCorRight, title = "a_before", 
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft, title = "b_before", 
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft)


datR <- a_post
for(i in 1:nrow(datR)){
  datR[i,1:i] <- b_post[i,1:i]
}
datR[1:3,1:3]

datR[datR>10]<-10
datR[datR<(-10)]<- -10
datR=datR/10
p1<-Heatmap(datR, rect_gp = gpar(type = "none"), 
            show_heatmap_legend = F,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.rect(x = x, y = y, width = width, height = height,
                        gp = gpar(col = "grey", fill = NA))
              if(i == j) {
                grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
              }else if(i > j) {
                grid.circle(x = x, y = y, r =  0.9*abs(datR[i, j])/2 * min(unit.c(width, height)), 
                            gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
              } else {
                grid.circle(x = x, y = y, r =  0.9*abs(datR[i, j])/2 * min(unit.c(width, height)), 
                            gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
              }
            },
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = T, show_column_names = T, 
            row_names_side = "right", 
            row_names_rot = 45,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8)
)

p1

lgdRight <- Legend(col_fun = colCorRight, title = "a_post", 
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft, title = "b_post", 
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft)

pdf("neighbor2.pdf", width = 5, height = 5.8)
draw(p1, annotation_legend_list = pd,
     annotation_legend_side = "top")
dev.off()
b_vs_a_before=b_before-a_before
b_vs_a_post=b_post-a_post
datR <- b_vs_a_before
for(i in 1:nrow(datR)){
  datR[i,1:i] <- b_vs_a_post[i,1:i]
}
datR[1:3,1:3]

datR[datR>10]<-10
datR[datR<(-10)]<- -10
datR=datR/10
p1<-Heatmap(datR, rect_gp = gpar(type = "none"), 
            show_heatmap_legend = F,
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.rect(x = x, y = y, width = width, height = height,
                        gp = gpar(col = "grey", fill = NA))
              if(i == j) {
                grid.circle(x = x, y = y, r = 0.5 * min(unit.c(width, height)), gp = gpar(fill = "grey", col = NA))
              }else if(i > j) {
                grid.circle(x = x, y = y, r =  0.9*abs(datR[i, j])/2 * min(unit.c(width, height)), 
                            gp = gpar(fill = colCorLeft(datR[i, j]), col = NA))
              } else {
                grid.circle(x = x, y = y, r =  0.9*abs(datR[i, j])/2 * min(unit.c(width, height)), 
                            gp = gpar(fill = colCorRight(datR[i, j]), col = NA))
              }
            },
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = T, show_column_names = T, 
            row_names_side = "right", 
            row_names_rot = 45,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8)
)

p1

lgdRight <- Legend(col_fun = colCorRight, title = "b_vs_a_before", 
                   direction = "horizontal")
lgdLeft <- Legend(col_fun = colCorLeft, title = "b_vs_a_post", 
                  direction = "horizontal")
pd = list(lgdRight, lgdLeft)

pdf("neighbor3.pdf", width = 5, height = 5.8)
draw(p1, annotation_legend_list = pd,
     annotation_legend_side = "top")
dev.off()
ST=qs::qread("ST_ann_fin3_with_main_type.qs")
library(reshape2)
allsam=unique(ST$sample_id)
long_dat=data.frame()
for (i in 1:length(score_lis)) {
  tmp_score=score_lis[[i]]
  tmp_score$rowname=rownames(tmp_score)
  tmp_score=melt(tmp_score)
  tmp_score$sam=allsam[i]
  long_dat=rbind(long_dat,tmp_score)
  
}

group=rep(unique(ST$group),c(441*5,441*5,441*5,441*5))
long_dat$group=group
p_A_before_B_before=matrix(nrow = length(celltype), ncol =  length(celltype))
rownames(p_A_before_B_before)=celltype
colnames(p_A_before_B_before)=celltype
p_A_post_B_post=p_A_before_B_before
p_A_before_A_post=p_A_before_B_before
p_B_before_B_post=p_A_before_B_before
for (i in celltype) {
  long_dat_i=long_dat[long_dat$rowname==i,]
  for (j in celltype) {   
    long_dat_ij=long_dat_i[long_dat_i$variable==j,]
    p1=t.test(long_dat_ij[long_dat_ij$group=="A_before",]$value,long_dat_ij[long_dat_ij$group=="B_before",]$value)
    p1=p1$p.value
    p2=t.test(long_dat_ij[long_dat_ij$group=="A_post",]$value,long_dat_ij[long_dat_ij$group=="B_post",]$value)
    p2=p2$p.value
    p3=t.test(long_dat_ij[long_dat_ij$group=="A_before",]$value,long_dat_ij[long_dat_ij$group=="A_post",]$value)
    p3=p3$p.value
    p4=t.test(long_dat_ij[long_dat_ij$group=="B_before",]$value,long_dat_ij[long_dat_ij$group=="B_post",]$value)
    p4=p4$p.value
    if (i==j) {
      p1=p2=p3=p4=1
    }
    p_A_before_B_before[i,j]=p1
    p_A_post_B_post[i,j]=p2
    p_A_before_A_post[i,j]=p3
    p_B_before_B_post[i,j]=p4
    
  }
  
  
}

ratio_table=b_post - b_before
p_use=p_B_before_B_post
p_use=as.data.frame(p_use)
p_use[is.na(p_use)]<-1
ratio_table$rowname=rownames(ratio_table)
ratio_table=melt(ratio_table)
colnames(ratio_table)=c("x","y","val")
ratio_table$x=factor(ratio_table$x,levels = celltype)
ratio_table$y=factor(ratio_table$y,levels =celltype )
p_use=-log10(p_use)
p_use$rowname=rownames(p_use)
p_use=melt(p_use)
colnames(p_use)=c("x","y","p")
p_use$x=factor(p_use$x,levels = celltype)
p_use$y=factor(p_use$y,levels =celltype )
red_palette <- colorRampPalette(c("#F5CBA7", "#F1948A", "#E74C3C", "#C0392B"))(10)
red_palette2 <- colorRampPalette(c("#3288BD", "white", "#C0392B"))(10)
ratio_table[ratio_table$val>10,]$val<-10

ratio_table[ratio_table$val<(-10),]$val<- (-10)

ggplot() +
  geom_tile(data=ratio_table[,], aes(x=x, y=y,fill=val,color =val), size = 0.8) +
  scale_fill_gradientn(
    colors = red_palette2,
    limits = c(min(ratio_table$val), max(ratio_table$val)),
    name = "Ratio"
  ) +
  scale_color_gradientn(
    colors = red_palette2,
    limits = c(min(ratio_table$val), max(ratio_table$val)),
    name = "Ratio"
  ) +
  geom_point(data =p_use,aes(x=x, y=y, size=p, ),color="#EE9B48" ) +
  theme_minimal() +
  theme(legend.position = "bottom")+
  theme(
    axis.text.x = element_text(angle=90, hjust=1, size=10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(x=" ", y=element_blank(), size=" ", color=" ")

#fig6c,d,e
ST=qs::qread("ST_ann_fin3_with_main_type.qs")
seu<-ST
seu[["RNA"]]<-JoinLayers(seu[['Spatial']])
spe=as.SingleCellExperiment(seu,assay ="RNA")
spe$cellID=spe$name
spe$Center_X=spe$x
spe$Center_Y=spe$y
spe$ImageNumber=spe$sample_id
meta=ST@meta.data
meta$cell_patch_raw=cell_patch
meta$cell_community_50_raw=cell_community_50
meta$cell_community_0_raw=cell_community_0
cell_patch=strsplit(as.vector(cell_patch),":")
cell_patch=lapply(cell_patch,FUN = function(x){x[2]})
cell_patch=unlist(cell_patch)
cell_community_50=strsplit(as.vector(cell_community_50),":")
cell_community_50=lapply(cell_community_50,FUN = function(x){x[2]})
cell_community_50=unlist(cell_community_50)
cell_community_0=strsplit(as.vector(cell_community_0),":")
cell_community_0=lapply(cell_community_0,FUN = function(x){x[2]})
cell_community_0=unlist(cell_community_0)  
meta$cell_patch=cell_patch
meta$cell_community_50=cell_community_50
meta$cell_community_0=cell_community_0
meta$Milieu_type="outside"
meta[meta$cell_community_50!="community_0",]$Milieu_type="boundary"
meta[meta$cell_community_0!="inside_0",]$Milieu_type="intra"
meta$cell_milieu=paste0(meta$fincell,"_",meta$Milieu_type)
cellratio=roc_input
cellratio$celltype=cellratio$cell
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

#fig6f

library(ggplot2)
library(ggmapcn)
library(sf)
library(ggspatial) 
allsam=unique(ST_ponit$sample_id)
for (sam in allsam[1:length(allsam)]) {
  sf=paste0(sam,".geojson")
  map= read_sf(sf)
  if (map$objectType[1]=="annotation") {
    map=map[2:nrow(map),]
  }
  tmp=ST_ponit@meta.data
  fib.cell <- rownames(tmp[tmp$sample_id%in%sam,])
  tmp=tmp[fib.cell,]
  p1<- ggplot(st_set_crs(map, 3006), aes(fill =niche,color=niche)) + 
    geom_sf(size=0.8) +
    scale_fill_manual(values = celltype_col_ponit)+
    scale_color_manual(values = celltype_col_ponit)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  ggsave(plot = p1,paste0("ann_",sam,"_type1.pdf"), width = 50, height = 50,limitsize = F)
  p2<- ggplot(st_set_crs(map, 3006), aes(fill =niche,color=niche)) + 
    geom_sf(size=0.8) +
    scale_fill_manual(values = celltype_col_ponit)+
    scale_color_manual(values = celltype_col_ponit)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black"))
  ggsave(plot = p2,paste0("ann_",sam,"_type2.pdf"), width = 50, height = 50,limitsize = F)
  p2<- ggplot(st_set_crs(map, 3006), aes(fill =niche,color=niche)) + 
    geom_sf(size=0.8) +
    scale_fill_manual(values = celltype_col_ponit)+
    scale_color_manual(values = celltype_col_ponit)+
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "grey80"))
  ggsave(plot = p2,paste0("ann_",sam,"_type3.pdf"), width = 50, height = 50,limitsize = F)
}


fig6g,h,i

library(tidyr)
library(data.table)
library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
require(ggpubr)
library(Seurat)
library(NMF)
library(ggradar)
library(corrr)
library(patchwork)

ranks <- 3:8
nrun <- 20  
meta_for_nmf$Patient_time_origin=paste0(meta_for_nmf$group,"-",meta_for_nmf$cell_community_50_raw)

sub_freq_df=data.frame()

for (i in unique(meta_for_nmf$Patient_time_origin)) {
  meta_sub=meta_for_nmf[meta_for_nmf$Patient_time_origin==i,]
  tmp=prop.table(table(meta_sub$Patient_time_origin,meta_sub$fincell),margin = 1)
  sub_freq_df=rbind(sub_freq_df,tmp)
}

colnames(sub_freq_df)<-c("Patient_time_origin","fincell","Freq")

nmf_input_long = sub_freq_df %>% 
  mutate(pt_stg_2 = Patient_time_origin) %>% separate(pt_stg_2, into = c('Patient','community'), sep = ':')# %>% 
nmf_input=reshape2::dcast(nmf_input_long,Patient_time_origin~fincell,value.var = "Freq")
rownames(nmf_input)=nmf_input[,1]
nmf_input_fin=nmf_input[2:ncol(nmf_input)]
nmf_input_fin[is.na(nmf_input_fin)]<-0
nmf_input_fin<-nmf_input_fin[,colSums(nmf_input_fin)!=0]
write.csv(nmf_input_fin,"nmf_input_fin.csv")
cophenetic_coeffs <- numeric(length(ranks))

for (i in seq_along(ranks)) {
  rank <- ranks[i]
  nmf_results <- nmf(nmf_input_fin, rank = rank, nrun = nrun, .options = "v", seed = 123456)
  cophenetic_coeffs[i] <- cophcor(nmf_results)
}

cophenetic_df <- data.frame(FactorizationRank = ranks, CopheneticCoefficient = cophenetic_coeffs)

ggplot(cophenetic_df, aes(x = FactorizationRank, y = CopheneticCoefficient)) +
  geom_point(color = "purple", size = 3) +
  geom_line(color = "purple") +
  labs(title = "Cophenetic correlation survey", x = "Factorization rank", y = "Coefficient") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 16))
nrank=4
mtx = nmf_input_fin[, 1:ncol(nmf_input_fin)]
nmf.res = nmf(t(nmf_input_fin[, 1:ncol(nmf_input_fin)]), rank=nrank, nrun=20, seed=123456)
w = nmf.res@fit@W
rownames(w) = colnames(mtx)
colnames(w) = paste0('NMF', seq(1,nrank))
w_ = w %>% t %>% scale
p=ComplexHeatmap::Heatmap(w_, width = 11, height = 4, column_km = nrank, 
                          cluster_rows = F, clustering_method_columns = 'single') 
heatmap = draw(p)
h = nmf.res@fit@H
rownames(h) = paste0('NMF', seq(1,nrank))
h_ = h %>% t %>% scale
p=ComplexHeatmap::Heatmap(h_) 
heatmap = draw(p)
library(patchwork)
library(ape)
library(plyr)
library(ggtree)
library(vegan)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(dplyr)
nrank=4
df_dist <- vegdist(as.data.frame(t(w_)),method = "euclidean")
df_hc <- hclust(df_dist,method="complete")
plot(as.dendrogram(df_hc),type="rectangle",horiz=F)
df_tree <- as.phylo(df_hc)
cut=cutree(df_hc,k=nrank)
cut_order=cut[rev(df_hc$order)]
cut=paste0("C",cut_order)
cut=factor(cut,levels = unique(cut))
cut_tab=table(cut)
cluster_num <- cut_tab
gro_num <- cumsum(cluster_num)
gro <- list(group1=names(cut_order[1:gro_num[1]]),
            group2=names(cut_order[(gro_num[1]+1):gro_num[2]]),
            group3=names(cut_order[(gro_num[2]+1):gro_num[3]]),
           group4=names(cut_order[(gro_num[3]+1):gro_num[4]])
            #group5=names(cut_order[(gro_num[4]+1):gro_num[5]])
)

tmp=lapply(gro, FUN = function(x){length(x)})
tmp=unlist(tmp)
tmp=names(sort(tmp))
gro=gro[tmp]
names(gro)=paste0("group",1:nrank)
tree<-groupOTU(df_tree,gro)
p2 <- ggtree(tree,size=1.2,aes(color=group),branch.length = "none")+
  theme(legend.position = 'none')+
  layout_dendrogram()+
  theme(plot.margin = margin(t = 0,  
                             r = 0,  
                             b = 0,  
                             l = 0,  
                             unit = "cm"))+
  scale_color_manual(values = c("#FF935F",rev(c("#A48AD3","#71BC7A","#1CC5FE","#FB7D80"))))+
  geom_tippoint(aes(color=group),
                show.legend = T,
                size=4)
p2
w_df=melt(as.matrix(w_))
colnames(w_df)<-c("NMF","celltype","value")
sam_order=unlist(gro)
names(sam_order)=NULL
w_df$celltype=factor(w_df$celltype,levels =sam_order)
w_df$NMF=factor(w_df$NMF,levels =c("NMF2","NMF1","NMF4","NMF3"))
p1=ggplot(w_df,aes(celltype,NMF,fill=value))+
  geom_tile(color="white",fill="white",size=0.7)+
  geom_point(pch=22,color="gray2",size=8.7)+
  theme(panel.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90,color = "black",size=13)
        
        # axis.text.y = element_blank()
  )+
  theme(plot.margin = margin(t = 0,  
                             r = 0,  
                             b = 0, 
                             l = 0,  
                             unit = "cm"))+
  scale_fill_gradientn(colours =colorRampPalette(c('white','white','gray88','#BB4444','#BB4444','#BB4444'))(111))+
  labs(x = NULL,y = NULL,fill="")
p2/p1+plot_layout(heights = c(1.3 ,1))
p2/p1+plot_layout(heights = c(1 ,0.8))box_input=as.data.frame(t(h))
box_input$id=rownames(box_input)
box_input = box_input %>% 
  separate(id, into = c('group','Patient_patch'), sep = '-') 
box_input=box_input
box_input$group=factor(box_input$group,levels = sort(unique(box_input$group)))
box_input$id=box_input$Patient
write.csv(box_input,"box_input.csv")
box_input_long=melt(box_input)
c <- ggplot(box_input_long,
            aes(x=group, y=value, 
                
                color = group
            )) + 
  geom_boxplot(notch = F, alpha = 0.95, 
               outlier.shape = 16,
               outlier.colour = "black", #outlier点用黑色
               outlier.size = 0.65) +
  facet_wrap(~ variable, scales = 'free_x', nrow = 2) +
  geom_point()+
  scale_color_manual(values = sample_color) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
        axis.text.y = element_text(angle = 90, size = 12),
        axis.title.y = element_text(angle = 90, size = 15)) +
  theme(legend.position = "top")
p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
ggsave('box1.pdf',plot = p, width = 6, height = 8)
library(scales)
library(ggradar)
library(tibble)
group="B"
radar_input <- apply(as.matrix(t(h)), 1, function(x) (x-min(x))/(max(x)-min(x))) %>% 
  data.frame(check.names = F)
radar_input=as.data.frame(t(radar_input))
radar_input_fin = radar_input %>% 
  rownames_to_column(var="id")%>% 
  separate(id, into = c('time_origin','Patient'), sep = '-')# %>% 
radar_input_fin$group2=radar_input_fin$Patient
radar_input_fin = radar_input_fin %>% 
  separate(group2, into = c('Patient','patch'), sep = ':')# %>% 
radar_input$group=radar_input_fin$time_origin
radar_input_sub=radar_input[radar_input$group==group,]
radar_input_sub=radar_input_sub[,1:(nrank)]
radar_input_sub <- radar_input_sub %>% 
  as_tibble(rownames = "group")

fig6j

library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)

LR_sec=tile
receptor <- unique(LR_sec$receptor)
receptor
LR_sec$Target <- LR_sec$receptor
for (i in 1:length(receptor)){
  LR_sec[, "receptor"][LR_sec[, "receptor"]==receptor[i]] = i
  
}
ligand <- unique(LR_sec$ligand)
LR_sec$ligand <- factor(LR_sec$ligand, levels = c(1:length(unique(LR_sec$ligand))))
LR_sec$receptor <- factor(LR_sec$receptor, levels = c(1:length(unique(LR_sec$receptor))))

if (length(unique(LR_sec$source))-length(unique(LR_sec$Target))>0) {
  p11 <- LR_sec %>% 
    ggparcoord(
      columns = 1:2,
      groupColumn = 3, 
      alphaLines = 0.6,scale="globalminmax",
      splineFactor=T)+
    scale_y_continuous(breaks = c(1:length(unique(LR_sec$source))),
                       labels =unique(LR_sec$source),
                       sec.axis = sec_axis(~./1, 
                                           breaks=c(1:length(unique(LR_sec$source))),
                                           labels=ylab_right))+
    theme(legend.position = 'none')+
    theme(plot.margin = margin(t = 0,  
                               r = 0,  
                               b = 0,  
                               l = 0,  
                               unit = "cm"))+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(colour =col))+
    scale_color_manual(values = col_unique)  
  
}else{
  p11 <- LR_sec %>% 
    ggparcoord(
      columns = 1:2, 
      groupColumn = 3, 
      showPoints = F, 
      alphaLines = 0.6,scale="globalminmax",
      splineFactor=T)+
    scale_y_continuous(breaks = c(1:length(unique(LR_sec$Target))),
                       labels =ylab_right,
                       sec.axis = sec_axis(~./1, 
                                           breaks=c(1:length(unique(LR_sec$Target))),
                                           labels=unique(LR_sec$Target)))+
    theme(legend.position = 'none')+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.background = element_blank(),
          legend.position = 'none',
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y.left = element_text(colour =col))+
    scale_color_manual(values =col_unique)#修改下连线的分组颜色，dittoColors()是dittoSeq包的一个颜色函数
  
}

p11

p1 <- LR_sec %>% 
  ggparcoord(
    columns = 1:2, 
    groupColumn = 3, 
    showPoints = F,
    alphaLines = 0.6,
    scale="globalminmax"
    )+
  annotate('text', x = 2, 
           y = c(1:length(unique(LR_sec$Target))) , 
           label = unique(LR_sec$Target),size=2 ,family="B",hjust=0)+
  annotate('text', x = 1, 
           y = c(1:length(unique(LR_sec$source))) , 
           label = unique(LR_sec$source),size=2 ,family="B",hjust=1,col=col)+
  theme(legend.position = 'none')+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        legend.position = 'none',
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())+
  scale_color_manual(values = col_unique )
library(Seurat)
data.combined=data.combined
Idents(data.combined)<-"group"
HC <-data.combined
ligand_dot <- DotPlot(HC, features = ligand,assay = "RNA")
ligand_exp <- ligand_dot$data
p2 <- ggplot(ligand_exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x=element_text(angle=90,hjust = 0,vjust=1,colour = 'black',size = 9),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 5))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#1A5592','white',"#B83D3D"),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(title = "",order = 1,nrow = 2))+
  theme(plot.margin=unit(c(1, 1, 1, 1),'cm'))
data.combined=data.combined_fib
Idents(data.combined)<-"group"
immune=data.combined

Target=LR_sec$Target
Target=strsplit(Target,"_")
target_fin=c()
for (i in 1:length(Target)) {
  tmp=Target[[i]]
  target_fin=c(target_fin,tmp[1])
}
p5 <- ggplot(receptor_exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "left")+
  theme(axis.text.x=element_text(angle=90,hjust = 0,vjust=1,colour = 'black',size = 9),
        axis.text.y = element_blank(),
        title = element_text(colour = 'black',size = 8),
        plot.title = element_text(hjust=0.5),
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 6))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#1A5592','white',"#B83D3D"),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  ggtitle("")+
  theme(plot.margin=unit(c(1, 1, 1, 1),'cm'))

p7 <- p5
p7

