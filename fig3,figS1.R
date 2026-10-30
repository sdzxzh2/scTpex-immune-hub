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


#fig3a
fib.cell=rownames(meta_all)
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined=sub_CT
data.combined$fincell=Idents(data.combined)
meta_all=meta_all[rownames(data.combined@meta.data),]
meta_all2=meta_all[order(meta_all$fincell),]
meta_all2=meta_all2[order(meta_all2$celltype_main),]
data.combined@meta.data=meta_all
data.combined$fincell=factor(data.combined$fincell,levels = unique(meta_all2$fincell))
Idents(data.combined)<-"fincell"

col=c("#FFE96C","#C6E0C0","#BE7EB2","#FFD0DF","#ADD167","#FFB95C","#FE846A","#BAB8D5","#7DC6C5","#FFFBAE","#ECA2B3",
      "#AFAFAF","#EBC596","#FFD306","#E988B7","#8B9DC7","#FF935F","#4BB39F",
      "#393874","#4F509D","#E97270","#DE4C00","#D0D597","#88672F","#F6B948","#F1CA94",
      "#803D34","#D96669","#ED969B","#A7508B","#CE67AE","#E29DCE","#A188BB","#D3CCE0","#E50065",
      "#BEADD3","#FFC489","#FFFA94","#B9B7D9","#6B68CB","#DE4C00","#0070BB","#71BC7A","#BC5300",
      "#962120","#BBD4E0","#39398A","#427256","#B96C27","#E97251","#6494B3",
      "#D1AA63","#989898","#59518A","#BBD4E0","#462672","#606F33","#511818","#050103","#B96C27","#6494B3",
      "#AF1E1F","#D9A833","#5F0C0C","#39398A","#60AFB3","#73398D","#F476BE","#A6D176","#B29BC9","#2D78AB",
      "#D5A9CD","#AAC5CF","#A75A36","#7071A6","#C5539F","#D2E0AC","#F59294","#97C3D7","#F2D3CA","#8CB4D2",
      "#957E9F","#608CC8","#F0F09B","#469D46","#D5D5D5","#989898","#674698","#F4B5A7","#CA9EC1","#E45777",
      "#F6B373","#F47C30","#E5838F","#916AA5"
      
)

pdf("Overview_Cluster.pdf",11,4.6)
DimPlot(data.combined, reduction = 'umap', label=F,raster=FALSE,cols =col ) 
dev.off()



#fig3b
#data.combined is the seuratObj is this study
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
## NMF program
# Constructing cellular abundance matrix

ranks <- 3:8  
nrun <- 20   
meta = data.combined@meta.data
meta$Patient_origin=paste0(meta$stim,"-",meta$group)
sub_freq_df=data.frame()

for (i in unique(meta$Patient_origin)) {
  meta_sub=meta[meta$Patient_origin==i,]
  tmp=prop.table(table(meta_sub$Patient_origin,meta_sub$fincell),margin = 1)
  sub_freq_df=rbind(sub_freq_df,tmp)
}

colnames(sub_freq_df)<-c("Patient_origin","fincell","Freq")

nmf_input_long = sub_freq_df %>% 
  mutate(pt_stg_2 = Patient_time_origin) %>% separate(pt_stg_2, into = c('Patient','time_origin'), sep = '-')  
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

nrank=5
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

#fig3d
df_igraph <- graph_from_data_frame(edge_list,node, directed = FALSE)
ggraph(df_igraph, layout = "linear",circular = TRUE)+
  geom_edge_fan(aes(color=color,edge_width=weight,edge_alpha = weight),show.legend = T)+
  geom_node_point(aes(fill=size,color=size),size=20,shape=21)+
  geom_node_text(aes(label=name, size = weight),
                  angle=0,hjust=0.1,size=3) +
  scale_edge_width_continuous(range = c(1,3),
                              limits = c(0.01,1),
                   #breaks = c(0.1),
                   guide = guide_legend(title = "Jaccard index",
                                        order=1))+
  scale_alpha_continuous(range = c(0.4,0.7))+
  scale_color_gradientn( colours =colorRampPalette(c("#FFED99","#FF7600"))(10))+
  scale_fill_gradientn(  colours =colorRampPalette(c("#FFED99","#FF7600"))(10),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"))+
   scale_edge_colour_manual(values="grey")+
  theme_graph()+
  theme(legend.position = 'none')+
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))

ggsave('net1_group4.pdf', width = 6, height = 6)

#fig3c
ggraph(df_igraph, layout = "linear",circular = TRUE)+
  geom_edge_fan(aes(color=color,edge_width=weight,edge_alpha = weight),show.legend = T)+
  geom_node_point(aes(fill=name,color=name),size=20,shape=21,alpha=1)+
  geom_node_text(aes(label=name, size = weight),
                 angle=0,hjust=0.1,size=3) +
  scale_edge_width_continuous(range = c(1,3),
                              limits = c(0.01,1),
                              #breaks = c(0.1),
                              guide = guide_legend(title = "Jaccard index",
                                                   order=1))+
  scale_alpha_continuous(range = c(0.4,0.7))+
  scale_color_manual(values = cols_for_net)+
  scale_fill_manual(values = cols_for_net)+
    scale_edge_colour_manual(values="grey")+
  theme_graph()+
  theme(legend.position = 'none')+
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))
ggsave('net2_group4.pdf', width = 6, height =6)

box_input=as.data.frame(t(h))
box_input$id=rownames(box_input)
box_input = box_input %>% 
  separate(id, into = c('Patient','group'), sep = '-') 
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
               outlier.colour = "black", 
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


#fig3e
celltype=c("Malignant_AMTN","Malignant_CSMD1","Malignant_FBN2",
           "Malignant_CFL1","Malignant_FTL","Malignant_S100A8","Malignant_KYNU","Malignant_NOTCH2","Malignant_FABP5",
           "Malignant_EMT","Malignant_EMT_MKI67","Malignant_EMT_TOP2A",
           "Malignant_AGR2","Malignant_SAA1","Malignant_ERBB4","Malignant_XIST"
           
)

col=c("#BBD4E0","#8CB4D2","#AAC5CF",
      "#FFB95C","#FFC489","#FF935F","#E97270","#D62C2D","#AF1E1F",
      "#73398D","#B29BC9","#A7508B",#"#A6D176","#A7508B",
      "#79992F","#71BC7A","#A6D176","#D0D597"
)

celltype_col=col
names(col)= celltype
names(celltype_col)= celltype
pdf("plot.pdf",5,5)
DimPlot(data.combined, reduction = 'umap', label=T,raster=FALSE,cols =celltype_col )
dev.off()


#fig3g,i
library(Seurat)
library(ggplot2)
library(UCell)
library(patchwork)
library(tidyr)
library(dplyr)
library(RColorBrewer)
library(GeneNMF) # 用于NMF
library(Seurat)
library(msigdbr)
library(fgsea)
library(viridis)

data.combined_fib$all="all"
seu.list <- SplitObject(data.combined_fib, split.by = "orig.ident")
geneNMF.programs <- multiNMF(seu.list, 
                             assay="RNA", slot="data", 
                             k=4:8, L1=c(0,0),
                             do_centering=TRUE, 
                             nfeatures = 2000)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=10,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.2)
geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        nprograms=9,
                                        max.genes=50,
                                        hclust.method="ward.D2",
                                        min.confidence=0.1)
genes=geneNMF.metaprograms$metaprograms.genes
writeGmtPathways(genes,"genes.gmt")
ph <- plotMetaPrograms(geneNMF.metaprograms, jaccard.cutoff = c(0.1,0.8),
                       palette = viridis(100, option = "A", direction = -1),
                       show_rownames = F)
geneNMF.metaprograms$metaprograms.metrics
t(as.data.frame(lapply(geneNMF.metaprograms$metaprograms.genes, head)))
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(data.combined_fib), category = "C5",subcategory = "GO:BP", pval.thr = 0.05)
 })
head(top_p)
top_p_path=lapply(top_p,FUN = function(x){x[[1]]})
writeGmtPathways(top_p_path,"top_go.gmt")
mp.genes <- geneNMF.metaprograms$metaprograms.genes
data.combined_fib <- AddModuleScore_UCell(data.combined_fib, features = mp.genes, ncores=14, name = "")
VlnPlot(data.combined_fib, features=names(mp.genes), group.by = "group",
        pt.size = 0, ncol=5)
library(viridis)
FeaturePlot(data.combined_fib, features = names(mp.genes), ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())

library(Seurat)
library(scRepertoire)
library(plyr)
library(dplyr)
library(ggplot2)
library(gProfileR)
library(plyr)
library(getopt)    # in o
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyverse)
celltype_sub=celltype
sigScores <- data.combined_fib@meta.data[,names(mp.genes)]
data=sigScores
id=data.combined_fib@meta.data$name
type=factor(data.combined_fib@meta.data$fincell,levels=celltype_sub)
group=as.vector(data.combined_fib@meta.data$group)
gene_use=colnames(data)
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
    theme(panel.grid = element_blank(), 
    geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means(label = "p.signif")
  ggsave(paste0("cell&group2_",i,".pdf"),plot = c, width=20, height=5)
  c <- ggplot(pdata_melt,
              aes(x=group, y=Signature_score, 
                  fill = group )) + 
    geom_boxplot(notch = F, alpha = 0.95,# color = "white",
                 outlier.shape = 16,
                 outlier.colour = "black",
                 outlier.size = 0.65) +
    scale_fill_manual(values = sample_color) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
          axis.text.y = element_text(angle = 90, size = 12),
          axis.title.y = element_text(angle = 90, size = 15)) +
    theme(legend.position = "top")+ stat_compare_means(label = "p.signif",comparisons = my_comparisons)
  ggsave(paste0("group1_",i,".pdf"),plot = c, width=2.6, height=5)
  
  c <- ggplot(pdata_melt, aes(x = group, y = Signature_score, fill = group)) +
    geom_violin( aes(color=group))+
    scale_fill_manual(values = sample_color)+
    scale_color_manual(values = sample_color)+
    theme_classic() +
    theme(panel.grid = element_blank() 
    ) +
    geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means()
  p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
  ggsave(paste0("group2_",i,".pdf"),plot = p, width=3.9, height=5)
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
    theme(legend.position = "top")#+stat_compare_means()
  
  ggsave(paste0("cell1_",i,".pdf"),plot = c, width=8, height=5)
  
  c <- ggplot(pdata_melt, aes(x = type, y = Signature_score, fill = type)) +
    geom_violin( aes(color=type))+
    scale_fill_manual(values = celltype_col)+
    scale_color_manual(values = celltype_col)+
    theme_classic() +
    theme(panel.grid = element_blank())
          ggsave(paste0("cell2_",i,".pdf"),plot = c, width=9, height=4)  
  for (t in celltype) { 
    tmp <- pdata_melt[pdata_melt$type==t,]
    c <- ggplot(tmp,
                aes(x=group, y=Signature_score, 
                    fill = group,   
                )) + 
      geom_boxplot(notch = F, alpha = 0.95, 
                   outlier.shape = 16,
                   outlier.colour = "black", 
                   outlier.size = 0.65) +
      scale_fill_manual(values = sample_color) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, hjust = 1,size = 10), 
            axis.text.y = element_text(angle = 90, size = 12),
            axis.title.y = element_text(angle = 90, size = 15)) +
      theme(legend.position = "top")
    p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
    ggsave(paste0("split_cell1_",i,"_",t,".pdf"),plot = p, width=2.6, height=5)
        c <- ggplot(tmp, aes(x = group, y = Signature_score, fill = group)) +
      geom_violin( aes(color=group))+
      scale_fill_manual(values = sample_color)+
      scale_color_manual(values = sample_color)+
      theme_classic() +
      theme(panel.grid = element_blank()
      geom_boxplot( aes(color=group),fill="white",width=.8,outlier.size = 0,outlier.stroke = 0,alpha = 0.4)+stat_compare_means()
    p <- c + stat_compare_means(label = "p.signif",comparisons = my_comparisons)
    ggsave(paste0("split_cell2_",i,"_",t,".pdf"),plot = p, width=3.5, height=5)
    
  }
 
}


#fig3f,h
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(vistributions)
library(ggplot2)
library(ggpubr)
library(ggsci)

dir_main_all="D:/cnv/"
for (group_type in unique(data.combined$group)) {
  dir_main_set=paste0(dir_main_all,group_type)
  setwd(dir_main_set)
  infercnv_obj = readRDS("./run.final.infercnv_obj")
  expr <- infercnv_obj@expr.data
  normal_loc <- infercnv_obj@reference_grouped_cell_indices
  normal_loc <- normal_loc$cCON
  test_loc <- infercnv_obj@observation_grouped_cell_indices
  test_loc<-unlist(test_loc)
  pvalue_signif_cutoff=0.01
  expr=expr-1
  expr_normal=expr[,colnames(expr)[normal_loc]]
  expr_tumor=expr[,colnames(expr)[test_loc]]
  expr_normal=annotateGenes.hg20(expr_normal, ID.type = "S")
  expr_tumor=annotateGenes.hg20(expr_tumor, ID.type = "S")
  dat_fin=data.frame()
  for (i in 1:nrow(expr_tumor)) {
    bin_normal=as.numeric(expr_normal[i,8:ncol(expr_normal)])
    bin_tumor=as.numeric(expr_tumor[i,8:ncol(expr_tumor)]) 
    normal_mean = mean(bin_normal)
    normal_sd = sd(bin_normal)
    cut_top=qnorm((1-pvalue_signif_cutoff),mean=normal_mean, sd=normal_sd)
    cut_floor=qnorm((pvalue_signif_cutoff),mean=normal_mean, sd=normal_sd)
    Amp_pos=bin_tumor>cut_top
    Del_pos=bin_tumor<cut_floor
    bin_tumor[Amp_pos]<-"Amp"
    bin_tumor[Del_pos]<-"Del"
    bin_tumor[bin_tumor!="Amp"&bin_tumor!="Del"]<-"nosig"
    dat_fin=rbind(dat_fin,bin_tumor)
      }
  colnames(dat_fin)=colnames(expr_tumor[i,8:ncol(expr_tumor)])
  Amp_freq= apply(dat_fin, 1, function(x) sum(x=="Amp"))/ncol(dat_fin)
  Del_freq= -apply(dat_fin, 1, function(x)sum(x=="Del"))/ncol(dat_fin)
  data_all=cbind(expr_tumor[,1:7],Amp_freq,Del_freq)
  chrom=paste0("Chr",data_all$chromosome_name)
  data_all$chrom=chrom
  data_all$chrom=factor(data_all$chrom,levels = paste0("Chr",1:22))
  for (chorm in unique(data_all$chrom)) {
    data_all[data_all$chrom==chorm,]$abspos=data_all[data_all$chrom==chorm,]$abspos-min(data_all[data_all$chrom==chorm,]$abspos)
  }
  write.csv(data_all,"data_all_long_plot.csv")
  p<-ggplot(data_all) +     
    geom_point(mapping=aes(x=abspos, y=Amp_freq), colour="red", size=0.1) + 
    geom_line(mapping=aes(x=abspos, y=Amp_freq), colour="red") +
    geom_point(mapping=aes(x=abspos, y=Del_freq), colour="blue", size=0.1) + 
    geom_line(mapping=aes(x=abspos, y=Del_freq), colour="blue") +     
    facet_grid(~chrom, space="free", scale="free") +
    geom_hline(yintercept=0, colour="black") +    
    geom_vline(data=data.frame(f=2, x=1), mapping=aes(xintercept=x),
               xintercept=0, 
               linetype=2, size=0.5, col="black") + 
    xlab("Chromosome / Position") +
    theme_bw() +    
    theme(panel.spacing = unit(0, "mm"),          
          panel.border = element_blank(),           
          panel.background = element_blank(),           
          panel.grid = element_blank(),          
          strip.background = element_blank(),          
          strip.text.x = element_text(colour="black",size=8,face="plain", angle=90),          
          axis.text.x=element_blank(),          
          axis.text.y = element_text(colour="black",size=8,face="plain"),          
          axis.ticks.x = element_blank(),          
          axis.line.x=element_blank(),          
          axis.title = element_text(colour="black",size=8,face="plain"),          
          plot.title = element_text(colour="black",size=8,face="plain"))   
  ggsave(paste0("cnv_type1_",group_type,".pdf"), plot = p, width = 12, height = 2.5)
  p<-ggplot(data_all) +     
    geom_area(mapping=aes(x=abspos, y=Amp_freq),fill="#EE0000")+
    geom_area(mapping=aes(x=abspos, y=Del_freq),fill="#00448C")+
    facet_grid(~chrom, space="free", scale="free") +
    geom_vline(data=data.frame(f=2, x=1), mapping=aes(xintercept=x),
               xintercept=0, 
               linetype=2, size=0.5, col="black")+
    ylim(-0.5,0.5)+
    theme_minimal()+
    xlab("Chromosome / Position") + ylab("Cell Freq") + 
    theme(panel.spacing = unit(0, "mm"),          
          panel.border = element_blank(),           
          panel.background = element_blank(),           
          panel.grid = element_blank(),          
          strip.background = element_blank(),          
          strip.text.x = element_text(colour="black",size=8,face="plain", angle=90),          
          axis.text.x=element_blank(),          
          axis.text.y = element_text(colour="black",size=8,face="plain"),          
          axis.ticks.x = element_blank(),          
          axis.line.x=element_blank(),          
          axis.title = element_text(colour="black",size=8,face="plain"),          
          plot.title = element_text(colour="black",size=8,face="plain")) #,hjust = 0.5
  ggsave(paste0("cnv_type2_",group_type,".pdf"), plot = p, width = 12, height = 2.5)
  p<-ggplot(data_all) +     
       geom_line(mapping=aes(x=abspos, y=Amp_freq), colour="red",size=0.7) +
    geom_area(mapping=aes(x=abspos, y=Amp_freq),fill="#EE0000",alpha=0.4)+
    geom_line(mapping=aes(x=abspos, y=Del_freq), colour="#00448C",size=0.7) +  
    geom_area(mapping=aes(x=abspos, y=Del_freq),fill="#00448C",alpha=0.4)+
    facet_grid(~chrom, space="free", scale="free") +
    ylim(-0.9,0.9)+
    geom_hline(yintercept=0, colour="black") +    
    geom_vline(data=data.frame(f=2, x=1), mapping=aes(xintercept=x),
               xintercept=0, 
               linetype=2, size=0.5, col="black") + 
    xlab("Chromosome / Position") + #ylab(ylabel) +     
    theme_bw() +    
    theme(panel.spacing = unit(0, "mm"),          
          panel.border = element_blank(),           
          panel.background = element_blank(),           
          panel.grid = element_blank(),          
          strip.background = element_blank(),          
          strip.text.x = element_text(colour="black",size=8,face="plain", angle=90),          
          axis.text.x=element_blank(),          
          axis.text.y = element_text(colour="black",size=8,face="plain"),          
          axis.ticks.x = element_blank(),          
          axis.line.x=element_blank(),          
          axis.title = element_text(colour="black",size=8,face="plain"),          
          plot.title = element_text(colour="black",size=8,face="plain")) 
  ggsave(paste0("cnv_type3_",group_type,".pdf"), plot = p, width = 12, height = 2.5)
  
}



library(ggplot2)
library(ggrastr)
library(tidyr)
library(scatterpie)
#data.combined_raw=data.combined
#################################################################
cnv专用



infercnv_obj = readRDS("./run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
expr<-expr[,test_loc]
expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)
#kmeans_df_s$CB=rownames(kmeans_df_s)
meta=data.combined_fib@meta.data
meta$CB=rownames(meta)
CNV_score=CNV_score%>%inner_join(meta,by="CB")


###########################################################################
meta=CNV_score

df <- data.combined@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = data.combined@meta.data$fincell)
colnames(df)
df2=df
df2$CB=rownames(df2)
meta=meta%>%inner_join(df2,by="CB")
score=quantile(meta$CNV_score)
cnv_type=c()
for (i in 1:nrow(meta)) {
  tmp=meta$CNV_score[i]
  if (tmp>score[4]|tmp==score[4]) {
    cnv_type[i]="CNV_D"
  }
  if (tmp>score[3]&tmp<score[4]|tmp==score[3]) {
    cnv_type[i]="CNV_C"
  }
  if (tmp>score[2]&tmp<score[3]|tmp==score[2]) {
    cnv_type[i]="CNV_B"
  }
  if (tmp>score[1]&tmp<score[2]|tmp==score[1]) {
    cnv_type[i]="CNV_A"
  }
  
  
}
meta$group=cnv_type
df3 <- meta
Cellratio <- prop.table(table(meta$group,meta$fincell), margin = 2)
Cellratio <- as.data.frame(Cellratio)
library(tidyr)
freq <-spread(Cellratio, Var1, Freq)
colnames(freq)[1] <- 'celltype'
freq <- freq[sort(freq$celltype),]
label <- df %>%group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))%>%
  as.data.frame()
rownames(label) <- label$cell_type
cell_number <- as.data.frame(table(data.combined$fincell))
colnames(cell_number)[2]<-'cellnumber'
cell_number$cellnumber <- log2(cell_number$cellnumber)/15

data = cbind(freq,label[,c(2:3)], cell_number[,c(2)])
colnames(data)[length(colnames(data))]<- 'cellnumber'
ggplot()+
  geom_point_rast(data=meta, aes(x= UMAP_1 , y = UMAP_2 ,color = group),size = 0.8,shape=16) +
  scale_color_manual(values = c("#71BC7A","#7DC6C5","#FFFA94","#FF935F"))+ 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_scatterpie(data=data,
                     color = "grey50",
                  aes(x=UMAP_1,y=UMAP_2,
                      group=celltype,
                      r=cellnumber),
                  cols=names(freq)[2:(1+length(unique(meta$group)))]
  )+
  scale_fill_manual(values = c("#71BC7A","#7DC6C5","#FFFA94","#FF935F"),name='cnv')
ggsave("scatterpie_cnv3.pdf", width = 6, height = 5)

#figS1b
celltype=c("T_NK","T_prolif",
           "B_cell","pB",
           "Myeloid",
           "EC",
           "Fib","myFib",
           "Normal_epi",
           "Tumor","Tumor_emt"
)

col=c("#DE4C00","#F47C30",
      "#ADD167","#79992F",
      "#73398D",
      "#FB7D80",
      "#FFB95C","#D9A833",
      "#608CC8",
      "#BBD4E0","#8CB4D2")

celltype_col=col
names(col)= celltype
names(celltype_col)= celltype

pdf("Overview_Cluster.pdf",11,4.6)
DimPlot(data.combined, reduction = 'umap', label=F,raster=FALSE,cols =col ) 
dev.off()


for (name in unique(data.combined$group)) {
  meta=data.combined@meta.data
  meta_sub=meta[meta$group==name,]
  fib.cell1 <- rownames(meta_sub)
  sub_CT<- subset(data.combined, name %in% fib.cell1)
  data.combined1=sub_CT
  data.combined1$fincell=Idents(data.combined1)
  data.combined5=data.combined1
  library(LSD)
  pdf(paste0("heatscatter","_",name,".pdf"),4,4.5)
  heatscatter(data.combined5@reductions$umap@cell.embeddings[,1],
              data.combined5@reductions$umap@cell.embeddings[,2])
  dev.off()
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggpointdensity) 
  sce_cca=data.combined5
  data <- cbind(Embeddings(object=sce_cca[['umap']]),
                FetchData(sce_cca,c("fincell","group")))
  ggplot(data = data, mapping = aes(x = UMAP_1,y = UMAP_2)) + 
    geom_pointdensity(adjust=6) + 
    geom_density_2d(bins = 5, colour="black")+
    xlim(min(data.combined@reductions$umap@cell.embeddings[,1]),max(data.combined@reductions$umap@cell.embeddings[,1]))+
    ylim(min(data.combined@reductions$umap@cell.embeddings[,2]),max(data.combined@reductions$umap@cell.embeddings[,2]))+
    theme_bw()+
    theme(
      axis.text=element_blank(),
      axis.ticks = element_blank(),
      axis.title=element_blank(),
      legend.text = element_text(size =10),
      aspect.ratio=1,
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.x = element_text(size=12, color = "black",
                                  vjust = 0.5,margin = margin(b = 3,t=3)), 
      strip.background = element_blank(),
      plot.margin=unit(c(1, 1, 1, 1),'cm'),
      legend.margin = margin(-0.2,-0.2,0,0,'cm'),
      legend.key.height = unit(1,'cm'),
      legend.key.width = unit(0.4,'cm'))+
                          guide = guide_colorbar(frame.colour = "black", 
                                                 ticks.colour = NA),
                          name = "Density",
                          labels = c('low',"high"),
                          breaks = c(100,700))   
  ggsave(paste0("density","_",name,".pdf"), width = 4, height = 4)	
}

#figS1b

heatmap_gene <- c("CD3D", "CD3E","MKI67","CD79A","CD79B","IGHG2",
                  "CD68","C1QA","CD1C","KIT","GATA2",
                  "VWF","CDH5",
                  "COL6A3","DCN","MYH11","NRXN1","CDH19",
                  "KRT3","KRT5"
                  
)
data.combined@meta.data$cluster <- factor(data.combined@meta.data$cluster, levels = cluster_sort_by_celltype)
Idents(data.combined)="celltype"
heatmap_AveE <- AverageExpression(data.combined, assays = "RNA", features = heatmap_gene,verbose = TRUE) %>% .$RNA
gaps_row <- cumsum(gene_num)
gaps_col <- cumsum(cluster_num)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
annotation_row <- data.frame(row.names = rownames(heatmap_AveE),
                             `GeneType` = rep(as.factor(c(1:length(gene_num))),gene_num))
annotation_col <- data.frame(row.names = colnames(heatmap_AveE),
                             `CellType` = colnames(heatmap_AveE))
annotation_colors = list(`CellType` = celltype_col)
pdf("plot.pdf",8,5)
pheatmap(heatmap_AveE,cluster_cols = F,cluster_rows = F,show_colnames=F,show_rownames=T,
                border=F,#border_color = "white",
                color = c(colorRampPalette(colors = c("#2166ac","#f7fbff"))(length(bk)/2),
                          colorRampPalette(colors = c("#f7fbff","#b2182b"))(length(bk)/2)),
                breaks=bk,scale="row",legend_breaks=seq(-2,2,2),
                gaps_row = gaps_row,
                gaps_col = gaps_col,
                #annotation_row = annotation_row,
                annotation_col = annotation_col,
                annotation_colors = annotation_colors,
                annotation_names_row = F,annotation_names_col = T)
dev.off()

library(reshape2)
library(ggpubr)
pB2_df <- table(data.combined@meta.data$cluster,data.combined@meta.data$group) %>% melt()
colnames(pB2_df) <- c("Cluster","Sample","Number")
pB2_df$Cluster <- factor(pB2_df$Cluster,levels = cluster_sort_by_celltype)
for (i in unique(pB2_df$Sample)) {
  pB2_df[pB2_df$Sample==i,]$Number=pB2_df[pB2_df$Sample==i,]$Number/sum(pB2_df[pB2_df$Sample==i,]$Number)
}
sample_color <- c("#1CC5FE","#FB7D80")
pB2 <- ggplot(data = pB2_df, aes(x = Cluster, y = Number, fill = Sample)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=sample_color) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  rotate_x_text()+
  scale_y_continuous(position = "right")+
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))
stat_raw=data.frame(data.combined$celltype,data.combined$time)
table1=table(stat_raw$data.combined.celltype,stat_raw$data.combined.time)
table2=table1[,-5]
num_neot=sum(table1[,5])
num_other= sum(table2)  
ratio_neo=table1[,5]/num_neot
ratio_other=table2/num_other
fin_ratio=ratio_neo/ratio_other
AveExpression=data.frame(fin_ratio)
AveExpression$Gene <- "fin_ratio"
AveExpression$Cluster=rownames(AveExpression)
Ave_df=AveExpression
colnames(Ave_df) <- c("Expression", "Gene", "Cluster")
Ave_df$Group <- paste(Ave_df$Cluster,Ave_df$Gene,sep = "_")
Ave_df$Expression[which(Ave_df$Expression>3)] <- 3 Ave_df$Cluster=factor(Ave_df$Cluster,levels = celltype)
pB1 <- ggplot(Ave_df,aes(x=Gene, y=Expression))+
  geom_hline(yintercept = seq(0, 2, 0.5),linetype = 2, color = "lightgray",size=1)+
  geom_line()+
  geom_segment(aes(x=Gene,xend=Gene,y=0,yend=Expression),color="lightgray",size = 1.5)+
  geom_point(size=4,aes(color=Gene))+
  scale_color_manual(values=c("#FB7D80","#AF1E1F")) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="ratio")
pB1 <- facet(pB1, facet.by = "Cluster",ncol = length(unique(Ave_df$Cluster)),panel.labs.font = list(size = 12),panel.labs.background = list(fill = "#a6cee3"))
pB1 <- pB1 + scale_y_continuous(position = "right")+ 
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_blank())+ 
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "right",
        panel.border = element_blank(),
        axis.ticks.x = element_line(color =  NA))
require(ggplotify)
pB_1_2 <- pB1 + pB2 + plot_layout(ncol = 1, heights = c(1.2, 2))
ggsave(paste0("plot.pdf"), plot = pB_1_2, width = 8, height =5)

#figS1c
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
pdf("Roe1.pdf",width = 3,height = 6)
Heatmap(as.data.frame(Roe),
        show_heatmap_legend = TRUE, 
        cluster_rows = F, 
        cluster_columns = F, 
        row_names_side = 'right', 
        column_names_side = "top",
        show_column_names = TRUE,
        show_row_names = TRUE,
        #row_split = 3, #行分组
        col = col_fun,
        row_names_gp = gpar(fontsize = 12), 
        heatmap_legend_param = list(
          title = "Ro/e value",  
          at = c(0,max(Roe)), 
          labels = seq(0,max(Roe)) 
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black")) )
dev.off()

#figS1e,d
library(tidyverse)
meta = meta_all
meta$Patient_time_origin=paste0(meta$stim,"-",meta$group)
meta_sub=meta
cellratio <- prop.table(table(meta_sub$fincell, meta_sub$Patient_time_origin), margin = 2)%>%
  as.data.frame()%>%
  `colnames<-`(c('celltype','orig.ident','Freq'))
cellratio = cellratio %>% 
  separate(orig.ident, into = c('Patient','origin'), sep = '-')
cellratio=cellratio[order(cellratio$celltype),]
cellratio_a=cellratio[cellratio$origin=="A_post",]
cellratio_u=cellratio[cellratio$origin=="B_post",]
cellratio_a<-
  cellratio_a%>%group_by(celltype)%>%summarise(mean_Freq=mean(Freq))
cellratio_u<-
  cellratio_u%>%group_by(celltype)%>%summarise(mean_Freq=mean(Freq))
ratio=cellratio_u$mean_Freq/cellratio_a$mean_Freq
cellratio_u$ratio=ratio
cellratio1 <- as.data.frame(cellratio_u)
cellratio1$celltype=as.vector(cellratio1$celltype)
nmf_type=c(rep("NMF3",length(gro[[1]])),
           rep("NMF4",length(gro[[2]])),
           rep("NMF1",length(gro[[3]])),
           rep("NMF2",length(gro[[4]])),
           rep("NMF3",length(gro[[5]]))
           )

tmp=unlist(gro)
names(tmp)=NULL
names(nmf_type)=tmp
cellratio1=cellratio1[cellratio1$celltype%in%tmp,]
cellratio1$nmf_type=nmf_type[cellratio1$celltype]
cellratio1$ratio=log2(cellratio1$ratio)
cellratio1=cellratio1[order(cellratio1$ratio,decreasing = T),]
cellratio1$nmf_type=factor(cellratio1$nmf_type,levels = c("NMF1","NMF2","NMF3","NMF4"))
cellratio1=cellratio1[order(cellratio1$nmf_type),]
cellratio1$celltype=factor(cellratio1$celltype,levels = unique(cellratio1$celltype))
nmf_cols=c("#FF7F00", "#1F78B4", 
           "#33A02C", "#6A3D9A", 
           "#E31A1C")

nmf_cols=c("#69c0a2","#D96669","#FF935F","#E29DCE","#69c0a2","#8CB4D2")
nmf_cols=c("#D96669","#FF935F","#69c0a2","#E29DCE","#69c0a2","#8CB4D2")

ggplot(data = cellratio1, mapping = aes(x = celltype, y = ratio,
                                        color =nmf_type, fill = nmf_type)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values=nmf_cols)+
  scale_color_manual(values=nmf_cols)+
  geom_bar(stat = "identity", position=position_dodge())

#figS1f

df_igraph <- graph_from_data_frame(edge_list,node, directed = FALSE)
ggraph(df_igraph, layout = "linear",circular = TRUE)+
  geom_edge_fan(aes(color=color,edge_width=weight,edge_alpha = weight),show.legend = T)+
  geom_node_point(aes(fill=size,color=size),size=20,shape=21)+
  geom_node_text(aes(label=name, size = weight),
                  angle=0,hjust=0.1,size=3) +
  scale_edge_width_continuous(range = c(1,3),
                              limits = c(0.01,1),
                   guide = guide_legend(title = "Jaccard index",
                                        order=1))+
  scale_alpha_continuous(range = c(0.4,0.7))+
  scale_color_gradientn( colours =colorRampPalette(c("#FFED99","#FF7600"))(10))+
  scale_fill_gradientn(colours =colorRampPalette(c("#FFED99","#FF7600"))(10),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"))+
  scale_edge_colour_manual(values="grey")+
  theme_graph()+
  theme(legend.position = 'none')+
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2))




