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

#figS3a
#refer to fig4d

#figS3b
#refer to fig4e

#figS3c,d
#refer to fig4f

#figS3e
#refer to figS2f

#figS3f
#refer to figS2i

#figS3g

library(VISION)
sigScores <- as.matrix(t(getSignatureScores(vis)))
path=rownames(as.matrix(t(getSignatureScores(vis))))

library(VISION)
monocle_meta=monocle_cds@phenoData@data
sigScores <- as.matrix(t(getSignatureScores(vis)))
fib.cell=intersect(colnames(monocle_cds),colnames(sigScores))
sub_CT<- subset(data.combined, name %in% fib.cell)
data.combined=sub_CT
data.combined$fincell=Idents(data.combined)
monocle_meta=monocle_meta[colnames(data.combined),]
data.combined$Pseudotime=monocle_meta$Pseudotime
sigScores <- as.matrix(t(getSignatureScores(vis)))
sigScores=sigScores[,colnames(data.combined)]
plot.dat=data.combined@meta.data

for (i in path) {
  gene=sigScores[i,]
  plot.dat$gene=gene
  p1<-ggplot(data=plot.dat,aes(x=Pseudotime,y=gene))+
    geom_smooth(aes(color=group), 
                span = 1, size = 0.7,
                method = "loess",
                level=0.65,
                formula = 'y~x')+
    scale_color_manual(values=sample_color)+
    theme_bw(base_size = 15)+
    theme(panel.grid = element_blank())+
    labs(y="Expression")
  ggsave(plot=p1,filename=paste0(i,".pdf"),width = 4.8, height = 3.5)
}