library(Seurat)
library(CellChat)
library(ggplot2)
library(patchwork)

#figS5a
#refer to fig4d

#figS5b
#refer to fig6j

#figS5c
#refer to fig6j


tile <- read.delim("./tile.txt", row.names=NULL)
p_tile<-ggplot(data=tile,aes(x=L,y= R)) + 
  geom_tile(data=tile,aes(x=L,y= R,  fill=prob),color="white",size=0)+
  scale_fill_gradient2(high="#9733CF", low="white")+#对应修改有填充的点颜色
   theme_bw()+
  theme(legend.position = "none",
  plot.margin = margin(t = 0,  
                         r = 0,  
                         b = 0, 
                         l = 0,  
                         unit = "cm"),
    panel.grid = element_blank(),
    axis.text.x=element_text(angle=90,hjust = 0,vjust=1),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title=element_text(vjust = 4, hjust = 0.5,size=12))+
  coord_cartesian(clip = 'off')
p_tile
library(Seurat)
data.combined=data.combined
Idents(data.combined)<-"group"
HC <-data.combined
ligand_tile=unique(tile$R)
ligand_tile=sort(ligand_tile)
ligand_dot <- DotPlot(HC, features = ligand_tile,assay = "RNA")
ligand_exp <- ligand_dot$data
unique(ligand_exp$id)
ligand_exp$id <- factor(ligand_exp$id, levels = group_level)
p2_2 <- ggplot(ligand_exp,aes(x=id,y= features.plot))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x=element_text(angle=90,hjust = 0,vjust=1,colour = 'black',size = 9),
        axis.text.y = element_blank(),
        plot.margin = margin(t = 0,  
                             r = 0,
                             b = 0, 
                             l = 0,  
                             unit = "cm"),
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 5))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours =colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(title = "",order = 1,nrow = 2))+
  coord_cartesian(clip = 'off') 
p4_2 <- p2_2
p4_2
library(Seurat)
data.combined=data.combined_fib

Idents(data.combined)<-"group"
immune=data.combined
target_fin=tile$L
target_fin=unique(target_fin)
target_tile=sort(target_fin)
#immune <- subset(HC, celltype%in% c("Macrophages", "Lymphocytes"))
receptor_dot <- DotPlot(immune, features = target_tile,assay = "RNA")
receptor_exp <- receptor_dot$data
p5_2 <- ggplot(receptor_exp,aes(x=features.plot,y=id ))+
  geom_point(aes(size=`pct.exp`,
                 color=`avg.exp.scaled`))+
  theme_bw()+
  scale_x_discrete(position = "top")+
  scale_y_discrete(position = "left")+
  theme(axis.text.x=element_blank(),
        axis.text.y = element_text(angle=0,size = 9),
        title = element_text(colour = 'black',size = 8),
        plot.title = element_text(hjust=0.5),
        plot.margin = margin(t = 0,  
                             r = 0, 
                             b = 0,  
                             l = 0, 
                             unit = "cm"),
        panel.grid = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.3,'cm'),
        legend.key.width = unit(0.3,'cm'),
        legend.title = element_text(size=8,vjust = 1),
        legend.text = element_text(size = 6))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours =colorRampPalette(c("#1E3163","#00C1D4","#FFED99","#FF7600"))(10),
                       # colours = c('#1A5592','white',"#B83D3D"),
                        guide = guide_colorbar(ticks.colour = "black",
                                               frame.colour = "black"),
                        name = "Scaled expr")+
  labs(x=NULL,y=NULL)+
  coord_cartesian(clip = 'off') 
p7_2 <- p5_2
p7_2
