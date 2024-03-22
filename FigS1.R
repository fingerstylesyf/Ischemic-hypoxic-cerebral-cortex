library(Matrix)
library(Seurat)
library(tidyverse)
library(liger)
library(patchwork)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(DoubletFinder)
library(Matrix)
library(magrittr)
library(stringr)
library(tidyverse)
library(scRNAtoolVis)

load("AB_RENAME.RData")

#---------UMAP(Cluster)
p7 <- DimPlot(AB, reduction = "umap", group.by = "age", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_Cluster.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

#---------QC
VlnPlot(AB, group.by = "age",features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), cols = c("#ECE6A2","#AFB5C5","#b4b47b","#F5BEC1","#74AED4","#80C684","#CFAFD4","#E2AA4C"),pt.size = 0 ,ncol = 3)
ggsave("ALL_control(no dot).pdf", width = 70, height = 35, units = "cm")

#-----------Cluster-Dotplot
jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Myt1","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Mbp","Plp1"),
                    gene.order = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Myt1","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Mbp","Plp1"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot-cluster.pdf", plot = jjplot, device = 'pdf', width = 32, height = 23, units = 'cm')

####highlight
p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "Sham.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "Sham3d.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2),
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "Sham7d.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "1h.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "4h.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "12h.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "3d.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                ExN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents="ExN"),
                InN=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='InN'),
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='CR'),
                AST=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST'),
                OPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='OPC'),
                OL=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='OL'),
                MG=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG'),
                Endo=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Endo'),
                Peri=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Peri'),
                VLMC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='VLMC')
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#81A88D', '#E86A52', '#FFE27E','#FDBF6F','#E17D8E', '#FBB4AE', '#6E90D0', '#AE7669',   "#6f4fa3","#DECBE4") 
)
p1
ggsave(filename = "7d.pdf", plot = p1, device = 'pdf', width = 18.5, height = 14.5, units = 'cm')

