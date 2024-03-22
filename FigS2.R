#-----------------Astrocytes-----------------#
library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(scRNAtoolVis)
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
library(Scillus)

####highlight
p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents="AST7"),
                AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#2b5aab","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "Sham.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents="AST7"),
                #AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "Sham3d.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents="AST7"),
                #AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "Sham7d.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')




p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents="AST7"),
                AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#2b5aab","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "1h.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents="AST7"),
                AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#2b5aab","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "4h.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents="AST7"),
                AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='AST1')
                #EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#c7ad6c","#03776a","#2b5aab","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "12h.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents="AST7"),
                AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#2b5aab","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "3d.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.8,cols = "grey85",
              cells.highlight=list(
                AST8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST8'),
                AST7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents="AST7"),
                AST6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST6'),
                AST5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST5'),
                AST4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST4'),
                AST3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST3'),
                AST2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST2'),
                AST1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='AST1'),
                EPC=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='EPC')
                
              ),
              sizes.highlight = c(0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#0EB8C9","#c7ad6c","#03776a","#2b5aab","#Afd2ce","#71CB89","#AB4361","#6F7796","#fbd15c") 
              
)
p1
ggsave(filename = "7d.png", plot = p1, device = 'png', width = 13, height = 9.6, units = 'cm')


