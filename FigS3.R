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

###AST4
AST4 <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("AST4")]
cr <- c(paste('AST4 ',c('(Sham)','(Sham3d)','(Sham7d)','(1h)','(4h)','(12h)','(3d)','(7d)'),sep = ''))

jjplot <- jjDotPlot(object = AST4,
                    gene = c("Actn1","Ltbp1","Adam12","Arntl2","Kcnt1","Sema3c", "Cdh2","Camk2d","Pde10a","Vav3","Ano6","Pcbp3","Hmga2","Met","Cd44","Iqgap2","Sbno2","Dclk1","Gpd2","Rap1gds1","Dlg1","Piezo1","Osmr","Nav2","Fosl1","Capn2","Nfatc2"),
                    id = 'celltype',
                    split.by = 'age',
                    split.by.aesGroup = T,dot.col = c("#2F4858",'white',"#02C837"),
                    point.geom = F,
                    tile.geom = T,ytree = F,cluster.order = cr,
                    col.min = -0.5,col.max = 1.5)+ coord_flip()

jjplot
ggsave(filename = "AST4.pdf", plot = jjplot, device = 'pdf', width = 13.9, height = 18, units = 'cm')

###AST6
AST6 <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("AST6")]
cr <- c(paste('AST6 ',c('(Sham)','(1h)','(4h)','(12h)','(3d)','(7d)'),sep = ''))

jjplot <- jjDotPlot(object = AST6,
                    gene = c("Nr4a2","Nr4a3","Atf3","Timp3","Fos","Fosb","Vegfa","Hmgcr","Hspa5","Stat3","Hsph1","Nfatc2","Sik1","Mest","Aff1","Ell2","Hsp90aa1","Btg2","Jun","Nr4a1","Egr1"),
                    id = 'celltype',
                    split.by = 'age',
                    split.by.aesGroup = T,dot.col = c("#FFEDCB",'white',"#2763A2"),
                    point.geom = F,
                    tile.geom = T,ytree = F,cluster.order = cr,
                    col.min = -0.4,col.max = 1.6)+ coord_flip()

jjplot
ggsave(filename = "AST6.pdf", plot = jjplot, device = 'pdf', width = 13.9, height = 18, units = 'cm')


