library(Matrix)
library(Seurat)
library(tidyverse)
library(rliger)
library(SeuratDisk)
library(SeuratWrappers)
library(patchwork)
library(cowplot)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(DoubletFinder)
library(magrittr)
library(stringr)
library(tidyverse)
library(scRNAtoolVis)
library(homologene)
#PCA(up40gene)
scRNA <- NormalizeData(scRNA, verbose = FALSE)
genes_of_interest <- c("Txn1", "Ctsb", "Srgn", "Cd72", "Plin2", "Slfn2", "Fth1", "C3ar1", "Atp6v0e", "Cd68",
                       "Glipr1", "Cd84", "Tpm4", "Prdx1", "Emp3", "Lgals3", "Tpd52", "Pkm", "Anxa5", "Sdc4",
                       "Gadd45b", "Dab2", "Ccl9", "Id2", "Cadm1", "Ccl4", "Ccl3", "Cd9", "Cstb", "Slc11a1",
                       "Cd63", "Plek", "Capg", "Cybb", "Bcl2a1d", "Cd14", "Lilr4b", "Vim", "Cxcl16", "Bcl2a1b")
data_to_export <- FetchData(scRNA, vars = c("group", genes_of_interest))
library(dplyr)
data <- data_to_export
data <- data %>%
  mutate(Sample_id = rownames(.)) %>%
  dplyr::select(Sample_id, everything())
rownames(data) <- NULL
names(data)[2] <- "Subtype"
write_tsv(data, file = "output_file(up40gene).xls")
df <- data
df <- read_tsv("output_file(up40gene).xls")
pca <- df %>% column_to_rownames(var="Sample_id") %>% 
  select(-Subtype) %>% prcomp(.,scale. = TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
df$Subtype <- factor(df$Subtype, levels = c("Sham", "S3h", "S4h", "S12h", "S1d", "S2d", "S3d", "S7d", "S14d"))
colors <- c("#C0BC84",'#FFAF88','#FFA1DC',"#FF7D8C","#D8BFD8","#A8C4C2",'#00bf72',"#339793","#2781B6")
p1 <- fviz_pca_biplot(pca, axes = c(1, 2),geom.ind = c("point"),geom.var = c("point", "text"),alpha.var = 0.8,
                      pointshape = 20,pointsize=4,label ="var",repel = TRUE,col.var = "black",
                      labelsize=0.5,addEllipses=TRUE, ellipse.level=0.95,
                      col.ind = df$Subtype)+
  scale_color_manual(values = colors)+
  labs(x=paste0("(PC1: ",round(var_explained[1]*100,2),"%)"),
       y=paste0("(PC2: ",round(var_explained[2]*100,2),"%)"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour="black",size = 12,margin = margin(t=12)),
        axis.title.y = element_text(colour="black",size = 12,margin = margin(r=12)),
        axis.text=element_text(color="black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),   
        legend.text = element_text(color="black",size=9), 
        legend.spacing.x=unit(0.06,'cm'), 
        legend.key.width=unit(0.01,'cm'), 
        legend.key.height=unit(0.01,'cm'), 
        legend.background=element_blank(), 
        legend.position=c(1,0),legend.justification=c(1,0))
p1
ggsave(filename = "PCA(up40gene).pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')

#PCA(down21gene)
scRNA <- NormalizeData(scRNA, verbose = FALSE)
genes_of_interest <- c("St3gal6", "Hexb", "Nrip1",  "Ivns1abp", "Susd3",  "Lpcat2", "F11r", "Vsir", "Arhgap5",  "Mef2a",  "P2ry13", "Maf",  "Sparc",  "Mef2c",  "Rnase4", "Abi3", "Selplg", "Tmem119",  "Gpr34",  "Mtus1",  "Siglech")
data_to_export <- FetchData(scRNA, vars = c("group", genes_of_interest))
library(dplyr)
data <- data_to_export
data <- data %>%
  mutate(Sample_id = rownames(.)) %>%
  dplyr::select(Sample_id, everything())
rownames(data) <- NULL
names(data)[2] <- "Subtype"
write_tsv(data, file = "output_file(down21gene).xls")
df <- data
df <- read_tsv("output_file(down21gene).xls")
pca <- df %>% column_to_rownames(var="Sample_id") %>% 
  select(-Subtype) %>% prcomp(.,scale. = TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)
df$Subtype <- factor(df$Subtype, levels = c("Sham", "S3h", "S4h", "S12h", "S1d", "S2d", "S3d", "S7d", "S14d"))
colors <- c("#C0BC84",'#FFAF88','#FFA1DC',"#FF7D8C","#D8BFD8","#A8C4C2",'#00bf72',"#339793","#2781B6")
p1 <- fviz_pca_biplot(pca, axes = c(1, 2),geom.ind = c("point"),geom.var = c("point", "text"),alpha.var = 0.8,
                      pointshape = 20,pointsize=4,label ="var",repel = TRUE,col.var = "black",
                      labelsize=0.5,addEllipses=TRUE, ellipse.level=0.95,
                      col.ind = df$Subtype)+
  scale_color_manual(values = colors)+
  labs(x=paste0("(PC1: ",round(var_explained[1]*100,2),"%)"),
       y=paste0("(PC2: ",round(var_explained[2]*100,2),"%)"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour="black",size = 12,margin = margin(t=12)),
        axis.title.y = element_text(colour="black",size = 12,margin = margin(r=12)),
        axis.text=element_text(color="black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),   
        legend.text = element_text(color="black",size=9), 
        legend.spacing.x=unit(0.06,'cm'), 
        legend.key.width=unit(0.01,'cm'), 
        legend.key.height=unit(0.01,'cm'), 
        legend.background=element_blank(), 
        legend.position=c(1,0),legend.justification=c(1,0))
p1
ggsave(filename = "PCA(down21gene).pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')

#PCA(up37gene)
genes_of_interest <- c("Hip1",  "Kif26b", "Slc44a2",  "Mitf", "Daglb",  "Fam241a",  "Gabbr2", "Agap1",  "Olr1", "Asap1",  "Rasal2", "Dhrs9",  "Frmd4b", "Rab31",  "Urb1", "Gab2", "Prkch",  "Map3k1", "Plekhm3",  "Nav3", "Dscam",  "Rasgrp3",  "Svil", "Arhgap24", "Lhfp", "Cass4",  "Wapl", "Abr",  "Niban2", "Niban1", "Ptpra",  "Gsg1l",  "Rnf150", "Pacsin1",  "Rerg", "Tec",  "Arsb")
data_to_export <- FetchData(scRNA, vars = c("group", genes_of_interest))
library(dplyr)
data <- data_to_export
data <- data %>%
  mutate(Sample_id = rownames(.)) %>%
  dplyr::select(Sample_id, everything())
rownames(data) <- NULL
names(data)[2] <- "Subtype"
write_tsv(data, file = "output_file(up37gene).xls")
df <- data
df <- read_tsv("output_file(up37gene).xls")
pca <- df %>% column_to_rownames(var="Sample_id") %>% 
  select(-Subtype) %>% prcomp(.,scale. = TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)

df$Subtype <- factor(df$Subtype, levels = c("Sham", "1h", "4h", "12h", "3d", "7d"))
colors <- c("#C0BC84",'#FFAF88',"#FF7D8C","#A8C4C2",'#00bf72',"#2781B6")
p1 <- fviz_pca_biplot(pca, axes = c(1, 2),geom.ind = c("point"),geom.var = c("point", "text"),alpha.var = 0.8,
                      pointshape = 20,pointsize=4,label ="var",repel = TRUE,col.var = "black",
                      labelsize=0.5,addEllipses=TRUE, ellipse.level=0.95,
                      col.ind = df$Subtype)+
  scale_color_manual(values = colors)+
  labs(x=paste0("(PC1: ",round(var_explained[1]*100,2),"%)"),
       y=paste0("(PC2: ",round(var_explained[2]*100,2),"%)"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour="black",size = 12,margin = margin(t=12)),
        axis.title.y = element_text(colour="black",size = 12,margin = margin(r=12)),
        axis.text=element_text(color="black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),   
        legend.text = element_text(color="black",size=9), 
        legend.spacing.x=unit(0.06,'cm'), 
        legend.key.width=unit(0.01,'cm'), 
        legend.key.height=unit(0.01,'cm'), 
        legend.background=element_blank(), 
        legend.position=c(1,0),legend.justification=c(1,0))
p1
ggsave(filename = "PCA(up37gene).pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')

#PCA(down14gene)
scRNA <- AB
genes_of_interest <- c("Gpr63", "Nfatc2", "Arhgap15", "Slc10a6",  "Zfp622", "C9orf72",  "Dtnb", "Adgre1", "Man1c1", "Fchsd2", "Prune2", "Fgfr1",  "Apobec1",  "Rapgef5")
data_to_export <- FetchData(scRNA, vars = c("group", genes_of_interest))
library(dplyr)
data <- data_to_export
data <- data %>%
  mutate(Sample_id = rownames(.)) %>%
  dplyr::select(Sample_id, everything())
rownames(data) <- NULL
names(data)[2] <- "Subtype"
write_tsv(data, file = "output_file(down14gene).xls")
df <- data
df <- read_tsv("output_file(down14gene).xls")
pca <- df %>% column_to_rownames(var="Sample_id") %>% 
  select(-Subtype) %>% prcomp(.,scale. = TRUE)

var_explained <- pca$sdev^2/sum(pca$sdev^2)

df$Subtype <- factor(df$Subtype, levels = c("Sham", "1h", "4h", "12h", "3d", "7d"))

colors <- c("#C0BC84",'#FFAF88',"#FF7D8C","#A8C4C2",'#00bf72',"#2781B6")
p1 <- fviz_pca_biplot(pca, axes = c(1, 2),geom.ind = c("point"),geom.var = c("point", "text"),alpha.var = 0.8,
                      pointshape = 20,pointsize=4,label ="var",repel = TRUE,col.var = "black",
                      labelsize=0.5,addEllipses=TRUE, ellipse.level=0.95,
                      col.ind = df$Subtype)+
  scale_color_manual(values = colors)+
  labs(x=paste0("(PC1: ",round(var_explained[1]*100,2),"%)"),
       y=paste0("(PC2: ",round(var_explained[2]*100,2),"%)"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour="black",size = 12,margin = margin(t=12)),
        axis.title.y = element_text(colour="black",size = 12,margin = margin(r=12)),
        axis.text=element_text(color="black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),   
        legend.text = element_text(color="black",size=9), 
        legend.spacing.x=unit(0.06,'cm'), 
        legend.key.width=unit(0.01,'cm'), 
        legend.key.height=unit(0.01,'cm'), 
        legend.background=element_blank(), 
        legend.position=c(1,0),legend.justification=c(1,0))
p1

ggsave(filename = "PCA(down14gene).pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')
#---------add neo
Q5 <- AB
Q5$group=str_replace(Q5$group,"1h","S1h")
Q5$group=str_replace(Q5$group,"4h","S4h")
Q5$group=str_replace(Q5$group,"12h","S12h")
Q5$group=str_replace(Q5$group,"3d","S3d")
Q5$group=str_replace(Q5$group,"7d","S7d")
###
gene_set1 <- data.frame(gene=row.names(Q1))
gene_set2 <- data.frame(gene=row.names(Q2))
gene_set3 <- data.frame(gene=row.names(Q3))
gene_set4 <- data.frame(gene=row.names(Q4))
gene_set5 <- data.frame(gene=row.names(Q5))
write_csv(gene_set5,file = "Rat.csv")

Rat_Mm = homologene(gene_set5$gene, inTax = 10116, outTax = 10090)
duplicated(Rat_Mm[,1])
Rat_Mm[duplicated(Rat_Mm[,1]) == F,]
Rat_Mm_single <- Rat_Mm[duplicated(Rat_Mm[,1]) == F,]

RenameGenesSeurat <- function(obj,newnames,gene.use=NULL,de.assay="RNA") { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

Q5 <- RenameGenesSeurat(obj=Q5,newnames=Rat_Mm_single[,2],gene.use=Rat_Mm_single[,1])
gene_set5 <- data.frame(gene=row.names(Q5))

gene_intersect <- Reduce(intersect, list(gene_set1[,1], gene_set2[,1], gene_set3[,1], gene_set4[,1], gene_set5[,1]))

Q1_data <- Q1@assays$RNA@data[Q1@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Q1@meta.data[,c("orig.ident","group")]
Q1 = CreateSeuratObject(counts = Q1_data, assay = 'RNA', meta.data = meta, project = "Q1")
Q1$Set="Q1"
Q1$Set.group <- paste(Q1$Set, Q1$group, sep = "_")

Q2_data <- Q2@assays$RNA@data[Q2@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Q2@meta.data[,c("orig.ident","group")]
Q2 = CreateSeuratObject(counts = Q2_data, assay = 'RNA', meta.data = meta, project = "Q2")
Q2$Set="Q2"
Q2$Set.group <- paste(Q2$Set, Q2$group, sep = "_")

Q3_data <- Q3@assays$RNA@data[Q3@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Q3@meta.data[,c("orig.ident","group")]
Q3 = CreateSeuratObject(counts = Q3_data, assay = 'RNA', meta.data = meta, project = "Q3")
Q3$Set="Q3"
Q3$Set.group <- paste(Q3$Set, Q3$group, sep = "_")

Q4_data <- Q4@assays$RNA@data[Q4@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Q4@meta.data[,c("orig.ident","group")]
Q4 = CreateSeuratObject(counts = Q4_data, assay = 'RNA', meta.data = meta, project = "Q4")
Q4$Set="Q4"
Q4$Set.group <- paste(Q4$Set, Q4$group, sep = "_")

Q5_data <- Q5@assays$RNA@data[Q5@assays$RNA@data@Dimnames[[1]] %in% gene_intersect,]
meta <- Q5@meta.data[,c("orig.ident","group")]
Q5 = CreateSeuratObject(counts = Q5_data, assay = 'RNA', meta.data = meta, project = "Q5")
Q5$Set="Q5"
Q5$Set.group <- paste(Q5$Set, Q5$group, sep = "_")

scRNA <- merge(Q1,y=c(Q2,Q3,Q4,Q5))
save(scRNA,file = "Merge(5Set).Rdata")
#---PCA(up5gene)
genes_of_interest <- c("Cd83",  "Lpl",  "Spp1", "Csf1", "Rab7b")
data_to_export <- FetchData(scRNA, vars = c("group", genes_of_interest))
library(dplyr)
data <- data_to_export
data <- data %>%
  mutate(Sample_id = rownames(.)) %>%
  dplyr::select(Sample_id, everything())
rownames(data) <- NULL
names(data)[2] <- "Subtype"
write_tsv(data, file = "output_file(up5gene).xls")
df <- data
df <- read_tsv("output_file(up5gene).xls")

pca <- df %>% column_to_rownames(var="Sample_id") %>% 
  select(-Subtype) %>% prcomp(.,scale. = TRUE)

var_explained <- pca$sdev^2/sum(pca$sdev^2)

df$Subtype <- factor(df$Subtype, levels = c("Sham", "S1h", "S3h", "S4h", "S12h", "S1d", "S2d", "S3d", "S7d", "S14d"))

colors <- c("#C0BC84",'#FFAF88','#FFA1DC',"#FF7D8C","#D8BFD8","#A8C4C2",'#00bf72',"#339793","#2781B6",'#6E71B9')
p1 <- fviz_pca_biplot(pca, axes = c(1, 2),geom.ind = c("point"),geom.var = c("point", "text"),alpha.var = 0.8,
                      pointshape = 20,pointsize=4,label ="var",repel = TRUE,col.var = "black",
                      labelsize=0.5,addEllipses=TRUE, ellipse.level=0.95,
                      col.ind = df$Subtype)+
  scale_color_manual(values = colors)+
  labs(x=paste0("(PC1: ",round(var_explained[1]*100,2),"%)"),
       y=paste0("(PC2: ",round(var_explained[2]*100,2),"%)"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour="black",size = 12,margin = margin(t=12)),
        axis.title.y = element_text(colour="black",size = 12,margin = margin(r=12)),
        axis.text=element_text(color="black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),   
        legend.text = element_text(color="black",size=9), 
        legend.spacing.x=unit(0.06,'cm'), 
        legend.key.width=unit(0.01,'cm'), 
        legend.key.height=unit(0.01,'cm'), 
        legend.background=element_blank(), 
        legend.position=c(1,0),legend.justification=c(1,0))
p1
ggsave(filename = "PCA(up5gene).pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')

#PCA(down9gene)
genes_of_interest <- c("Csf1r", "Srgap2", "P2ry12", "Cfh",  "Ifngr1", "Qk", "Cx3cr1", "Ptgs1",  "Mertk")
data_to_export <- FetchData(scRNA, vars = c("group", genes_of_interest))
library(dplyr)
data <- data_to_export
data <- data %>%
  mutate(Sample_id = rownames(.)) %>%
  dplyr::select(Sample_id, everything())
rownames(data) <- NULL
names(data)[2] <- "Subtype"
write_tsv(data, file = "output_file(down9gene).xls")
df <- data
df <- read_tsv("output_file(down9gene).xls")
pca <- df %>% column_to_rownames(var="Sample_id") %>% 
  select(-Subtype) %>% prcomp(.,scale. = TRUE)
var_explained <- pca$sdev^2/sum(pca$sdev^2)

df$Subtype <- factor(df$Subtype, levels = c("Sham", "S1h", "S3h", "S4h", "S12h", "S1d", "S2d", "S3d", "S7d", "S14d"))
colors <- c("#C0BC84",'#FFAF88','#FFA1DC',"#FF7D8C","#D8BFD8","#A8C4C2",'#00bf72',"#339793","#2781B6",'#6E71B9')
p1 <- fviz_pca_biplot(pca, axes = c(1, 2),geom.ind = c("point"),geom.var = c("point", "text"),alpha.var = 0.8,
                      pointshape = 20,pointsize=4,label ="var",repel = TRUE,col.var = "black",
                      labelsize=0.5,addEllipses=TRUE, ellipse.level=0.95,
                      col.ind = df$Subtype)+
  scale_color_manual(values = colors)+
  labs(x=paste0("(PC1: ",round(var_explained[1]*100,2),"%)"),
       y=paste0("(PC2: ",round(var_explained[2]*100,2),"%)"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.x = element_text(colour="black",size = 12,margin = margin(t=12)),
        axis.title.y = element_text(colour="black",size = 12,margin = margin(r=12)),
        axis.text=element_text(color="black"),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.key=element_blank(),   
        legend.text = element_text(color="black",size=9), 
        legend.spacing.x=unit(0.06,'cm'), 
        legend.key.width=unit(0.01,'cm'), 
        legend.key.height=unit(0.01,'cm'), 
        legend.background=element_blank(), 
        legend.position=c(1,0),legend.justification=c(1,0))
p1
ggsave(filename = "PCA(down9gene).pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')

