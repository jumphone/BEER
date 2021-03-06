library(Seurat)
library(cowplot)

##########
load('TSNE.RData')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }
sc_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$ori))    
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
D1=exp_sc_mat
dim(D1)
#19972  3005
###Load DevOlig#######
EXP=readRDS('GSE75330.RDS')
D2=as.matrix(EXP@raw.data)
dim(D2)
###############





ctrl.data<- D1
stim.data <- D2

ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
ctrl$stim <- "CTRL"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)
stim$stim <- "STIM"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)


DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

DimPlot(immune.combined, reduction = "umap", group.by = "stim")

immune.combined@meta.data$type=immune.combined@meta.data$stim

LABEL=readRDS('LABEL.RDS')
immune.combined@meta.data$type=LABEL
DimPlot(immune.combined, reduction = "umap", group.by = "type",label=T)






UTARGET_LABEL=LABEL
UTARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='Oligdend'
UTARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='Oligdend'

USE=which(!UTARGET_LABEL %in% c('NA'))
BATCH=UTARGET_LABEL[USE]

COM_UUMAP=immune.combined@reductions$umap@cell.embeddings[USE,]


library(vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3])


boxplot(out1[,3],ylim=c(-1,1), names=c('NEW'),las=2,main='Oligodend Merged' ,ylab='Silhouette Coefficient')


#########
UTARGET_LABEL=LABEL
UTARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='Oligdend'
UTARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='Oligdend'

USE=which(UTARGET_LABEL %in% c('astrocytes_ependymal_batch1','OPC_batch2','microglia_batch1'))
BATCH=UTARGET_LABEL[USE]

COM_UUMAP=immune.combined@reductions$umap@cell.embeddings[USE,]
library(vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3])
boxplot(out1[,3],ylim=c(-1,1), names=c('NEW'),las=2,main='Oligodend Merged' ,ylab='Silhouette Coefficient')



##########
UTARGET_LABEL=LABEL
UTARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='Oligdend'
UTARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='Oligdend'

USE=which(UTARGET_LABEL %in% c('Oligdend','interneurons_batch1','pyramidal SS_batch1'))
BATCH=UTARGET_LABEL[USE]

COM_UUMAP=immune.combined@reductions$umap@cell.embeddings[USE,]
library(vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3])
boxplot(out1[,3],ylim=c(-1,1), names=c('NEW'),las=2,main='Oligodend Merged' ,ylab='Silhouette Coefficient')

