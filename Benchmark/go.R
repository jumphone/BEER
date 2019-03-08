
library('Seurat')
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#####Load Zeisel#######
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
#23556  5069
#############

mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)

pbmc=mybeer$seurat


LABEL=c(paste0(as.character(ori_label[,2]),'_batch1'), paste0(EXP@meta.data$label,'_batch2') )
LABEL[which(LABEL %in% c('MOL1_batch2','MOL2_batch2','MOL3_batch2','MOL4_batch2','MOL5_batch2','MOL6_batch2'))]='Mature Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('MFOL1_batch2','MFOL2_batch2'))]='Myelin-forming Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('NFOL1_batch2','NFOL2_batch2'))]='Newly-formed Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('COP_batch2'))]='Differentiation-committed oligodendrocyte precursors_batch2'
pbmc@meta.data$label=LABEL



#NONE
ALLPC <- 1:length(mybeer$cor)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
NONE_DR=pbmc@dr$pca@cell.embeddings
NONE_UMAP=pbmc@dr$umap@cell.embeddings





#BEER
PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
BEER_DR=pbmc@dr$pca@cell.embeddings[,PCUSE]
BEER_UMAP=pbmc@dr$umap@cell.embeddings
saveRDS(pbmc,file='BEER.RDS')

#CCA
ctrl <- CreateSeuratObject(raw.data = D1, project = "D1", min.cells = 0)
ctrl@meta.data$stim <- "D1"
ctrl <- FilterCells(ctrl, subset.names = "nGene", low.thresholds = 0, high.thresholds = Inf)
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)
stim <- CreateSeuratObject(raw.data = D2, project = "D2", min.cells = 9)
stim@meta.data$stim <- "D2"
stim <- FilterCells(stim, subset.names = "nGene", low.thresholds = 0, high.thresholds = Inf)
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))
immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30)
immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim",  dims.align = 1:20)
immune.combined <- RunUMAP(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)

immune.combined@meta.data$label=LABEL
DimPlot(immune.combined, reduction.use='umap', group.by='stim', pt.size=0.1)
DimPlot(immune.combined, reduction.use='umap', group.by='label', pt.size=0.1)
saveRDS(immune.combined,file='CCA.RDS')

CCA_DR=immune.combined@dr$cca.aligned@cell.embeddings
CCA_UMAP=immune.combined@dr$umap@cell.embeddings



#MNN
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("scran", version = "3.8")
library(scran)


EXP=.simple_combine(D1,D2)

#CD1=EXP$exp_sc_mat1[which(rownames(CD1) %in% pbmc@var.genes),]
#CD2=EXP$exp_sc_mat2[which(rownames(CD2) %in% pbmc@var.genes),]
CD1=EXP$exp_sc_mat1[which(rownames(CD1) %in% pbmc@var.genes),]
CD2=EXP$exp_sc_mat2[which(rownames(CD2) %in% pbmc@var.genes),]

gene.counts1=CD1
sce1 <- SingleCellExperiment(list(counts=gene.counts1))
sce1 <- normalize(sce1)

gene.counts2=CD2
sce2 <- SingleCellExperiment(list(counts=gene.counts2))
sce2 <- normalize(sce2)

b1 <- sce1
b2 <- sce2
Sys.time()
out <- fastMNN(b1, b2)
Sys.time()


dim(out$corrected)

pbmc_mnn=pbmc
pbmc_mnn@dr$pca@cell.embeddings=as.matrix(out$corrected)

rownames(pbmc_mnn@dr$pca@cell.embeddings)=rownames(pbmc@dr$pca@cell.embeddings)
colnames(pbmc_mnn@dr$pca@cell.embeddings)=colnames(pbmc@dr$pca@cell.embeddings)

ALLPC <- 1:length(mybeer$cor)
pbmc_mnn <- RunUMAP(object = pbmc_mnn, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
pbmc_mnn@meta.data$label=LABEL
DimPlot(pbmc_mnn, reduction.use='umap', group.by='batch', pt.size=0.1)
DimPlot(pbmc_mnn, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)
MNN_DR=pbmc_mnn@dr$pca@cell.embeddings
MNN_UMAP=pbmc_mnn@dr$umap@cell.embeddings

DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)

##############



TARGET_LABEL=rep('NA',length(LABEL))
TARGET_LABEL[which(LABEL=='oligodendrocytes')]='D1'
TARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes','Newly-formed Oligodendrocytes','Mature Oligodendrocytes'))]='D2'


PCH=rep(20,length(TARGET_LABEL))
PCH[which(TARGET_LABEL=='D1')]=3
PCH[which(TARGET_LABEL=='D2')]=4

COL=rep('grey90',length(TARGET_LABEL))
COL[which(TARGET_LABEL=='D1')]='red'
COL[which(TARGET_LABEL=='D2')]='blue'


TOTAL=length(which(PCH==3))

CEX=0.4
par(mfrow=c(2,2))

plot(NONE_UMAP, col=COL,pch=PCH,cex=CEX, main='Original')
points(NONE_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)


plot(BEER_UMAP, col=COL,pch=PCH,cex=CEX, main='BEER')
points(BEER_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)


plot(CCA_UMAP, col=COL,pch=PCH,cex=CEX, main='Seurat (CCA alignment)')
points(CCA_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)

plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)









XL=0;XR=5;YB=0;YU=5
rect(XL,YB,XR,YU)
RNUM=length(which(NONE_UMAP[which(PCH==3),1]>XL & NONE_UMAP[which(PCH==3),1]<XR & NONE_UMAP[which(PCH==3),2]>YB & NONE_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL



par(mfrow=c(4,2))
plot(density(NONE_UMAP[which(PCH==3),1]),ylim=c(0,0.5),col='red')
lines(density(NONE_UMAP[which(PCH==4),1]),col='blue')
plot(density(NONE_UMAP[which(PCH==3),2]),ylim=c(0,0.5),col='red')
lines(density(NONE_UMAP[which(PCH==4),2]),col='blue')

plot(density(BEER_UMAP[which(PCH==3),1]),ylim=c(0,0.5),col='red')
lines(density(BEER_UMAP[which(PCH==4),1]),col='blue')
plot(density(BEER_UMAP[which(PCH==3),2]),ylim=c(0,0.5),col='red')
lines(density(BEER_UMAP[which(PCH==4),2]),col='blue')


plot(density(CCA_UMAP[which(PCH==3),1]),ylim=c(0,0.5),col='red')
lines(density(CCA_UMAP[which(PCH==4),1]),col='blue')
plot(density(CCA_UMAP[which(PCH==3),2]),ylim=c(0,0.5),col='red')
lines(density(CCA_UMAP[which(PCH==4),2]),col='blue')


plot(density(MNN_UMAP[which(PCH==3),1]),ylim=c(0,0.5),col='red')
lines(density(MNN_UMAP[which(PCH==4),1]),col='blue')
plot(density(MNN_UMAP[which(PCH==3),2]),ylim=c(0,0.5),col='red')
lines(density(MNN_UMAP[which(PCH==4),2]),col='blue')



#ks.test(MNN_UMAP[which(PCH==3),1],MNN_UMAP[which(PCH==4),1])



