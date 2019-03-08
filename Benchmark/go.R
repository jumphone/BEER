
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

################################################################################################
Sys.time()
mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)
Sys.time()

#BEER Start Time: 2019-03-08 13:56:14 EST
#BEER End Time: 2019-03-08 14:01:39 EST
#BEER Time: 5 min 


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


################################################################################################

#CCA
Sys.time()
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
Sys.time()

#CCA Start Time: 2019-03-08 14:05:30 EST
#CCA End Time: 2019-03-08 14:09:30 EST
#CCA Time: 4 min 
immune.combined <- RunUMAP(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)

immune.combined@meta.data$label=LABEL
DimPlot(immune.combined, reduction.use='umap', group.by='stim', pt.size=0.1)
DimPlot(immune.combined, reduction.use='umap', group.by='label', pt.size=0.1)
#saveRDS(immune.combined,file='CCA.RDS')

CCA_DR=immune.combined@dr$cca.aligned@cell.embeddings
CCA_UMAP=immune.combined@dr$umap@cell.embeddings



################################################################################################
#MNN
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("scran", version = "3.8")
library(scran)


EXP=.simple_combine(D1,D2)

CD1=EXP$exp_sc_mat1 #[which(rownames(EXP$exp_sc_mat1) %in% pbmc@var.genes),]
CD2=EXP$exp_sc_mat2 #[which(rownames(EXP$exp_sc_mat2) %in% pbmc@var.genes),]

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

#MNN Start Time: 2019-03-08 14:18:02 EST
#MNN End Time: 2019-03-08 14:53:31 EST
#MNN Time: 35 min

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
################################################################################################


TARGET_LABEL=rep('NA',length(LABEL))
TARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='D1'
TARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='D2'


PCH=rep(20,length(TARGET_LABEL))
PCH[which(TARGET_LABEL=='D1')]=3
PCH[which(TARGET_LABEL=='D2')]=4

COL=rep('grey90',length(TARGET_LABEL))
COL[which(TARGET_LABEL=='D1')]='red'
COL[which(TARGET_LABEL=='D2')]='blue'


TOTAL=length(which(PCH==3)) #820
#length(which(PCH==4) #4543
LWD=1.2
CEX=0.05
par(mfrow=c(2,2))
###############
plot(NONE_UMAP, col=COL,pch=PCH,cex=CEX, main='None')
points(NONE_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=5;YB=0;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(NONE_UMAP[which(PCH==3),1]>XL & NONE_UMAP[which(PCH==3),1]<XR & NONE_UMAP[which(PCH==3),2]>YB & NONE_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4073171

###############
plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=-5;XR=0;YB=-8;YU=-2
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(MNN_UMAP[which(PCH==3),1]>XL & MNN_UMAP[which(PCH==3),1]<XR & MNN_UMAP[which(PCH==3),2]>YB & MNN_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4902439

###############
plot(CCA_UMAP, col=COL,pch=PCH,cex=CEX, main='Seurat (CCA alignment)')
points(CCA_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=7;YB=1;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(CCA_UMAP[which(PCH==3),1]>XL & CCA_UMAP[which(PCH==3),1]<XR & CCA_UMAP[which(PCH==3),2]>YB & CCA_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.3195122

###############
plot(BEER_UMAP, col=COL,pch=PCH,cex=CEX, main='BEER')
points(BEER_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
XL=-8;XR=0;YB=-9;YU=0
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(BEER_UMAP[which(PCH==3),1]>XL & BEER_UMAP[which(PCH==3),1]<XR & BEER_UMAP[which(PCH==3),2]>YB & BEER_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.5341463




#save.image('SAVE.RData')

saveRDS(NONE_DR,'NONE_DR.RDS')
saveRDS(NONE_UMAP,'NONE_UMAP')

saveRDS(MNN_DR,'MNN_DR')
saveRDS(MNN_UMAP,'MNN_UMAP')

saveRDS(BEER_DR,'BEER_DR')
saveRDS(BEER_UMAP,'BEER_UMAP')

saveRDS(CCA_DR,'CCA_DR')
saveRDS(CCA_UMAP,'CA_UMAP')




########################################################
DimPlot(pbmc_mnn, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)





TOTAL=length(which(PCH==3)) #820
#length(which(PCH==4) #4543
LWD=1.2
CEX=0.05
par(mfrow=c(2,3))
###############
plot(COM_UMAP, col=COL,pch=PCH,cex=CEX, main='Combat')
points(CO,_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=5;YB=0;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(NONE_UMAP[which(PCH==3),1]>XL & NONE_UMAP[which(PCH==3),1]<XR & NONE_UMAP[which(PCH==3),2]>YB & NONE_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4073171

###############
plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=-5;XR=0;YB=-8;YU=-2
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(MNN_UMAP[which(PCH==3),1]>XL & MNN_UMAP[which(PCH==3),1]<XR & MNN_UMAP[which(PCH==3),2]>YB & MNN_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4902439


###############
plot(BEER_UMAP, col=COL,pch=PCH,cex=CEX, main='BEER')
points(BEER_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
XL=-8;XR=0;YB=-9;YU=0
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(BEER_UMAP[which(PCH==3),1]>XL & BEER_UMAP[which(PCH==3),1]<XR & BEER_UMAP[which(PCH==3),2]>YB & BEER_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.5341463


###############
plot(CCA_UMAP, col=COL,pch=PCH,cex=CEX, main='Seurat (CCA alignment)')
points(CCA_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=7;YB=1;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(CCA_UMAP[which(PCH==3),1]>XL & CCA_UMAP[which(PCH==3),1]<XR & CCA_UMAP[which(PCH==3),2]>YB & CCA_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.3195122



NCOL=rep('grey90',nrow(BEER_DR))
NCOL[which(LABEL=='astrocytes_ependymal_batch1')]='red'
NCOL[which(LABEL=='OPC_batch2')]='blue'
NCOL[which(LABEL=='microglia_batch1')]='darkgreen'

NCEX=0.1
plot(MNN_UMAP, col=NCOL,pch=19,cex=NCEX, main='fastMNN')
#points(MNN_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(MNN_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)


plot(BEER_UMAP, col=NCOL,pch=19,cex=NCEX, main='BEER')
#points(BEER_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(BEER_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)


###############@@@@@@@@@@@@@@@@@@@


TOTAL=length(which(PCH==3)) #820
#length(which(PCH==4) #4543
LWD=1.2
CEX=0.05
par(mfrow=c(2,3))
###############
plot(NONE_UMAP, col=COL,pch=PCH,cex=CEX, main='None')
points(NONE_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=5;YB=0;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(NONE_UMAP[which(PCH==3),1]>XL & NONE_UMAP[which(PCH==3),1]<XR & NONE_UMAP[which(PCH==3),2]>YB & NONE_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4073171

###############
plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=-5;XR=0;YB=-8;YU=-2
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(MNN_UMAP[which(PCH==3),1]>XL & MNN_UMAP[which(PCH==3),1]<XR & MNN_UMAP[which(PCH==3),2]>YB & MNN_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4902439



NCOL=rep('grey90',nrow(BEER_DR))
NCOL[which(LABEL=='astrocytes_ependymal_batch1')]='red'
NCOL[which(LABEL=='OPC_batch2')]='blue'
NCOL[which(LABEL=='microglia_batch1')]='darkgreen'

NCEX=0.1
plot(MNN_UMAP, col=NCOL,pch=19,cex=NCEX, main='fastMNN')
#points(MNN_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(MNN_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)



###############
plot(CCA_UMAP, col=COL,pch=PCH,cex=CEX, main='Seurat (CCA alignment)')
points(CCA_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=7;YB=1;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(CCA_UMAP[which(PCH==3),1]>XL & CCA_UMAP[which(PCH==3),1]<XR & CCA_UMAP[which(PCH==3),2]>YB & CCA_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.3195122


###############
plot(BEER_UMAP, col=COL,pch=PCH,cex=CEX, main='BEER')
points(BEER_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
XL=-8;XR=0;YB=-9;YU=0
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(BEER_UMAP[which(PCH==3),1]>XL & BEER_UMAP[which(PCH==3),1]<XR & BEER_UMAP[which(PCH==3),2]>YB & BEER_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.5341463

plot(BEER_UMAP, col=NCOL,pch=19,cex=NCEX, main='BEER')
#points(BEER_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(BEER_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)





######Combat-sva############
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("sva", version = "3.8")

library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)


Sys.time()
pheno = data.frame(batch=as.matrix(pbmc@meta.data$batch))
edata = EXP$combine
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

pbmc_com=CreateSeuratObject(raw.data = combat_edata, project = "combat", min.cells = 0)
pbmc_com <- NormalizeData(object = pbmc_com, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_com <- FindVariableGenes(object = pbmc_com, mean.function = ExpMean, dispersion.function = LogVMR,do.plot=F, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc_com <- ScaleData(object = pbmc_com)
PCNUM=50
pbmc_com <- RunPCA(object = pbmc_com, pc.genes = pbmc@var.genes, pcs.compute=PCNUM,do.print = TRUE, pcs.print = 1:5, genes.print = 5)
Sys.time()

#Combat start: 2019-03-08 17:46:21 EST
#Combat End: 2019-03-08 17:48:25 EST
#Combat time:2

ALLPC=1:50
pbmc_com <- RunUMAP(object = pbmc_com, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)

COM_DR=pbmc_com@dr$pca@cell.embeddings
COM_UMAP=pbmc_com@dr$umap@cell.embeddings


###############



TOTAL=length(which(PCH==3)) #820
#length(which(PCH==4) #4543
LWD=1.2
CEX=0.05
par(mfrow=c(2,3))
###############
plot(COM_UMAP, col=COL,pch=PCH,cex=CEX, main='Combat')
points(COM,_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=-3;XR=3;YB=3;YU=7
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(COM_UMAP[which(PCH==3),1]>XL & COM_UMAP[which(PCH==3),1]<XR & COM_UMAP[which(PCH==3),2]>YB & COM_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.7682927



###############
plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=-5;XR=0;YB=-8;YU=-2
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(MNN_UMAP[which(PCH==3),1]>XL & MNN_UMAP[which(PCH==3),1]<XR & MNN_UMAP[which(PCH==3),2]>YB & MNN_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.4902439


###############
NCOL=rep('grey90',nrow(BEER_DR))
NCOL[which(LABEL=='astrocytes_ependymal_batch1')]='red'
NCOL[which(LABEL=='OPC_batch2')]='blue'
NCOL[which(LABEL=='microglia_batch1')]='darkgreen'

NCEX=0.1
plot(MNN_UMAP, col=NCOL,pch=19,cex=NCEX, main='fastMNN')
#points(MNN_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(MNN_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)



###############
plot(CCA_UMAP, col=COL,pch=PCH,cex=CEX, main='Seurat (CCA alignment)')
points(CCA_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=7;YB=1;YU=5
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(CCA_UMAP[which(PCH==3),1]>XL & CCA_UMAP[which(PCH==3),1]<XR & CCA_UMAP[which(PCH==3),2]>YB & CCA_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.3195122


###############
plot(BEER_UMAP, col=COL,pch=PCH,cex=CEX, main='BEER')
points(BEER_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
XL=-8;XR=0;YB=-9;YU=0
rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
RNUM=length(which(BEER_UMAP[which(PCH==3),1]>XL & BEER_UMAP[which(PCH==3),1]<XR & BEER_UMAP[which(PCH==3),2]>YB & BEER_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL #0.5341463

plot(BEER_UMAP, col=NCOL,pch=19,cex=NCEX, main='BEER')
#points(BEER_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(BEER_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)

























###############



XL=-5;XR=0;YB=0;YU=5
rect(XL,YB,XR,YU)
RNUM=length(which(NONE_UMAP[which(PCH==3),1]>XL & NONE_UMAP[which(PCH==3),1]<XR & NONE_UMAP[which(PCH==3),2]>YB & NONE_UMAP[which(PCH==3),2]<YU))
RNUM/TOTAL
NCOL=rep()


plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)

plot(BEER_UMAP, col=NCOL,pch=NPCH,cex=CEX, main='BEER')
points(BEER_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)









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



