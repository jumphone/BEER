
NONE_DR=readRDS('NONE_DR.RDS')
NONE_UMAP=readRDS('NONE_UMAP.RDS')
pca=NONE_DR
batch=as.matrix(c(rep('D1',3005) ,rep('D2',5069)))


library(reticulate)
use_python("/usr/bin/python3")

anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy.api",convert=FALSE)

adata = anndata$AnnData(X=pca, obs=batch)
sc$tl$pca(adata)
adata$obsm$X_pca = pca
bbknn$bbknn(adata,batch_key=0)
sc$tl$umap(adata)
umap = py_to_r(adata$obsm$X_umap)

BBK_UMAP=umap
saveRDS(BBK_UMAP,'BBK_UMAP.RDS')

#################

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




LABEL=c(paste0(as.character(ori_label[,2]),'_batch1'), paste0(EXP@meta.data$label,'_batch2') )
LABEL[which(LABEL %in% c('MOL1_batch2','MOL2_batch2','MOL3_batch2','MOL4_batch2','MOL5_batch2','MOL6_batch2'))]='Mature Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('MFOL1_batch2','MFOL2_batch2'))]='Myelin-forming Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('NFOL1_batch2','NFOL2_batch2'))]='Newly-formed Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('COP_batch2'))]='Differentiation-committed oligodendrocyte precursors_batch2'






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
CEX=0.3
#par(mfrow=c(2,3))
###############
plot(BBK_UMAP, col=COL,pch=PCH,cex=CEX, main='BBKNN', xlab='UMAP1',ylab='UMAP2')
points(BBK_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
#XL=-3;XR=3;YB=3;YU=7
#rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
#RNUM=length(which(COM_UMAP[which(PCH==3),1]>XL & COM_UMAP[which(PCH==3),1]<XR & COM_UMAP[which(PCH==3),2]>YB & COM_UMAP[which(PCH==3),2]<YU))
#RNUM/TOTAL #0.7682927


###############
NCOL=rep('grey90',nrow(NONE_DR))
NCOL[which(LABEL=='astrocytes_ependymal_batch1')]='red'
NCOL[which(LABEL=='OPC_batch2')]='blue'
NCOL[which(LABEL=='microglia_batch1')]='darkgreen'

NCEX=0.3
plot(BBK_UMAP, col=NCOL,pch=19,cex=NCEX, main='fastMNN', xlab='UMAP1',ylab='UMAP2')










