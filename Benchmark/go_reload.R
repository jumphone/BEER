
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

#####################################


LABEL=c(paste0(as.character(ori_label[,2]),'_batch1'), paste0(EXP@meta.data$label,'_batch2') )
LABEL[which(LABEL %in% c('MOL1_batch2','MOL2_batch2','MOL3_batch2','MOL4_batch2','MOL5_batch2','MOL6_batch2'))]='Mature Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('MFOL1_batch2','MFOL2_batch2'))]='Myelin-forming Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('NFOL1_batch2','NFOL2_batch2'))]='Newly-formed Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('COP_batch2'))]='Differentiation-committed oligodendrocyte precursors_batch2'



#mybeer <- BEER(D1, D2, CNUM=100, PCNUM=50, CPU=2)
pbmc=mybeer$seurat


####################################


COM_DR=readRDS('COM_DR.RDS')
COM_UMAP=readRDS('COM_UMAP.RDS')

NONE_DR=readRDS('NONE_DR.RDS')
NONE_UMAP=readRDS('NONE_UMAP.RDS')

MNN_DR=readRDS('MNN_DR.RDS')
MNN_UMAP=readRDS('MNN_UMAP.RDS')

BEER_DR=readRDS('BEER_DR.RDS')
BEER_UMAP=readRDS('BEER_UMAP.RDS')

CCA_DR=readRDS('CCA_DR.RDS')
CCA_UMAP=readRDS('CCA_UMAP.RDS')

BBK_UMAP=readRDS('BBK_UMAP.RDS')


pbmc@dr$umap=pbmc@dr$pca


ALLPC <- 1:10
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)

pbmc@meta.data$label=LABEL

BBK_UMAP=as.matrix(BBK_UMAP)
rownames(BBK_UMAP)=rownames(pbmc@dr$umap@cell.embeddings)
colnames(BBK_UMAP)=colnames(pbmc@dr$umap@cell.embeddings)
BBK_UMAP[which(BBK_UMAP[,1] < -8),1]= -8

pbmc@dr$umap@cell.embeddings=BBK_UMAP

DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)

########################################################


TARGET_LABEL=rep('NA',length(LABEL))
TARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='D1'
TARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='D2'


PCH=rep(20,length(TARGET_LABEL))
PCH[which(TARGET_LABEL=='D1')]=3
PCH[which(TARGET_LABEL=='D2')]=4

COL=rep('grey90',length(TARGET_LABEL))
COL[which(TARGET_LABEL=='D1')]='red'
COL[which(TARGET_LABEL=='D2')]='blue'






tiff("COMPARE.tif", width = 8, height = 8, units = 'in',res = 500)

TOTAL=length(which(PCH==3)) #820
#length(which(PCH==4) #4543
LWD=1.2
CEX=0.3
par(mfrow=c(3,3))
###############
plot(COM_UMAP, col=COL,pch=PCH,cex=CEX, main='Combat')
points(COM_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
#XL=-3;XR=3;YB=3;YU=7
#rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
#RNUM=length(which(COM_UMAP[which(PCH==3),1]>XL & COM_UMAP[which(PCH==3),1]<XR & COM_UMAP[which(PCH==3),2]>YB & COM_UMAP[which(PCH==3),2]<YU))
#RNUM/TOTAL #0.7682927



###############
plot(MNN_UMAP, col=COL,pch=PCH,cex=CEX, main='fastMNN')
points(MNN_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
#XL=-5;XR=0;YB=-8;YU=-2
#rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
#RNUM=length(which(MNN_UMAP[which(PCH==3),1]>XL & MNN_UMAP[which(PCH==3),1]<XR & MNN_UMAP[which(PCH==3),2]>YB & MNN_UMAP[which(PCH==3),2]<YU))
#RNUM/TOTAL #0.4902439


###############
NCOL=rep('grey90',nrow(BEER_DR))
NCOL[which(LABEL=='astrocytes_ependymal_batch1')]='red'
NCOL[which(LABEL=='OPC_batch2')]='blue'
NCOL[which(LABEL=='microglia_batch1')]='darkgreen'

NCEX=0.3
plot(MNN_UMAP, col=NCOL,pch=19,cex=NCEX, main='fastMNN')
#points(MNN_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(MNN_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)



###############
plot(CCA_UMAP, col=COL,pch=PCH,cex=CEX, main='Seurat (CCA alignment)')
points(CCA_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)
###############
XL=0;XR=7;YB=1;YU=5
#rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
#RNUM=length(which(CCA_UMAP[which(PCH==3),1]>XL & CCA_UMAP[which(PCH==3),1]<XR & CCA_UMAP[which(PCH==3),2]>YB & CCA_UMAP[which(PCH==3),2]<YU))
#RNUM/TOTAL #0.3195122


###############
plot(BEER_UMAP, col=COL,pch=PCH,cex=CEX, main='BEER')
used=which(PCH==3 & BEER_UMAP[,1]> -8 & BEER_UMAP[,1]< 0 & BEER_UMAP[,2]< 0 & BEER_UMAP[,2] > -9)
points(BEER_UMAP[used,], col=COL[used],pch=PCH[used],cex=CEX)
#XL=-8;XR=0;YB=-9;YU=0
#rect(XL,YB,XR,YU,border='black',lwd=LWD,lty='longdash')
#RNUM=length(which(BEER_UMAP[which(PCH==3),1]>XL & BEER_UMAP[which(PCH==3),1]<XR & BEER_UMAP[which(PCH==3),2]>YB & BEER_UMAP[which(PCH==3),2]<YU))
#RNUM/TOTAL #0.5341463

plot(BEER_UMAP, col=NCOL,pch=19,cex=NCEX, main='BEER')
#points(BEER_UMAP[which(NCOL=='red'),],pch=20, col=NCOL[which(NCOL=='red')],cex=CEX)
#points(BEER_UMAP[which(NCOL=='darkgreen'),],pch=20, col=NCOL[which(NCOL=='darkgreen')],cex=CEX)


plot(BBK_UMAP, col=COL,pch=PCH,cex=CEX, main='BBKNN')
points(BBK_UMAP[which(PCH==3),], col=COL[which(PCH==3)],pch=PCH[which(PCH==3)],cex=CEX)


NNCOL=rep('grey90',nrow(BEER_DR))
NNCOL[which(LABEL=='interneurons_batch1')]='red'
NNCOL[which(LABEL=='pyramidal SS_batch1')]='darkgreen'
NNCOL[which(TARGET_LABEL=='D1')]='blue'
NNCOL[which(TARGET_LABEL=='D2')]='blue'


plot(BBK_UMAP, col=NNCOL,pch=19,cex=0.3, main='BBKNN')
points(BBK_UMAP[which(NNCOL=='red'),], col=NNCOL[which(NNCOL=='red')],pch=19,cex=0.3)

plot(BEER_UMAP, col=NNCOL,pch=19,cex=0.3, main='BEER')
points(BEER_UMAP[which(NNCOL=='red'),], col=NNCOL[which(NNCOL=='red')],pch=19,cex=0.3)


dev.off()













