
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

##########################################################

#library(devtools)
#install_github('theislab/kBET')
library(kBET)

COM_D=readRDS('COM_DR.RDS')
MNN_D=readRDS('MNN_DR.RDS')
BEER_D=readRDS('BEER_DR.RDS')
CCA_D=readRDS('CCA_DR.RDS')


USE=which(TARGET_LABEL %in% c("D1","D2"))
BATCH=TARGET_LABEL[USE]
SIZE=length(BATCH)

COM_KBET=COM_D[USE,]
MNN_KBET=MNN_D[USE,]
BEER_KBET=BEER_D[USE,]
CCA_KBET=CCA_D[USE,]


pca.data=a
pca.data$rotation=COM_KBET
pca.data$x=pca.data$rotation
pca.data$sdev=c(ncol(pca.data$rotation):1)

COMP <- batch.silhouette <- batch_sil(pca.data, BATCH, do.PCA=FALSE)

COMV <- kBET(COM_KBET, BATCH, do.pca=FALSE,testSize=SIZE) #mean kBET rejection rate: 0.9932873
MNNV <- kBET(MNN_KBET, BATCH, do.pca=FALSE,testSize=SIZE) #mean kBET rejection rate: 0.8719001
BEERV <- kBET(BEER_KBET, BATCH, do.pca=FALSE,testSize=SIZE) #mean kBET rejection rate: 0.9168376
CCAV <- kBET(CCA_KBET, BATCH, do.pca=FALSE,testSize=SIZE) #mean kBET rejection rate:0.9285847


