
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

COM_UMAP=readRDS('COM_UMAP.RDS')
MNN_UMAP=readRDS('MNN_UMAP.RDS')
BEER_UMAP=readRDS('BEER_UMAP.RDS')
CCA_UMAP=readRDS('CCA_UMAP.RDS')
BBK_UMAP=readRDS('BBK_UMAP.RDS')

USE=which(TARGET_LABEL %in% c("D1","D2"))
BATCH=TARGET_LABEL[USE]

COM_KBET=COM_UMAP[USE,]
MNN_KBET=MNN_UMAP[USE,]
BEER_KBET=BEER_UMAP[USE,]
CCA_KBET=CCA_UMAP[USE,]
BBK_KBET=BBK_UMAP[USE,]

COMV <- kBET(COM_KBET, BATCH)
MNNV <- kBET(MNN_KBET, BATCH)
BEERV <- kBET(BEER_KBET, BATCH)
CCAV <- kBET(CCA_KBET, BATCH)
BBKV <- kBET(BBK_KBET, BATCH)



