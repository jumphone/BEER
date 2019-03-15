
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



COM_D=readRDS('COM_DR.RDS')
MNN_D=readRDS('MNN_DR.RDS')
BEER_D=readRDS('BEER_DR.RDS')
CCA_D=readRDS('CCA_DR.RDS')


COM_UMAP=readRDS('COM_UMAP.RDS')
NONE_UMAP=readRDS('NONE_UMAP.RDS')
MNN_UMAP=readRDS('MNN_UMAP.RDS')
BEER_UMAP=readRDS('BEER_UMAP.RDS')
CCA_UMAP=readRDS('CCA_UMAP.RDS')
BBK_UMAP=readRDS('BBK_UMAP.RDS')


OUT=c()

####################################################################################################




par(mfrow=c(1,3))


UTARGET_LABEL=LABEL
UTARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='Oligdend'
UTARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='Oligdend'

USE=which(!UTARGET_LABEL %in% c('NA'))
BATCH=UTARGET_LABEL[USE]

COM_UUMAP=COM_UMAP[USE,]
MNN_UUMAP=MNN_UMAP[USE,]
BEER_UUMAP=BEER_UMAP[USE,]
CCA_UUMAP=CCA_UMAP[USE,]
BBK_UUMAP=BBK_UMAP[USE,]

library(vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3])

dis = dist(CCA_UUMAP)
out2=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out2[,3]) 

dis = dist(BBK_UUMAP)
out3=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out3[,3]) 

dis = dist(MNN_UUMAP)
out4=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out4[,3]) 

dis = dist(BEER_UUMAP)
out5=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out5[,3])  

boxplot(out1[,3],out2[,3],out3[,3],out4[,3],out5[,3],ylim=c(-1,1), names=c('Combat','CCA','BBKAA','fastMNN','BEER'),las=2,main='Oligodend Merged' ,ylab='Silhouette Coefficient')
####################################################################################################
####################################################################################################

UTARGET_LABEL=LABEL
UTARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='Oligdend'
UTARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='Oligdend'

USE=which(UTARGET_LABEL %in% c('astrocytes_ependymal_batch1','OPC_batch2','microglia_batch1'))
BATCH=UTARGET_LABEL[USE]

COM_UUMAP=COM_UMAP[USE,]
MNN_UUMAP=MNN_UMAP[USE,]
BEER_UUMAP=BEER_UMAP[USE,]
CCA_UUMAP=CCA_UMAP[USE,]
BBK_UUMAP=BBK_UMAP[USE,]

library (vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3]) 
dis = dist(MNN_UUMAP)
out2=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out2[,3]) 
dis = dist(BEER_UUMAP)
out3=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out3[,3])  
dis = dist(CCA_UUMAP)
out4=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out4[,3])
dis = dist(BBK_UUMAP)
out5=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out5[,3]) 

library(vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3])

dis = dist(CCA_UUMAP)
out2=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out2[,3]) 

dis = dist(BBK_UUMAP)
out3=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out3[,3]) 

dis = dist(MNN_UUMAP)
out4=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out4[,3]) 

dis = dist(BEER_UUMAP)
out5=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out5[,3])  

boxplot(out1[,3],out2[,3],out3[,3],out4[,3],out5[,3],ylim=c(-1,1), names=c('Combat','CCA','BBKAA','fastMNN','BEER'),las=2,main='Astro & OPC & Microglia' ,ylab='Silhouette Coefficient')
####################################################################################################
####################################################################################################




UTARGET_LABEL=LABEL
UTARGET_LABEL[which(LABEL=='oligodendrocytes_batch1')]='Oligdend'
UTARGET_LABEL[which(LABEL %in% c('Myelin-forming Oligodendrocytes_batch2','Newly-formed Oligodendrocytes_batch2','Mature Oligodendrocytes_batch2'))]='Oligdend'

USE=which(UTARGET_LABEL %in% c('Oligdend','interneurons_batch1','pyramidal SS_batch1'))
BATCH=UTARGET_LABEL[USE]

COM_UUMAP=COM_UMAP[USE,]
MNN_UUMAP=MNN_UMAP[USE,]
BEER_UUMAP=BEER_UMAP[USE,]
CCA_UUMAP=CCA_UMAP[USE,]
BBK_UUMAP=BBK_UMAP[USE,]

library (vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3]) 
dis = dist(MNN_UUMAP)
out2=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out2[,3]) 
dis = dist(BEER_UUMAP)
out3=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out3[,3])  
dis = dist(CCA_UUMAP)
out4=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out4[,3]) 
dis = dist(BBK_UUMAP)
out5=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out5[,3])  

library(vegan)
library(cluster)

dis = dist(COM_UUMAP)
out1=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out1[,3])

dis = dist(CCA_UUMAP)
out2=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out2[,3]) 

dis = dist(BBK_UUMAP)
out3=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out3[,3]) 

dis = dist(MNN_UUMAP)
out4=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out4[,3]) 

dis = dist(BEER_UUMAP)
out5=silhouette(as.numeric(as.factor(BATCH)),dis)
summary(out5[,3])  

boxplot(out1[,3],out2[,3],out3[,3],out4[,3],out5[,3], ylim=c(-1,1), names=c('Combat','CCA','BBKAA','fastMNN','BEER'),las=2,main='Oligodend & Interneuron & Pyramidal_SS' ,ylab='Silhouette Coefficient')



#####################################












