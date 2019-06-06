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


#length(x = pbmc@var.genes)

mybeer5 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=5)
mybeer50 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=50)
mybeer100 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=100)
saveRDS(mybeer5, 'mybeer5.RDS')
saveRDS(mybeer50, 'mybeer50.RDS')
saveRDS(mybeer100, 'mybeer100.RDS')

mybeer3010 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=30)
mybeer3050 <- BEER(D1, D2, CNUM=50, PCNUM=50, CPU=2,PP=30)
mybeer30100 <- BEER(D1, D2, CNUM=100, PCNUM=50, CPU=2,PP=30)

saveRDS(mybeer3010, 'mybeer3010.RDS')
saveRDS(mybeer3050, 'mybeer3050.RDS')
saveRDS(mybeer30100, 'mybeer30100.RDS')


LABEL=c(paste0(as.character(ori_label[,2]),'_batch1'), paste0(EXP@meta.data$label,'_batch2') )
LABEL[which(LABEL %in% c('MOL1_batch2','MOL2_batch2','MOL3_batch2','MOL4_batch2','MOL5_batch2','MOL6_batch2'))]='Mature Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('MFOL1_batch2','MFOL2_batch2'))]='Myelin-forming Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('NFOL1_batch2','NFOL2_batch2'))]='Newly-formed Oligodendrocytes_batch2'
LABEL[which(LABEL %in% c('COP_batch2'))]='Differentiation-committed oligodendrocyte precursors_batch2'



mybeer=mybeer5

pbmc <- mybeer$seurat
pbmc@meta.data$label=LABEL
PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
pdf('mybeer5.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)
dev.off()


mybeer=mybeer50

pbmc <- mybeer$seurat
pbmc@meta.data$label=LABEL
PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
pdf('mybeer50.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)
dev.off()


mybeer=mybeer100

pbmc <- mybeer$seurat
pbmc@meta.data$label=LABEL
PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
pdf('mybeer100.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)
dev.off()

mybeer=mybeer3050

pbmc <- mybeer$seurat
pbmc@meta.data$label=LABEL
PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
pdf('mybeer3050.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)
dev.off()

mybeer=mybeer30100

pbmc <- mybeer$seurat
pbmc@meta.data$label=LABEL
PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
pdf('mybeer30100.pdf',width=15,height=12)
DimPlot(pbmc, reduction.use='umap', group.by='label', pt.size=0.1,do.label=T)
dev.off()





