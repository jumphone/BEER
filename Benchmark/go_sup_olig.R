MGH36=read.table('MGH36.txt',header=T,row.names=1,sep='\t')
MGH53=read.table('MGH53.txt',header=T,row.names=1,sep='\t')
MGH54=read.table('MGH54.txt',header=T,row.names=1,sep='\t')
MGH60=read.table('MGH60.txt',header=T,row.names=1,sep='\t')
MGH93=read.table('MGH93.txt',header=T,row.names=1,sep='\t')
MGH97=read.table('MGH97.txt',header=T,row.names=1,sep='\t')


library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

############################
D1 <- MGH36
D2 <- MGH53
mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)
pbmc <- mybeer$seurat


ALLPC <- 1:length(mybeer$cor)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)


PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)

##############################



D1 <- MGH54
D2 <- MGH60
mybeer <- BEER(D1, D2, CNUM=5, PCNUM=50, CPU=2)

pbmc <- mybeer$seurat


ALLPC <- 1:length(mybeer$cor)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)


PCUSE <- which(mybeer$cor> 0.8  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)

##############################

D1 <- MGH93
D2 <- MGH97

mybeer <- BEER(D1, D2, CNUM=2, PCNUM=50, CPU=2)

pbmc <- mybeer$seurat

ALLPC <- 1:length(mybeer$cor)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)


PCUSE <- which(mybeer$cor> 0.8  & mybeer$fdr<0.05)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)


##############################





