setwd('C:/Users/cchmc/Desktop/BEER_IMP/')

library(Seurat)

source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
#source('BEER.R')

#Read 10X data: pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")

#Load Demo Data (subset of GSE70630: MGH53 & MGH54)
#Download: https://github.com/jumphone/BEER/raw/master/DATA/demodata.zip

D1 <- read.table(unz("demodata.zip","DATA1_MAT.txt"), sep='\t', row.names=1, header=T)
D2 <- read.table(unz("demodata.zip","DATA2_MAT.txt"), sep='\t', row.names=1, header=T)

# "D1" & "D2" are UMI matrix (or FPKM, RPKM, TPM, PKM ...; Should not be gene-centric scaled data)
# Rownames of "D1" & "D2" are gene names
# Colnames of "D1" & "D2" are cell names 

# There shouldn't be duplicated colnames in "D1" & "D2":
colnames(D1)=paste0('D1_', colnames(D1))
colnames(D2)=paste0('D2_', colnames(D2))

DATA=.simple_combine(D1,D2)$combine
BATCH=rep('D2',ncol(DATA))
BATCH[c(1:ncol(D1))]='D1'


mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   

pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1) 




#####################################################################
DATA=as.matrix(pbmc@assays$RNA@data)
BATCH=pbmc@meta.data$batch
VEC=pbmc@reductions$umap@cell.embeddings
print_step=100
CUTOFF=1

DATA=as.matrix(DATA)
NC=ncol(DATA)

DIST=dist(VEC)
DIST=as.matrix(DIST)

DIST.NUM=as.numeric(DIST)
E.DIST.NUM=ecdf(DIST.NUM)
PV=apply(DIST, 2, E.DIST.NUM)
LOG.PV=-log(PV,2)
LOG.PV[which(PV>CUTOFF)]=0



NEW.DATA=matrix(0,ncol=ncol(DATA),nrow=nrow(DATA))
rownames(NEW.DATA)=rownames(DATA)
colnames(NEW.DATA)=colnames(DATA)



i=1
while(i<=NC){
    NEW.DATA[,i]= DATA %*% as.matrix(LOG.PV[,i] ,ncol=1)
    NEW.DATA[,i] = NEW.DATA[,i] / sum(LOG.PV[,i])
    if(i %% print_step==1){print(i)}
    i=i+1}


NEW.DATA
##############################################################################################



source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

NEW.DATA=BEER.IMP(as.matrix(pbmc@assays$RNA@data),pbmc@reductions$umap@cell.embeddings)





pbmc_batch=CreateSeuratObject(counts = NEW.DATA, min.cells = 0, min.features = 0, project = "ALL") 
pbmc_batch@meta.data$batch=BATCH
pbmc_batch=FindVariableFeatures(object = pbmc_batch, selection.method = "vst", nfeatures = 2000)   
VariableFeatures(object = pbmc_batch)
pbmc_batch <- NormalizeData(object = pbmc_batch, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_batch <- ScaleData(object = pbmc_batch, features = VariableFeatures(object = pbmc_batch))
pbmc_batch <- RunPCA(object = pbmc_batch, seed.use=123, npcs=50, features = VariableFeatures(object = pbmc_batch), ndims.print=1,nfeatures.print=1)


pbmc_batch <- RunUMAP(pbmc_batch, dims = 1:50, seed.use = 123,n.components=2)
DimPlot(pbmc_batch, reduction.use='umap', group.by='batch', pt.size=1) 


VariableFeatures(pbmc_batch)
FeaturePlot(pbmc_batch, features=c('AQP4','SPP1','SOX2'))




