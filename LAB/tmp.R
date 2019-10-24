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



DATA=DATA
BATCH=BATCH
VEC=pbmc@reductions$umap@cell.embeddings
















