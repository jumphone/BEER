
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
#############

mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)

pbmc=mybeer$seurat


LABEL=c(as.character(ori_label[,2]), EXP@meta.data$label )
LABEL[which(LABEL %in% c('MOL1','MOL2','MOL3','MOL4','MOL5','MOL6'))]='Mature Oligodendrocytes'
LABEL[which(LABEL %in% c('MFOL1','MFOL2'))]='Myelin-forming Oligodendrocytes'
LABEL[which(LABEL %in% c('NFOL1','NFOL2'))]='Newly-formed Oligodendrocytes'
LABEL[which(LABEL %in% c('COP'))]='Differentiation-committed oligodendrocyte precursors'
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
saveRDS(pbmc,file='BEER.RDS')

#CCA
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
immune.combined <- RunUMAP(immune.combined, reduction.use = "cca.aligned", dims.use = 1:20, do.fast = T)

immune.combined@meta.data$label=LABEL
DimPlot(immune.combined, reduction.use='umap', group.by='stim', pt.size=0.1)
DimPlot(immune.combined, reduction.use='umap', group.by='label', pt.size=0.1)
saveRDS(immune.combined,file='CCA.RDS')

CCA_DR=immune.combined@dr$cca.aligned@cell.embeddings
CCA_UMAP=immune.combined@dr$umap@cell.embeddings







