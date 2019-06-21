
pca.all=pbmc@reductions$pca@cell.embeddings
pca.use=pbmc@reductions$pca@cell.embeddings[,PCUSE]
batch=as.character(pbmc@meta.data$batch)

saveRDS(pca.all,file='pca.all.RDS')
saveRDS(pca.use,file='pca.use.RDS')
saveRDS(batch,file='batch.RDS')




pca.all=readRDS('pca.all.RDS')
pca.use=readRDS('pca.use.RDS')
batch=readRDS('batch.RDS')



library(reticulate)
#use_python("C:\Users\cchmc\Anaconda3\python")

anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy.api",convert=FALSE)

adata = anndata$AnnData(X=pca.all, obs=batch)
PCNUM=ncol(pca.use)

sc$tl$pca(adata, n_comps=as.integer(PCNUM))
adata$obsm$X_pca = pca.use
NB=50
NT=10

print(NB)
bbknn$bbknn(adata,batch_key=0,neighbors_within_batch=as.integer(NB),n_pcs=as.integer(PCNUM), n_trees =as.integer(NT))
sc$tl$umap(adata)
umap = py_to_r(adata$obsm$X_umap)


plot(umap)



saveRDS(umap,file='umap.RDS')



umap=readRDS('umap.RDS')
rownames(umap)=rownames(pbmc@reductions$umap@cell.embeddings)
colnames(umap)=colnames(pbmc@reductions$umap@cell.embeddings)


pbmc@reductions$umap@cell.embeddings=umap

