


library(reticulate)
use_python("/usr/bin/python3")

anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy.api",convert=FALSE)

adata = anndata$AnnData(X=pca, obs=batch)
sc$tl$pca(adata)
adata$obsm$X_pca = pca
bbknn$bbknn(adata,batch_key=0)
sc$tl$umap(adata)
umap = py_to_r(adata$obsm$X_umap)
