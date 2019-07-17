#####################
# Author: Feng Zhang
#####################
# Please adjust following parameters

PCA='used.pca.txt'
BATCH='batch.txt'
OUTPUT='bbknn_umap.txt'
NB=3
NT=10


#####################
#####################
#####################
print('PCA=',PCA)
print('BATCH=',BATCH)
print('OUTPUT=',OUTPUT)
print('NB=',NB)
print('NT=',NT)
#####################
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import bbknn
#####################
print('Start')
fi=open(BATCH)
batch=[]
for line in fi:
    seq=line.rstrip().split(',')
    batch=batch+seq
fi.close
batch=batch[1:]
used_pca=sc.read_csv(PCA) 
adata=anndata.AnnData(X=used_pca.X, obs=batch)
PCNUM=used_pca.X.shape[1] 
sc.tl.pca(adata, n_comps=PCNUM)
adata.obsm.X_pca = used_pca.X
bbknn.bbknn(adata,batch_key=0,neighbors_within_batch=NB,n_pcs=PCNUM, n_trees =NT)
sc.tl.umap(adata)
umap = adata.obsm.X_umap

fo=open(OUTPUT,'w')
for one in umap:
    fo.write(str(one[0])+'\t'+str(one[1])+'\n')
fo.close()
print('Finished !')
#####################
