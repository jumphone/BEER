<img src="https://github.com/jumphone/BEER/raw/master/DATA/BEER_LOGO.png" width="200">

### BEER: Batch EffEct Remover for single-cell data

Author: Feng Zhang

Date: Mar. 7, 2019


# Requirement:

    #R >=3.5
    install.packages('Seurat')
    install.packages('pcaPP')

# Usage:

* [I. Combine Two Batches](#I-Combine-Two-Batches)
* [II. Combine Multiple Batches](#II-Combine-Multiple-Batches)
* [III. UMAP-based Clustering](#III-UMAP-based-Clustering)

</br>
</br>

# I. Combine Two Batches

Please use the function named "BEER" to combine two batches.

Download demo data: https://github.com/jumphone/BEER/raw/master/DATA/demodata.zip 

Please do basic quality control before using BEER (e.g. remove low-quality cells & genes).

### Step1. Load Data

    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    
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


### Step2. Detect Batch Effect

    mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)
    
    #CNUM: the number of cells in each group
    #PCNUM: the number of computated PCA subspaces 
    
    #If the gene number per cell is small, please try a larger "CNUM", e.g. 20.
    
    par(mfrow=c(1,2))
    plot(mybeer$cor, xlab='PCs', ylab="COR", pch=16)
    plot(-log(mybeer$fdr,10), xlab='PCs', ylab='-log10(FDR)', pch=16)
    
### Step3. Visualization 
    
#### Keep batch effect:
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/KeepBatchEffect.png" width="400">
    
    ALLPC <- 1:length(mybeer$cor)
    pbmc <- mybeer$seurat
    
    # UMAP:
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
    #DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
    
    # tSNE:
    pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = ALLPC, do.fast = TRUE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='tsne', group.by='batch', pt.size=0.1)
    #DimPlot(pbmc, reduction.use='tsne', group.by='map', pt.size=0.1)

#### Remove batch effect:

<img src="https://github.com/jumphone/BEER/raw/master/DATA/RemoveBatchEffect.png" width="400">
  
    PCUSE <- which(mybeer$cor> min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05)
    # Users can set the cutoff of "mybeer$cor" based on the distribution of "mybeer$cor".
    pbmc <- mybeer$seurat
    
    # UMAP:
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
    #DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
    
    # tSNE:
    pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='tsne', group.by='batch', pt.size=0.1)
    #DimPlot(pbmc, reduction.use='tsne', group.by='map', pt.size=0.1)
    
    
</br>
</br>
    
# II. Combine Multiple Batches

Please use the function named "MBEER" to combine multiple batches (n>=3).

This function implements the iteration of "Combine Two Batches".

MBEER compares each batch with the batch having the largest cell number.

The assumption is that the batch having the largest cell number has almost all cell-types within all batches.

If you want to define a batch having almost all cell-types, please set "MAXBATCH" to the label of that batch.

Download demo data: https://sourceforge.net/projects/beergithub/files/
   
### Step1. Load Data
    
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    
    #Load Demo Data (Oligodendroglioma, GSE70630)
    #Download: https://sourceforge.net/projects/beergithub/files/
    
    D1=readRDS('MGH36.RDS')
    D2=readRDS('MGH53.RDS')
    D3=readRDS('MGH54.RDS')
    D4=readRDS('MGH60.RDS')
    D5=readRDS('MGH93.RDS')
    D6=readRDS('MGH97.RDS')
    
    BATCH=c(rep('D1',ncol(D1)),
            rep('D2',ncol(D2)),
            rep('D3',ncol(D3)),
            rep('D4',ncol(D4)),
            rep('D5',ncol(D5)),
            rep('D6',ncol(D6)) )
            
    D12=.simple_combine(D1,D2)$combine
    D34=.simple_combine(D3,D4)$combine
    D56=.simple_combine(D5,D6)$combine
    D1234=.simple_combine(D12,D34)$combine
    D123456=.simple_combine(D1234,D56)$combine
    
    DATA=D123456   
   
    rm(D1);rm(D2);rm(D3);rm(D4);rm(D5);rm(D6)
    rm(D12);rm(D34);rm(D56);rm(D1234);rm(D123456)
    
### Step2. Detect Batch Effect

    mybeer=MBEER(DATA, BATCH, MAXBATCH="", CNUM=10, PCNUM=20, CPU=2, SEED=1 )

    par(mfrow=c(1,2))
    plot(mybeer$cor, xlab='PCs', ylab='COR', pch=16)
    plot(-log(mybeer$fdr,10), xlab='PCs', ylab='-log10(FDR)', pch=16)
    
### Step3. Visualization 
    
#### Keep batch effect:
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/MBEER1.png" width="400">
    
    ALLPC <- 1:length(mybeer$cor)
    pbmc <- mybeer$seurat
    
    # UMAP:
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
   

#### Remove batch effect:

<img src="https://github.com/jumphone/BEER/raw/master/DATA/MBEER2.png" width="400">

    PCUSE <- which(mybeer$cor> 0.7  & mybeer$fdr<0.05)
    # Users can set the cutoff of "mybeer$cor" based on the distribution of "mybeer$cor".
    pbmc <- mybeer$seurat
    
    # UMAP:
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)  
    
</br>   

### *Another Demo of MBEER (GSE102130):

Here, we show the final UMAP figures (All parameters are the same with that of the first demo).

#### Keep batch effect:
<img src="https://github.com/jumphone/BEER/raw/master/DATA/MBEER3.png" width="400">

#### Remove batch effect:
<img src="https://github.com/jumphone/BEER/raw/master/DATA/MBEER4.png" width="400">

</br>
</br>

# III. UMAP-based Clustering
    
    #Demo Data (GSE102130)
    VEC=pbmc@dr$umap@cell.embeddings

### Step1. Clustering:

<img src="https://github.com/jumphone/BEER/raw/master/DATA/CLUST1.png" width="400">    
    
    # Here, we use the "dbscan" function to do clustering.
    library("fpc")
    set.seed(123)
    df=VEC
    db <- fpc::dbscan(df, eps = 0.5, MinPts = 5)
    DC=db$cluster
    pbmc@meta.data$DC=DC
    DimPlot(pbmc, reduction.use='umap', group.by='DC', pt.size=0.5)
    
    
    
### Step2. Find marker genes & draw heatmap:

* Details are in the instruction page of Seurat: https://satijalab.org/seurat/get_started.html

<img src="https://github.com/jumphone/BEER/raw/master/DATA/CLUST2.png" width="400">    

    library(dplyr)
    tmp=pbmc@ident
    pbmc@ident=as.factor(DC)
    names(pbmc@ident)=names(tmp)
    pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
    top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
    DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
    
    
    
</br>   
</br> 
    
    
# License
    
    MIT License
    
    Copyright (c) 2019 Zhang, Feng

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

  
