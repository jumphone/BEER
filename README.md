<img src="https://github.com/jumphone/BEER/raw/master/DATA/BEER_LOGO.png" width="200">

### BEER: Batch EffEct Remover for single-cell data

Author: Feng Zhang

* BEER's latest version: https://github.com/jumphone/BEER/releases

* BEER's manuscript version: https://github.com/jumphone/BEER/archive/0.0.2.zip

* Test version: https://github.com/jumphone/BEER/tree/master/LAB


# News:

* June 2019 ( v0.1 ): Batch-effect removal enhancement with ComBat & BBKNN

* ...

* June 2019 ( v0.0.4 ): "MBEER" is integrated into "BEER". Please directly use BEER to integrate multiple batches

* June 2019: stop updating the BEER source code for Seurat2. New feature is only for Seurat3.

# Requirement:

    #R >=3.5
    install.packages('Seurat') # >=3.0     

For batch-effect removal enhancement, please install Combat & BBKNN:

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("sva")
    BiocManager::install("limma")

    Install bbknn in python: https://github.com/Teichlab/bbknn

# Usage:

* [I. Combine Two Batches](#I-Combine-Two-Batches)
* [II. Combine Multiple Batches](#II-Combine-Multiple-Batches)
* [III. UMAP-based Clustering](#III-UMAP-based-Clustering)
* [IV. Combine scATAC-seq & scRNA-seq](#iv-combine-scatac-seq--scrna-seq)
* [V. Batch-effect Removal Enhancement](#v-batch-effect-removal-enhancement)
* [VI. Transfer Labels](#vi-transfer-labels)
* [VII. Biological Interpretation](#vii-biological-interpretation)

</br>

# I. Combine Two Batches

Download demo data: https://github.com/jumphone/BEER/raw/master/DATA/demodata.zip 

Please do basic quality control before using BEER (e.g. remove low-quality cells & genes). 

For QC, please see: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html 

### Step1. Load Data

    library(Seurat)
  
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    #source('BEER.R')
    
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
    

### Step2. Detect Batch Effect

    mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, CPU=2, GN=2000, SEED=1, MTTAG='^MT-')   

    # GNUM: the number of groups in each batch (default: 30)
    # PCNUM: the number of computated PCA subspaces (default: 50)
    # ROUND: batch-effect removal strength, positive integer (default: 1)
    # GN: the number of variable genes in each batch (default: 2000)

    # If you are facing "huge" batch effect, please try "REGBATCH":
    # mybeer <- BEER(DATA, BATCH, REGBATCH=TRUE, GNUM=30, PCNUM=50, ROUND=1, CPU=2, GN=2000, SEED=1, MTTAG='^MT-') 
    
    # Users can use "ReBEER" to adjust GNUM, PCNUM, and ROUND (it's faster than directly using BEER).
    # mybeer <- ReBEER(mybeer, GNUM=30, PCNUM=50, ROUND=1, CPU=2, SEED=1) 
    
    # Check selected PCs
    PCUSE=mybeer$select
    COL=rep('black',length(mybeer$cor))
    COL[PCUSE]='red'
    plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
    
    # Users can select PCA subspaces based on the distribution of "Rank Correlation" and "Linear Correlation". 
    # PCUSE=.selectUSE(mybeer, CUTR=0.7, CUTL=0.7, RR=0.5, RL=0.5)
    
    
    
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT1.png" width="400">
    
### Step3. Visualization 
    
#### Keep batch effect:
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT2.png" width="400">
    
    pbmc <- mybeer$seurat
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)  
    

    

#### Remove batch effect:

<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT3.png" width="400">

    pbmc <- mybeer$seurat
    PCUSE <- mybeer$select
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
    
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1) 
    
    
    
    
</br>
</br>
    
# II. Combine Multiple Batches

Download demo data: https://sourceforge.net/projects/beergithub/files/
   
### Step1. Load Data
    
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    #source('BEER.R')
    
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

    mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, CPU=2, GN=2000, SEED=1, MTTAG='^MT-' )

    # Check selected PCs
    PCUSE=mybeer$select
    COL=rep('black',length(mybeer$cor))
    COL[PCUSE]='red'
    plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
    

<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT4.png" width="400">

    
### Step3. Visualization 
        
#### Keep batch effect:
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT5.png" width="400">
    
    pbmc <- mybeer$seurat
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    
    



#### Remove batch effect:

<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT6.png" width="400">

    pbmc <- mybeer$seurat
    PCUSE <- mybeer$select
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
    
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)   
    
    
    
</br>   
</br>


# III. UMAP-based Clustering
   

    VEC=pbmc@reductions$umap@cell.embeddings

    # Here, we use K-means to do the clustering
    N=20
    set.seed(123)
    K=kmeans(VEC,centers=N)

    CLUST=K$cluster
    pbmc@meta.data$clust=CLUST
    DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
    

<img src="https://github.com/jumphone/BEER/raw/master/DATA/CLUST1.png" width="400">    


    # Or, manually select some cells

    ppp=DimPlot(pbmc, reduction.use='umap', pt.size=0.5)
    used.cells <- CellSelector(plot = ppp)
    

<img src="https://github.com/jumphone/BEER/raw/master/DATA/CLUST2.png" width="400">    

    # Press "ESC"
    
<img src="https://github.com/jumphone/BEER/raw/master/DATA/CLUST3.png" width="400">    
    
    markers <- FindMarkers(pbmc, ident.1=used.cells,only.pos=T)    
    head(markers, n=20)
    
    
</br>   
</br>

# IV. Combine scATAC-seq & scRNA-seq

Please go to the website of Seurat to download DEMO data: https://satijalab.org/seurat/v3.0/atacseq_integration_vignette.html

The main difference between BEER and Seurat (combine scRNA-seq & scATAC-seq): BEER doesn't use "imputation".

### Step1. Load Data

    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    #source('BEER.R')
    
    library(Seurat)
    
    peaks <- Read10X_h5("../data/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")

    activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, 
        annotation.file = "../data/Homo_sapiens.GRCh37.82.gtf", 
        seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)
         
    pbmc.rna <- readRDS("../data/pbmc_10k_v3.rds")
    
    D1=as.matrix(activity.matrix)
    D2=as.matrix(pbmc.rna@assays$RNA@counts)
    colnames(D1)=paste0('ATAC_', colnames(D1))
    colnames(D2)=paste0('RNA_', colnames(D2))
    DATA=.simple_combine(D1,D2)$combine
    BATCH=rep('RNA',ncol(DATA))
    BATCH[c(1:ncol(D1))]='ATAC'
    
 
 
### Step2. Detect Batch Effect

    mybeer <- BEER(DATA, BATCH, REGBATCH=TRUE, GNUM=30, PCNUM=50, ROUND=1, CPU=2, GN=5000, SEED=1, MTTAG='^MT-')
    saveRDS(mybeer, file='mybeer')
    
    # Users can use "ReBEER" to adjust parameters
    mybeer <- ReBEER(mybeer, GNUM=100, PCNUM=100, ROUND=3, SEED=1)
    
    PCUSE=mybeer$select
    COL=rep('black',length(mybeer$cor))
    COL[PCUSE]='red'
    plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))
    

<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT7.png" width="400">   


### Step3. Visualization 
    
#### Keep batch effect:
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT8.png" width="400">
    
    pbmc <- mybeer$seurat
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    
   

#### Remove batch effect:

    pbmc <- mybeer$seurat  
    PCUSE=mybeer$select
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
    
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    
    
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT9.png" width="400">
   
    pbmc@meta.data$celltype=rep(NA,length(pbmc@meta.data$batch))
    pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='RNA')]=pbmc.rna@meta.data$celltype
    
    DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)
    
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT10.png" width="400">

    saveRDS(mybeer, file='mybeer.final.RDS')
    

# It's not good enough !

### For further enhancement, please see [V. Batch-effect Removal Enhancement](#v-batch-effect-removal-enhancement).

</br>
</br>

# V. Batch-effect Removal Enhancement

Please install ComBat & BBKNN.

    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("sva")
    BiocManager::install("limma")

    Install bbknn in python: https://github.com/Teichlab/bbknn
    
    
This DEMO follows [IV. Combine scATAC-seq & scRNA-seq](#iv-combine-scatac-seq--scrna-seq)
    
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    #source('BEER.R')
    mybeer=readRDS('mybeer.final.RDS')
    pbmc.rna <- readRDS("../data/pbmc_10k_v3.rds")
    
    
### Use ComBat&BBKNN without BEER:

    pbmc <- mybeer$seurat
    PCUSE=c(1:ncol(pbmc@reductions$pca@cell.embeddings))
    pbmc=BEER.combat(pbmc) 
    umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
    pbmc@reductions$umap@cell.embeddings=umap
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
    
    
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT12.png" width="400">     

### Use ComBat&BBKNN with BEER:

    pbmc <- mybeer$seurat
    PCUSE=mybeer$select
    pbmc=BEER.combat(pbmc) 
    umap=BEER.bbknn(pbmc, PCUSE, NB=3, NT=10)
    pbmc@reductions$umap@cell.embeddings=umap
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=F)
     
    saveRDS(pbmc, file='seurat.enh.RDS')
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT13.png" width="400"> 
  
    pbmc@meta.data$celltype=rep(NA,length(pbmc@meta.data$batch))
    pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='RNA')]=pbmc.rna@meta.data$celltype
    DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)
    
    
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT14.png" width="400"> 


# VI. Transfer labels

This DEMO follows [V. Batch-effect Removal Enhancement](#v-batch-effect-removal-enhancement)
   
    pbmc@meta.data$celltype=rep(NA,length(pbmc@meta.data$batch))
    pbmc@meta.data$celltype[which(pbmc@meta.data$batch=='RNA')]=pbmc.rna@meta.data$celltype
    #DimPlot(pbmc, reduction.use='umap', group.by='celltype', pt.size=0.1,label=T)
    
    #######
    VEC=pbmc@reductions$umap@cell.embeddings
    set.seed(123)
    N=150
    K=kmeans(VEC,centers=N)
    pbmc@meta.data$kclust=K$cluster   
    #DimPlot(pbmc, reduction.use='umap', group.by='kclust', pt.size=0.1,label=T)

    pbmc@meta.data$transfer=rep(NA, length(pbmc@meta.data$celltype))
    TMP=cbind(pbmc@meta.data$celltype, pbmc@meta.data$kclust)
    
    KC=unique(pbmc@meta.data$kclust)
    i=1
    while(i<=length(KC)){
        this_kc=KC[i]
        this_index=which(pbmc@meta.data$kclust==this_kc)
        this_tb=table(pbmc@meta.data$celltype[this_index])
        if(length(this_tb)!=0){
            this_ct=names(this_tb)[which(this_tb==max(this_tb))[1]]
            pbmc@meta.data$transfer[this_index]=this_ct}
        i=i+1}
        
    pbmc@meta.data$tf.ct=pbmc@meta.data$celltype
    NA.index=which(is.na(pbmc@meta.data$celltype))
    pbmc@meta.data$tf.ct[NA.index]=pbmc@meta.data$transfer[NA.index]
    
    ######
    RNA.cells=colnames(pbmc)[which(pbmc@meta.data$batch=='RNA')]
    ATAC.cells=colnames(pbmc)[which(pbmc@meta.data$batch=='ATAC')]
    
    library(ggplot2)
    
    plot.all <- DimPlot(pbmc, reduction.use='umap', group.by='batch', 
        pt.size=0.1,label=F) + labs(title = "Batches")
    
    plot.ct <- DimPlot(pbmc,reduction.use='umap', group.by='tf.ct', 
        pt.size=0.1,label=T) + labs(title = "CellType")
    
    plot.rna <- DimPlot(pbmc, cells=RNA.cells,reduction.use='umap', 
        group.by='tf.ct', pt.size=0.1,label=T,plot.title='RNA.transfer') + labs(title = "RNA")
    
    plot.atac <- DimPlot(pbmc, cells=ATAC.cells,reduction.use='umap', 
        group.by='tf.ct', pt.size=0.1,label=T,plot.title='ATAC.transfer') + labs(title = "ATAC")
    
    CombinePlots(list(all=plot.all, ct=plot.ct, rna=plot.rna, atac=plot.atac))
    

<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT15.png" width="900"> 

</br>
</br>

# VII. Biological Interpretation

Please install "RITANdata" and "RITAN".

RITAN: https://bioconductor.org/packages/devel/bioc/vignettes/RITAN/inst/doc/enrichment.html

This DEMO follows [IV. Combine scATAC-seq & scRNA-seq](#iv-combine-scatac-seq--scrna-seq)

    library(RITANdata)
    library(RITAN)
    
    PCUSE <- mybeer$select
    PCALL <- c(1:length(mybeer$cor))
    PCnotUSE <- PCALL[which(!PCALL %in% PCUSE)]
    
    LD=mybeer$seurat@reductions$pca@feature.loadings
    GNAME=rownames(LD)
    
    N=100
    getPosAndNegTop <- function(x){
        O=c(order(x)[1:N],order(x)[(length(x)-(N-1)):length(x)])
        G=GNAME[O]
        return(G)
        }
    
    GMAT=apply(LD,2,getPosAndNegTop)
    colnames(GMAT)=paste0(colnames(GMAT),'_R_',round(mybeer$cor,1),"_L_",round(mybeer$lcor,1))
    GMAT=toupper(GMAT)
    
    GMAT=GMAT[,PCnotUSE]
    #GMAT=GMAT[,PCUSE]
  
    study_set=list()
    TAG=colnames(GMAT)
    i=1
    while(i<=ncol(GMAT)){
         study_set=c(study_set,list(GMAT[,i]))
         i=i+1
         }  
         
    names(study_set)=TAG
    #names(geneset_list)
    resources=c('KEGG_filtered_canonical_pathways','MSigDB_Hallmarks')
    
    e <- term_enrichment_by_subset( study_set, q_value_threshold = 1e-5, 
                                resources = resources,
                                all_symbols = cached_coding_genes )
    
    plot( e, show_values = FALSE, label_size_y = 7, label_size_x = 7, cap=10 )
    
</br> 

<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOTEB.png" width="600"> 

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

  
