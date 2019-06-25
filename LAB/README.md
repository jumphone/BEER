
<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/BEERLAB_LOGO.png" width="300">

# Welcome to BEER's laboratory !

## We are making new BEER here.

Author: Feng Zhang


# Requirement:

    #R >=3.5
    install.packages('Seurat') # >=3.0     
    install.packages('pcaPP') 
    install.packages('mclust')
    install.packages('igraph')
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("sva")
    BiocManager::install("limma")
    
Install bbknn in python: https://github.com/Teichlab/bbknn



# Usage:

* [I. Combine Batches](#i-Combine-Batches)
* [II. Tune-Up](#ii-tune-up)
* [III. Biological Interpretation](#iii-biological-interpretation)
* [IV. UMAP-based Clustering](#iv-UMAP-based-Clustering)
 
</br>

# I. Combine Batches

Download demo data: https://sourceforge.net/projects/beergithub/files/
   
### Step1. Load Data
    
    source('https://raw.githubusercontent.com/jumphone/BEER/master/LAB/BEER.lab.R')
    #source('BEER.lab.R')
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

    mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=3, CPU=2, SEED=1 )
   
    # GNUM: the number of groups in each batch
    # PCNUM: the number of computated PCA subspaces  
    # ROUND: the strength of batch-effect removal
    
    # Users can use "ReBEER" to adjust GNUM, PCNUM, and ROUND (it's faster than directly using BEER).
    mybeer <- ReBEER(mybeer, GNUM=30, PCNUM=50, ROUND=1, CPU=2, SEED=1)
    
    # Check selected PCs
    PCUSE=mybeer$select
    COL=rep('black',length(mybeer$cor))
    COL[PCUSE]='red'
    plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))

<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB0.png" width="400">

    
### Step3. Visualization 
        
#### Keep batch effect:
  
<img src="https://github.com/jumphone/BEER/raw/master/DATA/PLOT5.png" width="400">
    
    pbmc <- mybeer$seurat
    ALLPC <- 1:length(mybeer$cor)   
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = ALLPC, check_duplicates=FALSE)
    
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    
    #DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
    


#### Remove batch effect:

<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB1.png" width="400">

    pbmc <- mybeer$seurat
    PCUSE <- mybeer$select
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
    
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)    
    #DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
    
    
</br>   
</br>

# II. Tune-Up

If you need a "Tune-Up", please install ComBat & BBKNN.

ComBat: 

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("sva")
    BiocManager::install("limma")

BBKNN: https://github.com/Teichlab/bbknn.
 
### Solution 1: ComBat  

    pbmc <- mybeer$seurat
    PCUSE=mybeer$select
    
    pbmc=BEER.combat(pbmc) 
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=T)

<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB2.png" width="400">

### Solution 2: BBKNN  
 
    pbmc <- mybeer$seurat
    PCUSE=mybeer$select

    umap=BEER.bbknn(pbmc, PCUSE, NB=10, NT=10)
    
    pbmc@reductions$umap@cell.embeddings=umap
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=T)

<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB3.png" width="400">

### Solution 3: ComBat + BBKNN 

    pbmc <- mybeer$seurat
    PCUSE=mybeer$select
    
    pbmc=BEER.combat(pbmc) 
    
    umap=BEER.bbknn(pbmc, PCUSE, NB=10, NT=10)
    
    pbmc@reductions$umap@cell.embeddings=umap
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1,label=T)
    
<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB4.png" width="400">


</br>
</br>

# III. Biological Interpretation

Please install "RITANdata" and "RITAN".

RITAN: https://bioconductor.org/packages/devel/bioc/vignettes/RITAN/inst/doc/enrichment.html

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
    GMAT=GMAT[,PCnotUSE]
    
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

<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB5.png" width="800">

</br>   
</br> 
    
# IV. UMAP-based Clustering
   
<img src="https://github.com/jumphone/BEER/blob/master/LAB/img/LAB6.png" width="400"> 

    VEC=pbmc@reductions$umap@cell.embeddings
    
    # Here, we use K-means to do the clustering
    N=30
    set.seed(123)
    K=kmeans(VEC,centers=20)
    
    CLUST=K$cluster
    pbmc@meta.data$clust=CLUST
    DimPlot(pbmc, reduction.use='umap', group.by='clust', pt.size=0.5,label=TRUE)
    
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

  
