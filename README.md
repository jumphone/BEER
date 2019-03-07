# BEER: Batch EffEct Remover of single-cell data

Author: Feng Zhang

Date: Mar. 7, 2019


# Requirement:

    #R >=3.5
    install.packages('Suerat')
    install.packages('pcaPP')

# Usage:

### Step1. Load Data

    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
    
    D1_address=https://raw.githubusercontent.com/jumphone/BEER/master/DATA/DATA1_MAT.txt
    D2_address=https://raw.githubusercontent.com/jumphone/BEER/master/DATA/DATA2_MAT.txt
    
    D1 <- read.table(D1_address, sep='\t', row.names=1, header=T)
    D2 <- read.table(D2_address, sep='\t', row.names=1, header=T)

    # "D1" & "D2" are UMI matrix (or FPKM, RPKM, TPM, PKM ...; Should not be gene-centric scaled data)
    # Rownames of "D1" & "D2" are gene names
    # Colnames of "D1" & "D2" are cell names 

### Step2. Detect Batch Effect

    mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)
    
    par(mfrow=c(1,2))
    plot(mybeer$cor, xlab='PCs', ylab='PCC', pch=16)
    plot(-log(mybeer$fdr,10), xlab='PCs', ylab='-log10(FDR)', pch=16)
    
### Step3. Visualization 
    
    pbmc <- mybeer$seurat
    
#### Keep batch effect:
    
    ALLPC <- 1:length(mybeer$cor)
    
    # UMAP:
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
    DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
    
    # tSNE:
    pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = ALLPC, do.fast = TRUE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
    DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)

#### Remove batch effect:

    PCUSE <- which(mybeer$cor>0.7 & mybeer$fdr<0.05)
    
    # UMAP:
    pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
    DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
    
    # tSNE:
    pbmc <- RunTSNE(object = pbmc, reduction.use='pca',dims.use = PCUSE, do.fast = TRUE, check_duplicates=FALSE)
    DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
    DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)
  
  
  
