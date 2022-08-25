
# Batch EffEct Remover for single-cell data (BEER)
# Version: 0.1.9
# Author: Feng Zhang
# Date: Mar. 9, 2021
# For Seurat 3
#
#source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

print('Welcome to BEER (v0.1.9)!')

library(Seurat)
library(sva)
library(limma)

############################################################################################
############################################################################################
CORMETHOD='spearman'

####################################

#.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
#    exp_sc_mat=exp_sc_mat1
#    exp_ref_mat=exp_sc_mat2
#    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
#    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
#    gene_sc=rownames(exp_sc_mat)
#    gene_ref=rownames(exp_ref_mat)
#    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
#    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
#    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
#    colname_sc=colnames(exp_sc_mat)
#    colname_ref=colnames(exp_ref_mat)
#    OUT=list()
#    OUT$exp_sc_mat1=exp_sc_mat
#    OUT$exp_sc_mat2=exp_ref_mat
#    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
#    return(OUT)
#    }
######################

####################################
# Problem too large
# 2021.9.17
as_matrix <- function(mat){
 
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
    
  for (i in seq_along(val)){
      tmp[row_pos[i],col_pos[i]] <- val[i]
  }
    
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
###############

#######################

# 2019.11.01

.check_rep <- function(MAT){
    MAT=as.matrix(MAT)
    RNAME=rownames(MAT)
    USED_NAME=names(which(table(RNAME)==1))
    MAT=MAT[which(rownames(MAT) %in% USED_NAME),]
    return(MAT)
    }


.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){    
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){ 
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(0,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(0,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################ 
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }

######################


.run_seurat <-function(DATA,PCNUM=50, GN=2000,  SEED=123, N=2){
    DATA=DATA
    PCNUM=PCNUM
    GN=GN
    SEED=SEED
    N=N
    ###############
    pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL")  
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc=FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = GN)
    pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc))
    print('Calculating PCs ...')
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)
    ################
    return(pbmc)
    }
    
################################    
    
    


.generate_ref <- function(exp_sc_mat, TAG, min_cell=1, refnames=FALSE){
    NewRef=c()
    TAG[,2]=as.character(TAG[,2])
    if(refnames==FALSE){
        refnames=names(table(TAG[,2]))}
        else{refnames=refnames}
    outnames=c()
    for(one in refnames){
        this_col=which(TAG[,2]==one)
        if(length(this_col)>= min_cell){
            outnames=c(outnames,one)
            if(length(this_col) >1){
                #this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
                this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
                }
                else{this_new_ref = exp_sc_mat[,this_col]}
            NewRef=cbind(NewRef,this_new_ref)
            }
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }



.generate_agg <- function(exp_sc_mat, TAG, print_step=100){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG
    
    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))
    
    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){   
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1       
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }

############
.generate_mean <- function(exp_sc_mat, TAG, print_step=100){
    print_step=print_step
    exp_sc_mat=exp_sc_mat
    TAG=TAG
    
    NewRef=matrix(0,ncol=length(unique(TAG)),nrow=nrow(exp_sc_mat))
    
    TAG=as.character(TAG)
    refnames=unique(TAG)
    total_num=length(refnames)
    outnames=c()
    i=1
    while(i<=length(refnames)){
        one=refnames[i]
        this_col=which(TAG==one)
        outnames=c(outnames,one)
        if(length(this_col) >1){   
            #this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            this_new_ref=apply(exp_sc_mat[,this_col],1,mean)
            }else{
            this_new_ref = exp_sc_mat[,this_col]
            }
        NewRef[,i]=this_new_ref
        if(i%%print_step==1){print(paste0(i,' / ' ,total_num ))}
        i=i+1       
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }
############################################################################################
############################################################################################



.getGroup <- function(X,TAG,GNUM){
    
    print('Get group for:')
    print(TAG)
    
    #D=dist(X)
    #H=hclust(D)
    #CLUST=cutree(H,k=GNUM)
    CLUST=kmeans(X,centers=GNUM,iter.max =100)$cluster
    GROUP=paste0(TAG,'_',as.character(CLUST))
    
    print('Group Number:')
    print(length(table(GROUP)))
    return(GROUP)
}




####################

.getVPall<- function(pbmc, ROUND){
    
    pbmc=pbmc
    ROUND=ROUND
    print('Finding MN pairs...')
    ################
    REF=.generate_agg(pbmc@assays$RNA@data, pbmc@meta.data$group)
    VREF=REF
    CVREF=cor(VREF,method=CORMETHOD)
    orig.CVREF=CVREF
    #ROUND=3

    UBATCH=unique(pbmc@meta.data$batch)

    .get_batch<-function(x){
        y=unlist(strsplit(x,'_'))[1]
        return(y)
        } 
    group_batch=apply(as.matrix(colnames(CVREF)),1,.get_batch)

    .getMN <- function(this_cor_mat){
        VP=c()
        i=1
        while(i<=nrow(this_cor_mat)){
            this_p1=rownames(this_cor_mat)[i]
            j=1
            while(j<=ncol(this_cor_mat)){
                this_p2=colnames(this_cor_mat)[j]  
                this_cor=this_cor_mat[i,j]
                if( this_cor!= -99999  & this_cor==max(this_cor_mat[i,]) & this_cor==max(this_cor_mat[,j])){                 
                    VP=cbind(VP,c(this_p1,this_p2))
                    }                
                j=j+1}        
            i=i+1}
        return(VP)
    
        }
    
    
    
    VP=c()
    I=1    
    while(I<=ROUND){
            
        if(length(VP)!=0){   
            i=1
            while(i<=ncol(VP)){
                p1=VP[1,i]
                p2=VP[2,i]
                b1=.get_batch(p1)
                b2=.get_batch(p2)
                b1i=which(group_batch==b1)
                b2i=which(group_batch==b2)
                #CVREF[which(rownames(CVREF)==p1), which(colnames(CVREF)==p2)]=-99999
                #CVREF[which(rownames(CVREF)==p2), which(colnames(CVREF)==p1)]=-99999    
                CVREF[which(rownames(CVREF)==p1),b2i]=-99999
                CVREF[which(rownames(CVREF)==p2),b1i]=-99999
                CVREF[b2i,which(rownames(CVREF)==p1)]=-99999
                CVREF[b1i,which(rownames(CVREF)==p2)]=-99999            
                i=i+1}
            }
    
        
        
        i=1
        while(i<length(UBATCH)){
            j=i+1
            while(j<=length(UBATCH)){
                b1=UBATCH[i]
                b2=UBATCH[j]
                b1i=which(group_batch==b1)
                b2i=which(group_batch==b2)
                this_cor_mat=CVREF[b1i,b2i]
                this_vp=.getMN(this_cor_mat)    
                VP=cbind(VP,this_vp)
                
                ########################
                                
                ########################    

                
                j=j+1}
            i=i+1
            }
        
        ########################
        
        ######################## 
        print('ROUND:')
        print(I)
        I=I+1
        }       
        
    
    VP=t(VP)
    VP=unique(VP)
    #######################
    
    print('Number of MN pairs:')
    print(nrow(VP))
    return(VP)
    }





################ 2019.06.18 #####

.evaluateProBEER <- function(DR, GROUP, VP){
    
    OUT=list() 
    VALID_PAIR=VP
    ALL_COR=c()   
    ALL_PV=c() 
    ALL_LCOR=c()
    ALL_LPV=c()
    ALL_LC1=c()
    ALL_LC2=c()
    print('Evaluating PCs ...')
    print('Start')
    THIS_DR=1
    while(THIS_DR<=ncol(DR)){
        THIS_PC = DR[,THIS_DR]
        
        lst1_quantile=c()
        lst2_quantile=c()
        i=1
        while(i<=nrow(VALID_PAIR)){
            this_pair=VALID_PAIR[i,]
            this_index1=which(GROUP %in% this_pair[1])
            this_index2=which(GROUP %in% this_pair[2])                   
            lst1_quantile=c(lst1_quantile,quantile(DR[this_index1,THIS_DR]))
            lst2_quantile=c(lst2_quantile,quantile(DR[this_index2,THIS_DR]))
                       
            i=i+1}
        
        this_test=cor.test(lst1_quantile, lst2_quantile, method=CORMETHOD)        
        this_cor=this_test$estimate
        this_pv=this_test$p.value

        this_test2=cor.test(lst1_quantile, lst2_quantile, method='pearson')        
        this_cor2=this_test2$estimate
        this_pv2=this_test2$p.value
        
        olddata=data.frame(lst1=lst1_quantile, lst2=lst2_quantile)
        fit=lm(lst1 ~lst2, data=olddata) 
        
        ALL_LC1=c(ALL_LC1, summary(fit)$coefficients[1,4])
        ALL_LC2=c(ALL_LC2, summary(fit)$coefficients[2,4])
        
        ALL_COR=c(ALL_COR, this_cor)
        ALL_PV=c(ALL_PV, this_pv) 
        ALL_LCOR=c(ALL_LCOR, this_cor2)
        ALL_LPV=c(ALL_LPV, this_pv2) 
        print(THIS_DR)
        
        THIS_DR=THIS_DR+1}
    
    
    OUT$cor=ALL_COR
    OUT$pv=ALL_PV
    OUT$fdr=p.adjust(ALL_PV,method='fdr')
    OUT$lc1=ALL_LC1
    OUT$lc2=ALL_LC2
    OUT$lcor=ALL_LCOR
    OUT$lpv=ALL_LPV
    OUT$lfdr=p.adjust(ALL_LPV,method='fdr')

    print('Finished!!!')
    return(OUT)
   }


BEER <- function(DATA, BATCH,  GNUM=30, PCNUM=50, GN=2000, CPU=4, COMBAT=TRUE, print_step=10, SEED=123, N=2, ROUND=1, RMG=NULL){

    set.seed( SEED)
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=DATA
    BATCH=BATCH
    RMG=RMG
    COMBAT=COMBAT
    COMBAT.EXP=NULL
    
    
    require(stringi)
    BATCH=stri_replace_all(BATCH, '.',fixed='_')
    CPU=CPU
    GNUM=GNUM
    PCNUM=PCNUM
    UBATCH=unique(BATCH)
    ROUND=ROUND
    
    GN=GN
    N=N
    print_step=print_step
    
    print('Group number (GNUM) is:')
    print(GNUM)
    print('Varible gene number (GN) of each batch is:')
    print(GN)
    print('ROUND is:')
    print(ROUND)
    
    VARG=c()
    i=1
    for(this_batch in UBATCH){
        print(i)
        i=i+1
        print(this_batch)
        this_pbmc=CreateSeuratObject(counts = DATA[,which(BATCH==this_batch)], min.cells = 0, 
                                 min.features = 0, project = this_batch)
        this_pbmc <- NormalizeData(object = this_pbmc, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
        this_pbmc <- FindVariableFeatures(object = this_pbmc, selection.method = "vst", nfeatures = GN)  
        this_varg=VariableFeatures(object = this_pbmc)
        VARG=c(VARG, this_varg)
        }
    VARG=unique(VARG)

    print('Total varible gene number (GN) is:')
    print(length(VARG))

    pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL") 
    pbmc@meta.data$batch=BATCH
    VariableFeatures(object = pbmc)=VARG

    ########
    if(!is.null(RMG)){
        print('Total removed gene number is:')
        print(length(RMG))
        VariableFeatures(object = pbmc)=VariableFeatures(object = pbmc)[which(! VariableFeatures(object = pbmc) %in% RMG)]
        print('Total used gene number is:')
        print(length(VariableFeatures(object = pbmc)))
        }
    ##########   

    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


    if(COMBAT==FALSE){
        pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc))
    }else{
        ##############
        library(sva)
        library(limma)
        pheno = data.frame(batch=as.matrix(BATCH))
        orig.data=pbmc@assays$RNA@data
        used.gene.index=which(rownames(orig.data) %in% VARG)
        #edata = as.matrix(orig.data)[used.gene.index,]
        edata = as_matrix(orig.data)[used.gene.index,]
        batch = pheno$batch
        modcombat = model.matrix(~1, data=pheno)
        combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
        rownames(combat_edata)=rownames(edata)
        colnames(combat_edata)=colnames(edata)
        combat_edata=as.matrix(combat_edata)
        combat_edata[which(combat_edata<0)]=0
        combat_edata[which(is.na(combat_edata))]=0
        pbmc@assays$RNA@data=combat_edata
        ######
        pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc))
        ######
        pbmc@assays$RNA@data=orig.data    
        COMBAT.EXP=combat_edata
        #################
        rm(edata)
        rm(combat_edata)
        rm(orig.data)
        gc()
    }
    print('Calculating PCs ...')
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)


    DR=pbmc@reductions$umap@cell.embeddings
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=DR[this_index,]
        #CNUM=max(c(5, round(length(this_index)/GNUM) ))
        this_gnum=min(GNUM, (length(this_index)-1))
        this_group=.getGroup(this_one,this_batch,this_gnum)
        GROUP[this_index]=this_group
    }


    pbmc@meta.data$group=GROUP
    
    ##########
    #VP=.getVPnet(pbmc, ROUND)
    VP=.getVPall(pbmc, ROUND)
    ##########
    
    DR=pbmc@reductions$pca@cell.embeddings  

    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='V1'
    MAP[which(GROUP %in% VP[,2])]='V2'
    pbmc@meta.data$map=MAP
    
    OUT=.evaluateProBEER(DR, GROUP, VP)

    RESULT=list()
    RESULT$seurat=pbmc
    RESULT$vp=VP
    #RESULT$vpcor=VP_OUT$cor
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=OUT$fdr
    RESULT$lcor=OUT$lcor
    RESULT$lpv=OUT$lpv
    RESULT$lc1=OUT$lc1
    RESULT$lc2=OUT$lc2
    RESULT$lfdr=OUT$lfdr
    
    ################
    RESULT$ROUND=ROUND
    RESULT$COMBAT=COMBAT
    RESULT$COMBAT.EXP=COMBAT.EXP
    RESULT$RMG=RMG
    RESULT$GNUM=GNUM
    RESULT$GN=GN
    RESULT$PCNUM=PCNUM
    RESULT$SEED=SEED
    RESULT$N=N
    RESULT$APP='BEER'   
    ###############
    
    
    PCUSE=which( (rank(RESULT$cor)>=length(RESULT$cor)/2 | RESULT$cor>0.7 )    & 
                (rank(RESULT$lcor) >=length(RESULT$cor)/2 | RESULT$lcor>0.7)   #&
                #p.adjust(RESULT$lc1,method='fdr') >0.05
               ) 
    
    RESULT$select=PCUSE
    
    print('############################################################################')
    print('BEER cheers !!! All main steps finished.')
    print('############################################################################')
    print(Sys.time())

    return(RESULT)
}

MBEER=BEER

.getUSE <-function(RESULT, CUTR=0.7,CUTL=0.7){
    PCUSE=which( (RESULT$cor>CUTR )    & (RESULT$lcor>CUTL) ) 
    return(PCUSE)
    }


.selectUSE <-function(RESULT, CUTR=0.7, CUTL=0.7, RR=0.5, RL=0.5){
    PCUSE=which( (rank(RESULT$cor)>=length(RESULT$cor)*RR | RESULT$cor>CUTR )    & 
                (rank(RESULT$lcor) >=length(RESULT$cor)*RR | RESULT$lcor>CUTL)  
               ) 
    return(PCUSE)
    }

####################



ReBEER <- function(mybeer,  GNUM=30, PCNUM=50,  CPU=4, print_step=10, SEED=123, N=2, ROUND=1, RMG=NULL){

    set.seed( SEED)
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=mybeer$seurat@assays$RNA@counts
    BATCH=mybeer$seurat@meta.data$batch
    GNUM=GNUM
    PCNUM=PCNUM
    RMG=RMG
    #MAXBATCH=MAXBATCH
    UBATCH=unique(BATCH)
    
    ROUND=ROUND
    N=N
    print_step=print_step
    
    pbmc=mybeer$seurat
    
    print('Group number (GNUM) is:')
    print(GNUM)
    print('Total varible gene number is:')
    print(length(VariableFeatures(object = pbmc)))
    print('ROUND is:')
    print(ROUND)
    ########
    if(!is.null(RMG)){
        print('Total removed gene number is:')
        print(length(RMG))
        VariableFeatures(object = pbmc)=VariableFeatures(object = pbmc)[which(! VariableFeatures(object = pbmc) %in% RMG)]
        print('Total used gene number is:')
        print(length(VariableFeatures(object = pbmc)))
        }
    ##########
    VARG = VariableFeatures(object = pbmc)
    print('Calculating PCs ...')
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)
    ########
    
    DR=pbmc@reductions$umap@cell.embeddings
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=DR[this_index,]
        #CNUM=max(c(5, round(length(this_index)/GNUM) ))
        this_gnum=min(GNUM, (length(this_index)-1))
        this_group=.getGroup(this_one,this_batch,this_gnum)
        GROUP[this_index]=this_group
    }

    pbmc@meta.data$group=GROUP
    
     
    ##########
    #VP=.getVPnet(pbmc, ROUND)
    VP=.getVPall(pbmc, ROUND)
    ##########
    
    DR=pbmc@reductions$pca@cell.embeddings  

    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='V1'
    MAP[which(GROUP %in% VP[,2])]='V2'
    pbmc@meta.data$map=MAP
    
       
    
    OUT=.evaluateProBEER(DR, GROUP, VP)

    RESULT=list()
    RESULT$seurat=pbmc
    RESULT$vp=VP
    #RESULT$mst=MST
    #RESULT$vpcor=VP_OUT$cor
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=OUT$fdr
    RESULT$lcor=OUT$lcor
    RESULT$lpv=OUT$lpv
    RESULT$lc1=OUT$lc1
    RESULT$lc2=OUT$lc2
    RESULT$lfdr=OUT$lfdr
    
    ################
    RESULT$ROUND=ROUND
    RESULT$GNUM=GNUM
    RESULT$PCNUM=PCNUM
    RESULT$SEED=SEED
    RESULT$N=N
    RESULT$APP='ReBEER'   
    ###############

    #########
       
    PCUSE=which( (rank(RESULT$cor)>=length(RESULT$cor)/2 | RESULT$cor>0.7 )    & 
                (rank(RESULT$lcor) >=length(RESULT$cor)/2 | RESULT$lcor>0.7)   #&
                #p.adjust(RESULT$lc1,method='fdr') >0.05
               ) 
    
    RESULT$select=PCUSE
    
    print('############################################################################')
    print('BEER cheers !!! All main steps finished.')
    print('############################################################################')
    print(Sys.time())

    return(RESULT)
}



#########################
BEER.combat <- function(pbmc){
    
    #mybeer=mybeer

    pbmc=pbmc#mybeer$seurat
    batch=as.character(pbmc@meta.data$batch)
    
    pca=pbmc@reductions$pca@cell.embeddings 
    library(sva)
    library(limma)
    pheno = data.frame(batch=as.matrix(batch))
    edata = t(pca)
    batch = pheno$batch
    modcombat = model.matrix(~1, data=pheno)
    combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    ttt=t(combat_edata)
    colnames(ttt)=colnames(pbmc@reductions$pca@cell.embeddings)
    rownames(ttt)=rownames(pbmc@reductions$pca@cell.embeddings)

    pca=ttt
    pbmc@reductions$pca@cell.embeddings=pca
    return(pbmc)
    }


BEER.bbknn <- function(pbmc, PCUSE, NB=3, NT=10, DM=2){
  
    NB=NB
    NT=NT
    DM=DM
    #mybeer=mybeer

    pbmc=pbmc#mybeer$seurat
    PCUSE=PCUSE
    
    #pbmc=mybeer$seurat
    batch=as.character(pbmc@meta.data$batch)
    
    pca.all=pbmc@reductions$pca@cell.embeddings
    pca.use=pbmc@reductions$pca@cell.embeddings[,PCUSE]
    
        
    library(reticulate)
    #use_python("C:\Users\cchmc\Anaconda3\python")
    
    
    anndata = reticulate::import("anndata",convert=FALSE) #anndata==0.7
    bbknn = reticulate::import("bbknn", convert=FALSE)
    #sc = reticulate::import("scanpy.api",convert=FALSE) #scanpy==1.5.1
    sc = reticulate::import("scanpy",convert=FALSE) #scanpy
    
    adata = anndata$AnnData(X=pca.all, obs=batch)
    PCNUM=ncol(pca.use)

    sc$tl$pca(adata, n_comps=as.integer(PCNUM))
    adata$obsm$X_pca = pca.use
    #NB=50
    #NT=10
    #print(NB)
    bbknn$bbknn(adata,batch_key=0,neighbors_within_batch=as.integer(NB),n_pcs=as.integer(PCNUM), annoy_n_trees =as.integer(NT))
    #bbknn$bbknn(adata,batch_key=0,neighbors_within_batch=as.integer(NB),n_pcs=as.integer(PCNUM), n_trees =as.integer(NT))
    #bbknn$bbknn(adata,batch_key=0, n_pcs=as.integer(PCNUM))
    
    #sc$tl$umap(adata)
    sc$tl$umap(adata, n_components=as.integer(DM))
    
    #umap = py_to_r(adata$obsm$X_umap)
    umap = py_to_r(adata$obsm['X_umap'])
    rownames(umap)=rownames(pbmc@reductions$umap@cell.embeddings)
    colnames(umap)=paste0('UMAP_',c(1:ncol(umap)))#colnames(pbmc@reductions$umap@cell.embeddings)
  
    return(umap)
    }

#######2019.07.17



.readTable <- function(PATH, SEP='\t', UP=FALSE){
    SEP=SEP
    PATH=PATH
    UP=UP
    DATA=read.table(file=PATH,sep=SEP,header=TRUE,row.names=NULL)
    DATA=apply(DATA,2,as.character)
    ###########
    if(UP==TRUE){DATA[,1]=toupper(DATA[,1])}
    ###########
    TAB=table(as.character(DATA[,1]))
    UNIQ=names(TAB)[which(TAB==1)]
    DATA=DATA[which(DATA[,1] %in% UNIQ),]
    RN=DATA[,1]
    DATA=DATA[,c(2:ncol(DATA))]
    DATA=apply(DATA,2,as.numeric)
    rownames(DATA)=RN
    return(DATA)
    }

.writeTable <- function(DATA, PATH, SEP='\t',TL='OUTPUT'){
    DATA=DATA
    PATH=PATH
    SEP=SEP
    TL=TL
    OUT=cbind(rownames(DATA),DATA)
    colnames(OUT)[1]=TL
    write.table(OUT, file=PATH,sep=SEP,row.names=FALSE,col.names=TRUE,quote=FALSE)
    }

###########





.getPos <- function(x){
    x=x
    y=length(which(x>0))
    return(y)
    }

.getBatchPos <- function(DATA,BATCH, N=50,SEED=123){
    #######
    DATA=DATA
    BATCH=BATCH
    SEED=SEED
    N=N
    ######
    UB=unique(BATCH)
    set.seed(SEED)
    BAT=c()
    POS=c()
    MED=c()
    for(this_batch in UB){
        batch_index=which(BATCH==this_batch)
        this_index=sample(batch_index, N,replace=TRUE)
        this_pos=apply(DATA[,this_index],2,.getPos)
        this_pos_m=median(this_pos)
        MED=c(MED,this_pos_m)
        POS=c(POS,this_pos)
        BAT=c(BAT,rep(this_batch,N))
        }  
    OUT=list()
    OUT$pos=POS
    OUT$bat=BAT
    OUT$ub=UB
    names(MED)=UB
    OUT$med=MED
    return(OUT)
   }



BEER.AGG <- function(DATA, BATCH, FOLD, PCNUM=50, GN=2000, CPU=4, print_step=10, SEED=123, N=2, RMG=NULL){
    DATA=DATA
    BATCH=BATCH
    FOLD=FOLD
    CPU=CPU
    SEED=SEED
    RMG=RMG
    require(stringi)
    BATCH=stri_replace_all(BATCH, '.',fixed='_')
    PCNUM=PCNUM
    UBATCH=unique(BATCH)
    GN=GN
    N=N
    print_step=print_step
    print('Varible gene number (GN) of each batch is:')
    print(GN)
     
    VARG=c()
    i=1
    for(this_batch in UBATCH){
        print(i)
        i=i+1
        print(this_batch)
        this_pbmc=CreateSeuratObject(counts = DATA[,which(BATCH==this_batch)], min.cells = 0, 
                                 min.features = 0, project = this_batch)
        this_pbmc <- NormalizeData(object = this_pbmc, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
        this_pbmc <- FindVariableFeatures(object = this_pbmc, selection.method = "vst", nfeatures = GN)  
        this_varg=VariableFeatures(object = this_pbmc)
        VARG=c(VARG, this_varg)
        }
    VARG=unique(VARG)
    
    print('Total varible gene number (GN) is:')
    print(length(VARG))

    pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL") 
    pbmc@meta.data$batch=BATCH
    VariableFeatures(object = pbmc)=VARG
    
    if(!is.null(RMG)){
        print('Total removed gene number is:')
        print(length(RMG))
        VariableFeatures(object = pbmc)=VariableFeatures(object = pbmc)[which(! VariableFeatures(object = pbmc) %in% RMG)]
        print('Total used gene number is:')
        print(length(VariableFeatures(object = pbmc)))
        }
    
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc))
    
    print('Calculating PCs ...')
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)
    #DimPlot(pbmc)
    VEC=pbmc@reductions$umap@cell.embeddings
    
    UB=unique(BATCH)
    TAG=c()
    
    ###########
    for(this_batch in UB){
         #this_batch=UB[31]
        
         this_index=which(BATCH == this_batch)
         this_fold=FOLD[which(names(FOLD)==this_batch)]
        
         if(this_fold==1){
             this_tag=paste0(this_batch,'...',c(1:length(this_index)))
             }else{      
             this_vec=VEC[this_index,]
             set.seed(SEED)
             this_n=round(length(this_index)/this_fold)
             this_km=kmeans(this_vec,centers=this_n)
             this_cl=this_km$cluster
             this_tag=paste0(this_batch,'...',this_cl)          
             }
         TAG=c(TAG,this_tag)
        
         }
    
    DATA.AGG=.generate_agg(DATA, TAG)
    
    .getAggBatch <- function(x){
        y=unlist(strsplit(x, "\\.\\.\\."))[1]
        return(y)
     }
    
    CN=colnames(DATA.AGG)
    DATA.AGG.BATCH=apply(matrix(CN,ncol=1),1,.getAggBatch)
    
    RESULT=list()
    RESULT$data.agg=DATA.AGG
    RESULT$data.agg.batch=DATA.AGG.BATCH
    RESULT$vec=VEC
    RESULT$tag=TAG
    RESULT$cell=colnames(pbmc)
    return(RESULT)
    
    }

#######################
#2019.07.25
BEER.SMOOTH<-function(EXP,VEC,N=3,print_step=10,SEED=123){
    EXP=as.matrix(EXP)
    VEC=as.matrix(VEC)
    SEED=SEED
    N=N
    print_step=print_step
    ########
    set.seed(SEED)
    ########
    EXP.SM=matrix(0,ncol=ncol(EXP),nrow=nrow(EXP))
    D=dist(VEC)
    D=as.matrix(D)
    i=1
    while(i<=ncol(EXP)){
        this_index= order(D[,i])[1:N]
        #EXP.SM[,i]=apply(EXP[,this_index], 1, weighted.mean, 2**c(N:1))
        EXP.SM[,i]=apply(EXP[,this_index], 1, mean)
        if(i %% print_step==1){print(paste0(i,' / ',ncol(EXP)))}
        i=i+1}
    ###########
    rownames(EXP.SM)=rownames(EXP)
    colnames(EXP.SM)=colnames(EXP)
    RESULT=list()
    RESULT$exp.smooth=EXP.SM
    RESULT$distance=D
    ###########
    return(RESULT)
    }



######
#2019.0726

.combat <- function(EXP, BATCH){
    library(sva)
    library(limma)
    pheno = data.frame(batch=as.matrix(BATCH))
    edata = as.matrix(EXP)
    batch = pheno$batch
    modcombat = model.matrix(~1, data=pheno)
    combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    return(combat_edata)
    }

.norm_exp<-function(x){
    y=x
    y[which(x<0)]=0
    y=x/sum(x)
    y=y*1000000
    return(y)
    }
    

#####
#2019.0806
.set_python <- function(PATH){
    library(reticulate)
    #use_python("C:/Users/cchmc/Anaconda3/python")
    use_python(PATH,required=TRUE)   
    }

####
#2019.10.24
BEER.IMP <- function(DATA, VEC, print_step=100, CUTOFF=0.2){
    DATA=DATA
    VEC=VEC
    print_step=print_step
    CUTOFF=CUTOFF
    #############################
    DATA=as.matrix(DATA)
    NC=ncol(DATA)

    DIST=dist(VEC)
    DIST=as.matrix(DIST)

    DIST.NUM=as.numeric(DIST)
    E.DIST.NUM=ecdf(DIST.NUM)
    PV=apply(DIST, 2, E.DIST.NUM)
    LOG.PV=-log(PV,2)
    LOG.PV[which(PV>CUTOFF)]=0
    ###############


    NEW.DATA=matrix(0,ncol=ncol(DATA),nrow=nrow(DATA))
    rownames(NEW.DATA)=rownames(DATA)
    colnames(NEW.DATA)=colnames(DATA)
    
    i=1
    while(i<=NC){
        NEW.DATA[,i]= DATA %*% as.matrix(LOG.PV[,i] ,ncol=1)
        NEW.DATA[,i] = NEW.DATA[,i] / sum(LOG.PV[,i])
        if(i %% print_step==1){print(i)}
        i=i+1}

    #########################
    
    return(NEW.DATA)
    }

##########################

########################
#2019.11.19

.getGSEAinput <- function( DATA, TAG, PATH ){
    DATA=as.matrix(DATA)
    TAG=TAG
    VAR=apply(DATA,1,var)
    DATA=DATA[which(VAR>0),]
    ########################
    EXP.FILE=paste0(PATH,'.EXP.txt')
    colnames(DATA)=paste0(TAG,'_',colnames(DATA))
    OUT=cbind(toupper(rownames(DATA)),rep('NO',nrow(DATA)),DATA)
    colnames(OUT)[c(1,2)]=c('GENE','DESCRIPTION')
    write.table(OUT, EXP.FILE, sep='\t',quote=F,row.names=F,col.names=T)
    
    #########################
    PT.FILE=paste0(PATH,'.PT.cls')
    PT=t(as.character(TAG))
    cat(paste0(length(TAG),' 2 1'),file=PT.FILE,sep="\n") 
    cat(paste(c('#',unique(TAG)),collapse=' '),file=PT.FILE,sep='\n',append=TRUE)
    cat(PT,file=PT.FILE,sep=' ',append=TRUE)
    ##########################  
    }

#######################################################

