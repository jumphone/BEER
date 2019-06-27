# Batch EffEct Remover for single-cell data (BEER)
# Author: Feng Zhang
# Date: June 18, 2019
# For Seurat 3
#
#source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


library(Seurat)
library(pcaPP)
library(igraph)
library(sva)
library(limma)

############################################################################################
############################################################################################

.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10, gene_check=FALSE){
    method=method
    CPU=CPU
    print_step=print_step
    #method = "pearson", "kendall", "spearman"
    ##################
    print('Gene number of exp_sc_mat1:')
    print(nrow(exp_sc_mat))
    print('Gene number of exp_sc_mat2:')
    print(nrow(exp_ref_mat))
    #################
    library(parallel)
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    ###############
    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))
    if(gene_check==TRUE){
    print('Press RETURN to continue:')
    scan();}
    ###################
    #Step 2. calculate prob
    SINGLE <- function(i){
        library('pcaPP')
        exp_sc = as.array(exp_sc_mat[,i])
        
        cor_list=rep(0,length(colname_ref))
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            
            if(method=='kendall'){this_cor=cor.fk(exp_sc,exp_ref)}
            else{
            this_cor=cor(exp_sc,exp_ref, method=method)}
            
            cor_list[j]=this_cor
            j=j+1}
        ################################
        if(i%%print_step==1){print(i)}
        return(cor_list)
        }
    #######################################
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    COR = c()
    for(cor_list in RUN){
        COR=cbind(COR, cor_list)}
    #######################################
    rownames(COR)=colname_ref
    colnames(COR)=colname_sc
    return(COR)
    }


.get_tag_max <- function(COR){
    RN=rownames(COR)
    CN=colnames(COR)
    TAG=cbind(CN,rep('NA',length(CN)))
    i=1
    while(i<=length(CN)){
        this_rn_index=which(COR[,i] == max(COR[,i]))[1]
        TAG[i,2]=RN[this_rn_index]
        i=i+1
        }
    colnames(TAG)=c('cell_id','tag')
    return(TAG)
    }

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
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
                this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
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



############################################################################################
############################################################################################



.data2one <- function(DATA, GENE, CPU=4, PCNUM=50, SEED=123,  PP=30){
    
    if(PP>5 & ncol(DATA)<=300 & ncol(DATA)>100){PP=5;print('The cell number is too small! The perplexity is changed to 5 !')} 
    if(PP>3 & ncol(DATA)<=100){PP=3;print('The cell number is too small! The perplexity is changed to 3 !')}
    
    
    PCUSE=1:PCNUM
    print('Start')
    library(Seurat)
    print('Step1.Create Seurat Object...')
    DATA =CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL")   
    print('Step2.Normalize Data...')
    DATA <- NormalizeData(object = DATA, normalization.method = "LogNormalize", scale.factor = 10000)
    print('Step3.Scale Data...')
    DATA <- ScaleData(object = DATA, features = GENE, vars.to.regress = c("nCount_RNA"), num.cores=CPU, do.par=TRUE)
    print('Step4.PCA...')
    DATA <- RunPCA(object = DATA, seed.use=SEED, npcs=PCNUM, features = GENE, ndims.print=1,nfeatures.print=1)
    print('Step5.One-dimention...')
    DATA <- RunTSNE(object = DATA, seed.use=SEED, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1,  perplexity= PP)
    DR=DATA@reductions$tsne@cell.embeddings
    print('Finished!!!')
    return(DR)
    }

.getGroupOld <- function(X,TAG,CNUM=10){
    DR=X
    RANK=rank(DR,ties.method='random')
    CUTOFF=CNUM 
    GROUP=rep('NA',length(RANK))
    i=1
    j=1
    while(i<=length(RANK)){
        GROUP[which(RANK==i)]=paste0(TAG,'_',as.character(j))
        #if(i%%CUTOFF==1){j=j+1;print(j)}
        ########################
        if(CUTOFF!=1 & i%%CUTOFF==1){j=j+1}
        ########################
        if(CUTOFF==1){j=j+1}  
        ########################
        i=i+1}
    print('Group Number:')
    print(j-1)
    return(GROUP)
}


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



.getValidpair <- function(DATA1, GROUP1, DATA2, GROUP2, CPU=4, method='kendall', print_step=10){
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print_step=print_step
    method=method
    CPU=CPU
    print('Start')
    print('Step1.Generate Reference...')
    REF1=.generate_ref(DATA1, cbind(GROUP1, GROUP1), min_cell=1) 
    REF2=.generate_ref(DATA2, cbind(GROUP2, GROUP2), min_cell=1) 
    print('Step2.Calculate Correlation Coefficient...')
    out = .get_cor( REF1, REF2, method=method,CPU=CPU, print_step=print_step)
    print('Step3.Analyze Result...')
    tag1=.get_tag_max(out)
    tag2=.get_tag_max(t(out))
    V=c()
    i=1
    while(i<=nrow(tag1)){
        t1=tag1[i,1]
        t2=tag1[i,2]
        if(tag2[which(tag2[,1]==t2),2]==t1){V=c(V,i)}           
        i=i+1}
    VP=tag1[V,]
    ##############################
    if(length(V)<=1){return(message("If BEER crashed, please try a different GNUM to get valid pair."))}
    ##############################
    C=c()
    t=1
    while(t<=nrow(VP)){
        this_c=out[which(rownames(out)==VP[t,2]),which(colnames(out)==VP[t,1])]
        C=c(C,this_c)
        t=t+1}
    #if(do.plot==TRUE){plot(C)}
    #VP=VP[which(C>=CUTOFF),]  
    print('Finished!!!')
    OUT=list()
    OUT$vp=VP
    OUT$cor=C
    return(OUT)
    }




.getVPnet<- function(pbmc, ROUND){
    library(igraph)
    pbmc=pbmc
    ROUND=ROUND

    ################
    REF=.generate_ref(pbmc@assays$RNA@data, cbind(pbmc@meta.data$group,pbmc@meta.data$group),min_cell=1)
    VREF=REF
    CVREF=cor(VREF,method='spearman')
    #ROUND=3
    VP=c()
    I=1    
    while(I<=ROUND){
        
        p1=c()
        p2=c()
        score=c()
        i=1
        while(i<=nrow(CVREF)){
            this_p1=rownames(CVREF)[i]
            j=i+1
            while(j<=ncol(CVREF)){
                this_p2=colnames(CVREF)[j]  
                p1=c(p1,this_p1)
                p2=c(p2,this_p2)
            
                vp_index=which(VP[1,]==this_p1 & VP[2,]==this_p2) 
                if(length(vp_index) >0){
                    this_score=999999 }else{
                    this_score=1-CVREF[i,j]}         
                score=c(score,this_score)
                j=j+1}        
            i=i+1}


        NET = cbind(p1,p2) 
        g <- make_graph(t(NET),directed = FALSE)
        MST=mst(g, weights = score, algorithm = NULL)
        E_MST=as_edgelist(MST, names = TRUE)

        i=1
        while(i<=nrow(E_MST)){
            t1=unlist(strsplit(E_MST[i,1],'_'))[1]
            t2=unlist(strsplit(E_MST[i,2],'_'))[1]
            if(t1!=t2){VP=cbind(VP,E_MST[i,])}
                i=i+1}

        print(I)
        I=I+1
        }       
        
        
    VP=t(VP)
    #######################
    return(VP)
    }





.getVPall<- function(pbmc, ROUND){
    
    pbmc=pbmc
    ROUND=ROUND

    ################
    REF=.generate_ref(pbmc@assays$RNA@data, cbind(pbmc@meta.data$group,pbmc@meta.data$group),min_cell=1)
    VREF=REF
    CVREF=cor(VREF,method='spearman')
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
                if(this_cor==max(this_cor_mat[i,]) & this_cor==max(this_cor_mat[,j])){                 
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
                CVREF[which(rownames(CVREF)==p1), which(rownames(CVREF)==p2)]=-99999
                CVREF[which(rownames(CVREF)==p2), which(rownames(CVREF)==p1)]=-99999    
                #b1=.get_batch(p1)
                #b2=.get_batch(p2)
                #b1i=which(group_batch==b1)
                #b2i=which(group_batch==b2)  
                #CVREF[b1i,b2i][which(rownames(CVREF[b1i,b2i])==p1),]= -99999
                #CVREF[b1i,b2i][,which(colnames(CVREF[b1i,b2i])==p2)]=   -99999                 
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
            #print(dim(this_cor_mat))
            #print(b1)
            #print(b2)
                j=j+1}
            i=i+1
            }
        print(I)
        I=I+1
        }       
        
        
    VP=t(VP)
    #######################
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
        
        this_test=cor.test(lst1_quantile, lst2_quantile, method='kendall')        
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


BEER <- function(DATA, BATCH,  GNUM=30, PCNUM=50, GN=2000, CPU=4, MTTAG="^MT-", REGBATCH=FALSE, print_step=10, SEED=123, N=2, ROUND=3){

    set.seed( SEED)
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=DATA
    BATCH=BATCH
    GNUM=GNUM
    PCNUM=PCNUM
    MTTAG=MTTAG
    #MAXBATCH=MAXBATCH
    UBATCH=unique(BATCH)
    REGBATCH=REGBATCH
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
    for(this_batch in UBATCH){
        this_pbmc=CreateSeuratObject(counts = DATA[,which(BATCH==this_batch)], min.cells = 0, 
                                 min.features = 0, project = this_batch)
        this_pbmc <- NormalizeData(object = this_pbmc, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
        this_pbmc <- FindVariableFeatures(object = this_pbmc, selection.method = "vst", nfeatures = GN)  
        this_varg=VariableFeatures(object = this_pbmc)
        VARG=c(VARG, this_varg)
        }
    VARG=unique(VARG)


    pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL") 
    pbmc@meta.data$batch=BATCH
    VariableFeatures(object = pbmc)=VARG


    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = MTTAG)

    CPU=4
    if(REGBATCH==FALSE){
    pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("nCount_RNA","percent.mt"), num.cores=CPU, do.par=TRUE)
    }else{
    pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("nCount_RNA", "batch", "percent.mt"), num.cores=CPU, do.par=TRUE)
    }
    
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)


    DR=pbmc@reductions$umap@cell.embeddings
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=DR[this_index,]
        #CNUM=max(c(5, round(length(this_index)/GNUM) ))
        this_group=.getGroup(this_one,this_batch,GNUM)
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


.selectUSE <-function(RESULT, CUTR=0.7, CUTL=0.7, RR=0.5, RL=0.5, CC=0.05){
    PCUSE=which( (rank(RESULT$cor)>=length(RESULT$cor)*RR | RESULT$cor>CUTR )    & 
                (rank(RESULT$lcor) >=length(RESULT$cor)*RR | RESULT$lcor>CUTL)   &
                 p.adjust(RESULT$lc1,method='fdr') > CC
               ) 
    return(PCUSE)
    }

####################



ReBEER <- function(mybeer,  GNUM=30, PCNUM=50, GN=2000, CPU=4, MTTAG="^MT-", print_step=10, SEED=123, N=2, ROUND=3){

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
    MTTAG=MTTAG
    #MAXBATCH=MAXBATCH
    UBATCH=unique(BATCH)
    GN=GN
    ROUND=ROUND
    N=N
    print_step=print_step
    
    pbmc=mybeer$seurat
    
    print('Group number (GNUM) is:')
    print(GNUM)
    print('Total varible gene number (GN) is:')
    print(length(VariableFeatures(object = pbmc)))
    print('ROUND is:')
    print(ROUND)
    
    VARG=VariableFeatures(object = pbmc)
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    
    
    
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)
    ########
    
    DR=pbmc@reductions$umap@cell.embeddings
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=DR[this_index,]
        #CNUM=max(c(5, round(length(this_index)/GNUM) ))
        this_group=.getGroup(this_one,this_batch,GNUM)
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
    
    #########

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


BEER.bbknn <- function(pbmc, PCUSE, NB=3, NT=10){
  
    NB=NB
    NT=NT
    #mybeer=mybeer

    pbmc=pbmc#mybeer$seurat
    PCUSE=PCUSE
    
    #pbmc=mybeer$seurat
    batch=as.character(pbmc@meta.data$batch)
    
    pca.all=pbmc@reductions$pca@cell.embeddings
    pca.use=pbmc@reductions$pca@cell.embeddings[,PCUSE]
    
        
    library(reticulate)
    #use_python("C:\Users\cchmc\Anaconda3\python")
    
    
    anndata = import("anndata",convert=FALSE)
    bbknn = import("bbknn", convert=FALSE)
    sc = import("scanpy.api",convert=FALSE)

    adata = anndata$AnnData(X=pca.all, obs=batch)
    PCNUM=ncol(pca.use)

    sc$tl$pca(adata, n_comps=as.integer(PCNUM))
    adata$obsm$X_pca = pca.use
    #NB=50
    #NT=10
    #print(NB)
    bbknn$bbknn(adata,batch_key=0,neighbors_within_batch=as.integer(NB),n_pcs=as.integer(PCNUM), n_trees =as.integer(NT))
    #bbknn$bbknn(adata,batch_key=0, n_pcs=as.integer(PCNUM))
    
    sc$tl$umap(adata)
    umap = py_to_r(adata$obsm$X_umap)
    rownames(umap)=rownames(pbmc@reductions$umap@cell.embeddings)
    colnames(umap)=colnames(pbmc@reductions$umap@cell.embeddings)
  
    return(umap)
    }




