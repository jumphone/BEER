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

################################################################################################



mybeer5 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=5)
mybeer50 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=50)
mybeer100 <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2,PP=100)
saveRDS(mybeer5, 'mybeer5.RDS')
saveRDS(mybeer50, 'mybeer50.RDS')
saveRDS(mybeer100, 'mybeer100.RDS')


