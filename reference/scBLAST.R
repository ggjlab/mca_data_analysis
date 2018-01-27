
# get_tissue_arv 
# this function return a single subs gene all-cell mean and 3 sample 100cells
get_tissue_arv <- function(cell_counts, subcount){
  single_tissue_names_per_subs = colnames(pbmc@data[,pbmc@ident==subcount])
  dge <- data.frame(as.matrix(pbmc@raw.data[,single_tissue_names_per_subs]))
  single_tissue_raw_data_per_subs <- as.matrix(t(t(dge)/colSums(dge))*1000000) 
  if(cell_counts > 100)	{
    single_sample_in_ref<-data.frame(matrix(nrow=nrow(single_tissue_raw_data_per_subs),ncol = 3))
    rownames(single_sample_in_ref)<-rownames(single_tissue_raw_data_per_subs)
         for (i in 1:3){
    single_sample_in_ref[,i] = rowMeans(single_tissue_raw_data_per_subs[,sample(1:length(single_tissue_raw_data_per_subs[1,]),100)])
              }             
    }
  else {single_sample_in_ref<-data.frame(matrix(nrow=nrow(single_tissue_raw_data_per_subs),ncol = 3))
  rownames(single_sample_in_ref)<-rownames(single_tissue_raw_data_per_subs)
  for (i in 1:3){
    single_sample_in_ref[,i] = rowMeans(single_tissue_raw_data_per_subs[,sample(1:length(single_tissue_raw_data_per_subs[1,]),length(single_tissue_raw_data_per_subs[1,])/2)])}
  }
  single_sample_in_ref<- floor(single_sample_in_ref/10)
  sample_tissue_mean = data.frame(   floor(apply(single_tissue_raw_data_per_subs,1   ,  mean)/10 )  )
  sumgene_sample_tissue <- sum(sample_tissue_mean>0)
  return(list(as.matrix(sample_tissue_mean),as.matrix(single_sample_in_ref),as.numeric(sumgene_sample_tissue)))   }
#get_tissue_sample_data
# get_3 sample 100cells data
get_tissue_sample_data <- function(pbmc){
  tissue_data=c()
  tissue_sample = c()
  tissue_sumgene =c()
  subs_count = length(table(pbmc@ident))
  cells_count_persub = as.numeric(table(pbmc@ident))
  for(i in 1:subs_count){
    xx=get_tissue_arv(cells_count_persub[i],i)	
    colnames(xx[[2]])<-paste0(i,"_",colnames(xx[[2]]))
    tissue_sample <- cbind(tissue_sample,xx[[1]])
    tissue_data <- cbind(tissue_data,xx[[2]])	
    tissue_sumgene <-c(tissue_sumgene,xx[[3]])
  }
  colnames(tissue_sample)<-1:length(cells_count_persub)
  rownames(tissue_sample)<-rownames(tissue_data)
  return(list(tissue_data,tissue_sample,tissue_sumgene)) }

####the seurat data including your clustering information is saved in the Rdata, you should set the pathway of Rdata as working pathway
#how to build a seurat data? please read the instructions of Seurat R package
setwd("/home/ggj/Documents/single_data/MCA_3.0/tissue_Rdata/")
tissuedata <- list.files(pattern="*.RData")
tissuenames <- reshape2::colsplit(tissuedata,pattern="[.]",names=c("tissue","c"))$tissue
total_tissue_gene = data.frame()
total_tissue_data = data.frame()
for(i in 1:length(tissuenames)){
  message(paste0("Loading ",tissuedata[i]))
  load(tissuedata[i])
  message("Finish Loading")
  xx = get_tissue_sample_data(pbmc)
  colnames(xx[[1]]) <- paste0(tissuenames[i],"_",colnames(xx[[1]]))
  colnames(xx[[2]]) <- paste0(tissuenames[i],"_",colnames(xx[[2]]))
  message("Staring mearge")
  if(i==1){total_tissue_gene = xx[[1]]; total_tissue_data = xx[[2]];tissue_sumgene <-xx[[3]]}
  else{ total_tissue_gene=merge(total_tissue_gene,xx[[1]],by="row.names",all=T,sort=T); rownames(total_tissue_gene)=total_tissue_gene[,1] ; total_tissue_gene=total_tissue_gene[,-1]; total_tissue_gene[is.na(total_tissue_gene)]=0 ; 
  total_tissue_data=merge(total_tissue_data,xx[[2]],by="row.names",all=T,sort=T); rownames(total_tissue_data)=total_tissue_data[,1] ; total_tissue_data=total_tissue_data[,-1]; total_tissue_data[is.na(total_tissue_data)]=0 ;
  tissue_sumgene <-c(tissue_sumgene,xx[[3]])
  }
  rm(pbmc)	
  message(paste0("Finish ",tissuenames[i]))	  }
tissue_sumgene<-data.frame(tissue_sumgene,row.names = colnames(tissuedata))
hist(tissue_sumgene$tissue_sumgene,breaks = 100,main="floor-ngene of celltype",xlab = "ngene")
summary(tissue_sumgene$tissue_sumgene)
#############################  
all_tissue_gene <- colnames(total_tissue_gene)
subs_group <- colnames(total_tissue_data)
names_id <- c()
#for(i in 1:length(subs_group)) { id=grep(pattern=paste0("^",subs_group[i],"_"),all_tissue_gene); ref.expr[,i]=rowMeans(total_tissue_gene_log[,id]); names_id = c(names_id,rep(subs_group[i],length(id)))  }
for(i in 1:length(subs_group)) { id=grep(pattern=paste0("^",subs_group[i],"_"),all_tissue_gene); names_id = c(names_id,rep(subs_group[i],length(id)))  }
ID_911<-data.frame(colnames(total_tissue_gene),names_id)
rownames(ID_911)<-ID_911$colnames.total_tissue_gene.
############################# run Seurat for the next step to calculate different gene test 
library(Seurat)
library(dplyr)
library(Matrix)
pbmc <- CreateSeuratObject(raw.data =total_tissue_gene, min.cells = 3, min.genes = 20, 
                           project = "nolog")
pbmc <- AddMetaData(object = pbmc, metadata =ID_911)
#pbmc<-StashIdent(pbmc,save.name = "Ident")
pbmc<-SetAllIdent(pbmc,id="names_id")
pbmc@ident
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                      scale.factor = 100000)
pbmc.markers <- FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 1,logfc.threshold = 1)
pbmc.markers %>% group_by(cluster) %>% top_n(10,avg_diff) -> top10
top10 <- top10$gene
top10 <- top10[!duplicated(top10)]
ref_top10<-total_tissue_data[top10,]
ref<-log(ref_top10+1)# finish build the reference 


#start scbalst
load("/home/ggj/Rdata/newFigure6/cds.RData")#load the matirx of new cells,row is genename ,col is cellname
cds<-log(cds+1)
tst <- data.frame(matrix(nrow =length(ref[,1]),ncol = length(cds[1,])))#the 104 is the test cell numbers 
rownames(tst)<-rownames(ref)
colnames(tst)<-colnames(cds)
for (i in rownames(ref)) {ref
  if(i%in%rownames(cds)) tst[i,]<- cds[i,]
}
tst[is.na(tst)]<-0
cors <- cor(ref,tst)
cor1<-cors
cor1m<-apply(cor1,2,max) #
cor1S<-apply(cors,2,function(x) rownames(cors)[which.max(x)])
cor1r<-cbind(cor1S,cor1m)
#scblast results is in scblast dataframe
scblast<-data.frame(cor1r)
colnames(scblast)<-c("scblast_result" ,  "cors_log" )







