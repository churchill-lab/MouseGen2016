## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(warning=FALSE, message=FALSE)

## ----R_package, results="hide"-------------------------------------------
library(DESeq2)
library(ggplot2)
library(dplyr)

## ----load_Robj, results="hide"-------------------------------------------
load("data/DO192_DataforSysGenCourse.Rdata")
exp.all = read.table("data/expected_read_counts_gene_level.txt", header=T)

## ----exp_data, results="hide"--------------------------------------------
geneIDs = exp.all[,1]
exp.all=exp.all[,-1]
rownames(exp.all)=geneIDs
exp.all[1:5,1:5]

## ----exp_design, results="hide"------------------------------------------
exp_design = data.frame(mouseIDs=colnames(exp.all),
                        diet=covariates.rna.192$Diet,
                        sex=covariates.rna.192$Sex,
                        coat_color=covariates.rna.192$Coat.Color)

## ----check_data, results="hide"------------------------------------------
all(colnames(exp.all)==exp_design$mouseIDs)

## ----check_xist, results="hide"------------------------------------------
geneID="ENSMUSG00000086503"
geneName="Xist"
gIndex = which(rownames(exp.all)==geneID)
data= data.frame(exp_design, 
                 exp=as.numeric(exp.all[gIndex,]))

## ----head_exp_design-----------------------------------------------------
head(data)


## ----plot_xist, results="hide"-------------------------------------------
p <- ggplot(data,aes(x=sex,y=exp)) 
p <- p + geom_point(position = position_jitter(width = 0.2),size=3,
                    aes(colour = factor(sex)))
p <- p + stat_summary(fun.y=mean, geom="point", shape=5, size=4)
p <- p + ylab("Gene Expression (Read Counts)")
p <- p + theme(axis.text=element_text(size=12),
               axis.title=element_text(size=12,face="bold", colour = "blue"),
               plot.title = element_text(size = rel(2)))
p <- p + ggtitle("Xist: ENSMUSG00000086503")
p

## ----sample_ind, results="hide"------------------------------------------
male_index = which(exp_design$sex=="M")
female_index = which(exp_design$sex=="F")
chow_index= which(exp_design$diet=="chow")
hf_index= which(exp_design$diet=="HF")
male_chow = intersect(male_index,chow_index)
male_hf = intersect(male_index,hf_index)

## ----sample_size, results="hide"-----------------------------------------
sampleSize = 3

## ----subset_exp_3--------------------------------------------------------
diet_DE = c(male_chow[1:sampleSize],male_hf[1:sampleSize])
exp_design_diet_DE= exp_design[diet_DE,]
head(exp_design_diet_DE)
exp_diet_DE=exp.all[,diet_DE]
all(colnames(exp_diet_DE)==as.vector(exp_design_diet_DE$mouseIDs))

## ----head_exp------------------------------------------------------------
head(exp_diet_DE)

## ----filter_exp, results="hide"------------------------------------------
thres= 5
nzIndex= as.vector(which(apply(exp_diet_DE,1,function(x){sum(x>thres)/length(x)})>=0.5))
head(nzIndex)
exp.dietDE = exp_diet_DE[nzIndex,]
dim(exp.dietDE)

## ----dataframe_deseq2, results="hide"------------------------------------
### colData contains the condition/group information for Differenetial expression analysis
colData <- DataFrame(group = factor(exp_design_diet_DE$diet))

## ----deseq2_obj, results="hide"------------------------------------------
### Create DESeq2 object using expression and colData
dds_3reps <- DESeqDataSetFromMatrix(countData = as.data.frame(round(exp.dietDE)),
         colData = colData, design = ~ group)
dds_3reps <- DESeq(dds_3reps)
res_3reps = results(dds_3reps)
resOrdered_3reps <- res_3reps[order(res_3reps$padj),]
head(resOrdered_3reps)

## ----res_summary---------------------------------------------------------
### summary of Differential Expression analysis
summary(res_3reps)

## ----pval_hist, results="hide"-------------------------------------------
hist(res_3reps$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value histogram: 3 Samples per group")

## ----r_helper_functions--------------------------------------------------
### helper functions to get gene information for a gene
get_gene_name<-function(ensID,gene_df){
  g_ind = which(as.vector(gene_df[,1])==ensID)
  g_name = gene_df[g_ind,2]
  return(g_name)
}
get_gene_info<-function(ensID,gene_df){
  g_ind = which(as.vector(gene_df[,1])==ensID)
  print(g_ind)
  ensID= as.vector(gene_df[g_ind,1])
  gName= as.vector(gene_df[g_ind,2])
  chro = as.vector(gene_df[g_ind,3])
  start = as.vector(gene_df[g_ind,5])
  end = as.vector(gene_df[g_ind,6])
  g_name = paste0(ensID,": ",gName," Chr",chro,"-",start,":",end)
  g_name = paste0(ensID,": ",gName," Chr",chro)
  return(g_name)
}

## ----top_gene_plots------------------------------------------------------
#par(mfrow=c(2,3),las=1)
top5_genes= rownames(resOrdered_3reps[1:5,])
gene_info=read.csv(url("ftp://ftp.jax.org/dgatti/ShortCourse2015/ENSMUSG-gene-info-R84.tsv"),header=FALSE,sep="\t")
print(head(gene_info))
for (i in 1:length(top5_genes)){
  g_id = top5_genes[i]
  g_name = get_gene_name(g_id,gene_info)
  g_info = get_gene_info(g_id,gene_info)
  data <- plotCounts(dds_3reps, gene=g_id, intgroup=c("group"), returnData=TRUE)
  p <- ggplot(data, aes(x=group, y=count, color=group))
  p <- p+ ggtitle(g_info)
  p <- p+ geom_point(position=position_jitter(width=.1,height=0), size=3)
  p <- p + theme(axis.text=element_text(size=12),axis.title=element_text(size=20,face="bold", colour = "blue"),  plot.title = element_text(size =rel(1.5)))
  print(p)
}

## ----sample_size_10, results="hide"--------------------------------------
sampleSize = 10
diet_DE = c(male_chow[1:sampleSize],male_hf[1:sampleSize])
exp_design_diet_DE= exp_design[diet_DE,]
head(exp_design_diet_DE)
exp_diet_DE=exp.all[,diet_DE]
all(colnames(exp_diet_DE)==as.vector(exp_design_diet_DE$mouseIDs))
head(exp_diet_DE)

## ----filter_exp_10, results="hide"---------------------------------------
thres= 5
nzIndex= as.vector(which(apply(exp_diet_DE,1,function(x){sum(x>thres)/length(x)})>=0.5))
head(nzIndex)
exp.dietDE = exp_diet_DE[nzIndex,]
dim(exp.dietDE)

## ----dataframe_deseq2_10, results="hide"---------------------------------
### colData contains the condition/group information for Differenetial expression analysis
colData <- DataFrame(group = factor(exp_design_diet_DE$diet))

## ----deseq2_obj_10, results="hide"---------------------------------------
### Create DESeq2 object using expression and colData
dds_10reps <- DESeqDataSetFromMatrix(countData = as.data.frame(round(exp.dietDE)),
         colData = colData, design = ~ group)
dds_10reps <- DESeq(dds_10reps)
res_10reps = results(dds_10reps)
resOrdered_10reps <- res_10reps[order(res_10reps$padj),]
head(resOrdered_10reps)


## ----res_summary_10------------------------------------------------------
### summary of Differential Expression analysis
summary(res_10reps)
summary(res_3reps)


## ----pval_hist_10, fig.width=12, fig.height=8----------------------------
par(mfrow=c(1,2))
hist(res_10reps$pvalue,breaks=100,col="grey", xlab="p-value",main="p-value Histogram: 10 Samples per group")
hist(res_3reps$pvalue,breaks=100,ylim=c(0,1200),col="grey", xlab="p-value",main="p-value histogram: 3 Samples per group")


## ----MA_plot_10----------------------------------------------------------
plotMA(res_10reps, main="M-A Plot: 10 Samples per group", ylim=c(-2,2))

## ----top_gene_plots_10---------------------------------------------------
top5_genes= rownames(resOrdered_10reps[1:5,])
par(mfrow=c(2,3),las=1)
for (i in 1:length(top5_genes)){
  g_id = top5_genes[i]
  g_name = get_gene_name(g_id,gene_info)
  g_info = get_gene_info(g_id,gene_info)
  data <- plotCounts(dds_10reps, gene=g_id, intgroup=c("group"), returnData=TRUE)
  p <- ggplot(data, aes(x=group, y=count, color=group))
  p <- p+ ggtitle(g_info)
  p <- p+ geom_point(position=position_jitter(width=.1,height=0), size=3)
  p <- p + theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold", colour = "blue"),  plot.title = element_text(size =rel(1.5)))
  print(p)
}

