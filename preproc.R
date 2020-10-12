#!/usr/bin/env Rscript 

library(pagoda2)
library(magrittr)
library(parallel)
library(Matrix)
require(rjson)
require(data.table)

args <- commandArgs(trailingOnly = TRUE)

args <- c("/home/rstudio/programming/gi_mok")

setwd("/home/rstudio/programming/gi_mok/outputs")



cat("processing 10x text format matrices from", args[1],'\n')



cd <- read.10x.matrices(args[1])
colnames(cd) <- gsub("^one_","",colnames(cd)) # strip out the prefix - don't need it since we're processing one dataset
# read gene IDs separately, since you want them
gene.map <- read.csv(paste(args[1],"genes.tsv",sep='/'),header=F,sep='\t',stringsAsFactors = F)
colnames(gene.map) <- c("id","name")
rownames(gene.map) <- gene.map[,1]
rownames(cd) <- gene.map[,2] # will use human readable names; if you want to use IDS, switch to [,1]


# Filter cells
cat("filtering ...")
cdf <- gene.vs.molecule.cell.filter(cd, min.cell.size = 1e3,plot=F)
cat("done\n")

## Filter for doublets

get.scrublet.scores <- function(mat,min.molecules.per.gene=10) {
  require(data.table)
  # write out a CSV file
  tf <- tempfile()
  dtf <- paste(tf,'doubletScores',sep='.')
  dt <- data.table(as.matrix(t(mat[Matrix::rowSums(mat)>=min.molecules.per.gene,])))
  data.table::fwrite(dt,file=tf)
  cmd <- paste("python -c 'import sys; import pandas; import scrublet; df = pandas.read_csv(\"",tf,"\"); scrub = scrublet.Scrublet(df); doublet_scores, predicted_doublets = scrub.scrub_doublets(); pandas.DataFrame(doublet_scores).to_csv(\"",dtf,"\");'",sep='')
  tmp <- system(cmd, intern=T)
  #system(cmd);
  x <- as.data.frame(data.table::fread(dtf,sep=',',header=F,skip=1))
  x <- as.numeric(as.data.frame(data.table::fread(dtf,sep=',',header=F,skip=1))[,2])
  names(x)<-colnames(mat)
  file.remove(tf);
  file.remove(dtf);
  x
}

#cat("scrublet scores ...")
#doublet.f <- get.scrublet.scores(cdf[Matrix::rowSums(cdf>0)>5,])
#cat("done\n")





cat("processing ...")

# Normalize

# gene names have to be unique
cd <- cdf[Matrix::rowSums(cdf>0)>5, ]; 
#cd <- cdf[,doublet.f[colnames(cdf)]<=0.2]  # doublet filtering
rownames(cd) <- make.unique(rownames(cd))
r <- Pagoda2$new(cd)


#Variance normalziation

r$adjustVariance(plot=F,gam.k=10)

# PCA
#Calculate PCs ... note we're taking 3000 top overdispersed genes (a lot more than the number that is significantly overdispersed)
r$calculatePcaReduction(nPcs=50,n.odgenes=3e3)


#kNN

r$makeKnnGraph(k=20,type='PCA',center=T,distance='cosine');

#get clusters
r$getKnnClusters(method=conos:::leiden.community,type='PCA')

#Embeddings
r$getEmbedding(type='PCA',embeddingType='tSNE')
cat("done\n")


#Output
cat("writing out files ...")
# p2 object
saveRDS(r,file='pagoda2.rds',compress=FALSE)

# Raw and scaled matrices
# Matrices


Matrix::writeMM(r$counts,file="normalized.mtx")
Matrix::writeMM(r$misc$rawCounts[rownames(r$counts),colnames(r$counts)],file="raw.mtx")
gdf <- gene.map[colnames(r$counts,),]
gdf$var <- r$misc$varinfo[as.character(gdf[,1]),'qv']
write.table(gdf,file='genes_output.tsv',quote=F,col.names=F,row.names=F,sep="\t")
write.table(rownames(r$counts),file='cells_output.tsv',quote=F,col.names=F,row.names=F,sep='\t')

# PCA 
## scores
write.table(r$reductions$PCA[rownames(r$counts),],file='X_pca_output.tsv',quote=F,row.names=F,col.names=F,sep='\t')
## loadings
write.table(r$misc$PCA$v,file='PCs_output.tsv',quote=F,row.names=T,col.names=F,sep='\t')

# embedding
write.table(r$embeddings$PCA$tSNE[rownames(r$counts),],file='embedding_output.tsv',quote=F,row.names=T,col.names=F,sep='\t')

# doublet scores
#df <- data.frame(cell=names(doublet.f),doubletScore=doublet.f)
#write.table(df,file='doublet.scores.tsv',col.names=F,row.names=F,quote=F,sep='\t')

cat("done\n")

