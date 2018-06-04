load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
require(pheatmap)
require(ggplot2)
require(quantreg)

MET500.tumor.content <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/MET500.tumor.content.csv", stringsAsFactors=FALSE)
flag                 <- MET500.tumor.content$sequecing.tumor.content == "<30%"
MET500.tumor.content[flag,'sequecing.tumor.content'] <- '-1%'
MET500.tumor.content$sequecing.tumor.content <- sub(x = MET500.tumor.content$sequecing.tumor.content,replacement = '',pattern = '%') %>% as.numeric
rownames(MET500.tumor.content) <- MET500.tumor.content$MET.500.id

###### Select marker genes ################## 
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
marker.gene <- names(sort(rank.sd,decreasing =TRUE))[1:1000]


####### Compute the correlation between MET500 samples and CCLE samples#######
common.gene          <- intersect(marker.gene,rownames(MET500.log2.fpkm.matrix))
cor.matrix           <- cor(MET500.log2.fpkm.matrix[common.gene,],CCLE.log2.read.count.matrix[common.gene,],method='spearman')

prostate.cancer.sample <- rownames(MET500.sample.meta)[MET500.sample.meta$cancer.type == 'Prostate Adenocarcinoma' | MET500.sample.meta$cancer.type == 'Prostate Neuroendocrine Carcinoma'] %>% as.character
prostate.cor.matrix    <- cor.matrix[prostate.cancer.sample,]
cell.line.median.cor   <- apply(prostate.cor.matrix,2,median)  %>% sort(decreasing = FALSE)
flag                   <- CCLE.sample.meta[names(cell.line.median.cor),'cell.line.tumor.site'] == 'PROSTATE'
prostate.cell.line     <- names(cell.line.median.cor)[flag] %>% as.character

plot(cell.line.median.cor,col=ifelse(flag == TRUE,'red','grey'),pch=19,cex=0.5,xlab='rank',ylab='median.corration.with.MET500',cex.lab=1.5)
points(rank(cell.line.median.cor)[flag],cell.line.median.cor[flag],col='red',pch=19,cex=1.5)
boxplot(prostate.cor.matrix[,prostate.cell.line],cex.lab=3.5)

good.prostate.cell.line <- c('MDAPCA2B_PROSTATE','VCAP_PROSTATE','LNCAPCLONEFGC_PROSTATE')
data          <- prostate.cor.matrix[,good.prostate.cell.line]
pheatmap(data,show_rownames = F,show_colnames = T,cluster_cols = T)

biopsy.site <- MET500.sample.meta[ rownames(data),'biopsy.site']
#pheatmap(data[biopsy.site == 'LIVER',],show_rownames = F,show_colnames = T,cluster_cols = T)
#pheatmap(data[biopsy.site == 'BONE',],show_rownames = F,show_colnames = T,cluster_cols = T)
#pheatmap(data[biopsy.site == 'LYMPH_NODE',],show_rownames = F,show_colnames = T,cluster_cols = T)


pca.rs         <- prcomp(data)
plot(pca.rs$x[,1:2])
pc1 <- pca.rs$x[,1] 
pc2 <- pca.rs$x[,2] 
draw.df <- data.frame(PC1=pc1,PC2=pc2,biopsy.site = MET500.sample.meta[names(pc1) %>% as.character,'biopsy.site'],mean.cor=apply(data,1,mean),library.selection=MET500.sample.meta[names(pc1) %>% as.character,'LibrarySelection'])

prostate.cancer.sample.MET500.id <- MET500.sample.meta[rownames(draw.df) %>% as.character,'MET500.id']
tumor.content                    <- MET500.tumor.content[prostate.cancer.sample.MET500.id,'sequecing.tumor.content']
draw.df$tumor.content            <- tumor.content


ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE','BONE'),]) + geom_point(aes(x=PC1,y=PC2,color=biopsy.site),size=3)
ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE','BONE'),]) + geom_boxplot(aes(y=mean.cor,x=biopsy.site,color=biopsy.site),size=3) + theme(axis.text=element_text(size=21),axis.title=element_text(size=21),legend.text=element_text(size=21))
ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE','BONE'),]) + geom_boxplot(aes(y=tumor.content,x=biopsy.site,color=biopsy.site),size=3)+ theme(axis.text=element_text(size=21),axis.title=element_text(size=21),legend.text=element_text(size=21))

t.test(draw.df$mean.cor[draw.df$biopsy.site == 'BONE'],draw.df$mean.cor[draw.df$biopsy.site == 'LYMPH_NODE'])
wilcox.test(draw.df$mean.cor[draw.df$biopsy.site == 'BONE'],draw.df$mean.cor[draw.df$biopsy.site == 'LYMPH_NODE'])

t.test(draw.df$mean.cor[draw.df$biopsy.site == 'BONE'],draw.df$mean.cor[draw.df$biopsy.site == 'LIVER'])
wilcox.test(draw.df$mean.cor[draw.df$biopsy.site == 'BONE'],draw.df$mean.cor[draw.df$biopsy.site == 'LIVER'])

t.test(draw.df$mean.cor[draw.df$biopsy.site == 'LYMPH_NODE'],draw.df$mean.cor[draw.df$biopsy.site == 'LIVER'])
wilcox.test(draw.df$mean.cor[draw.df$biopsy.site == 'LYMPH_NODE'],draw.df$mean.cor[draw.df$biopsy.site == 'LIVER'])


fit.data <- draw.df[draw.df$tumor.content >0  ,]
line.fit <- rq(formula=mean.cor ~ tumor.content,data=fit.data,tau=0.5)
plot(x=fit.data$tumor.content,y=fit.data$mean.cor,xlim=c(0,100),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='mean.cor',cex.axis=1.5,pch=19)
lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)

plot(x=line.fit$fitted.values,y=line.fit$residuals,xlab='fitted.values',ylab='residual',cex.axis=1.5,pch=19,cex.lab=2)
abline(h=0)
plot(sort(line.fit$residuals),cex.lab=2,xlab='rank',ylab='residual',cex.axis=1.5,pch=19)
abline(h=quantile(line.fit$residuals)[2] - 1.5 * IQR(line.fit$residuals))

color <- ifelse(line.fit$residuals <= sort(line.fit$residuals)[8] ,'purple','black')
plot(x=fit.data$tumor.content,y=fit.data$mean.cor,xlim=c(0,100),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='mean.cor',cex.axis=1.5,pch=19,col=color)
lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)

fit.data <- draw.df[draw.df$tumor.content >=0  & draw.df$library.selection=='PolyA',]
line.fit <- rq(formula=mean.cor ~ tumor.content,data=fit.data,tau=0.5)
plot(x=fit.data$tumor.content,y=fit.data$mean.cor,,xlim=c(0,100),ylim=c(0,1))
lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=3)

fit.data <- draw.df[draw.df$tumor.content >=0  & draw.df$library.selection=='Hybrid Selection',]
line.fit <- rq(formula=mean.cor ~ tumor.content,data=fit.data,tau=0.5)
plot(x=fit.data$tumor.content,y=fit.data$mean.cor,,xlim=c(0,100),ylim=c(0,1))
lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=3)


draw.df <- draw.df[prostate.cancer.sample,]
draw.df$cancer.type <- MET500.RNASeq.meta[match( rownames(draw.df) %>% as.character, MET500.RNASeq.meta$Run %>% as.character),'cancer.type']
#### check expression of prostate marker genes
GTEX.tissue.specific.gene <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/GTEX.tissue.specific.gene.csv", stringsAsFactors=FALSE)
prostate.specific.gene    <- GTEX.tissue.specific.gene$gene_id[GTEX.tissue.specific.gene$Tissue == 'Prostate' ] %>% as.character
remove.ensemble.version.id <- function(x){
  v <- strsplit(x = x,split = "\\.") %>% unlist  
  v[1]
}
prostate.specific.gene <- sapply(prostate.specific.gene,remove.ensemble.version.id) %>% as.character
CHGA <- 'ENSG00000100604'
PSA  <- 'ENSG00000142515'

plot(x=MET500.log2.fpkm.matrix[CHGA,prostate.cancer.sample],y=MET500.log2.fpkm.matrix[PSA,prostate.cancer.sample],xlab='CHGA',ylab='PSA',cex.lab=1.5,cex.axis=1.5,pch=19)
hist(MET500.log2.fpkm.matrix[CHGA,prostate.cancer.sample],breaks=30)
plot(MET500.log2.fpkm.matrix[CHGA,prostate.cancer.sample] %>% sort)


draw.df$CHGA <- MET500.log2.fpkm.matrix[CHGA,prostate.cancer.sample]
draw.df$PSA  <- MET500.log2.fpkm.matrix[PSA,prostate.cancer.sample]
hehe <- draw.df
ggplot(hehe) + geom_point(aes(x=CHGA,y=PSA,color=cancer.type),size=4)

low.CHGA.pop  <- draw.df[draw.df$CHGA <=0,]
high.CHGA.pop <- draw.df[draw.df$CHGA >= 5,]

plot(low.CHGA.pop$tumor.content,y=low.CHGA.pop$mean.cor)
points(high.CHGA.pop$tumor.content,y=high.CHGA.pop$mean.cor,col='red')


data   <- CCLE.log2.rpkm.matrix[prostate.specific.gene,]
pca.rs <- prcomp(data %>% t)
plot(pca.rs$x[,1:2],cex.axis=1.5,pch=19,cex.lab=1.5)



hehe <- prostate.specific.gene[prostate.specific.gene %in% rownames(MET500.log2.fpkm.matrix)]
data <- MET500.log2.fpkm.matrix[hehe,prostate.cancer.sample] 
pca.rs <- prcomp(data %>% t)
plot(pca.rs$x[,1:2])






hybrid.lib.id <- rownames(MET500.sample.meta)[MET500.sample.meta$LibrarySelection == 'Hybrid Selection']
hybrid.lib.id <- hybrid.lib.id[hybrid.lib.id %in% prostate.cancer.sample]
h <- MET500.sample.meta[hybrid.lib.id %>% as.character,'MET500.id']

data <- MET500.log2.fpkm.matrix[,hybrid.lib.id] 
tt <- MET500.tumor.content[h %>% as.character,'sequecing.tumor.content']
data <- data[,tt>=70]

pca.rs <- prcomp(data %>% t)
plot(pca.rs$x[,1:2])
pc1 <- pca.rs$x[,1]

# CCLE.cell.line.idx     <- apply(prostate.cor.matrix,1,function(x) which(x == max(x))[1])
# CCLE.max.cor           <- apply(prostate.cor.matrix,1,function(x) max(x))
# rs.df                  <- cbind(sample.name = rownames(prostate.cor.matrix),CCLE.max.cor=CCLE.max.cor,CCLE.sample.meta[CCLE.cell.line.idx,])
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# PCA.explore <- function(expr.matrix,color.label,pc.idx=1:5,feature.number=1000,rank.based=TRUE){
#     rank.matrix <- apply(expr.matrix,2,rank)
#     if(rank.based == FALSE){
#         rank.matrix <- expr.matrix
#     }
#     rank.mean   <- apply(rank.matrix,1,mean)
#     rank.sd     <- apply(rank.matrix,1,sd)
#     plot(rank.mean,rank.sd)
#     lowess(rank.mean,rank.sd) %>% lines(col='red',lwd=5)
#     marker.gene <- names(sort(rank.sd,decreasing = TRUE))[1:1000]
#     pca.rs      <- prcomp(rank.matrix[marker.gene,] %>% t)
#     #df <- data.frame(pc1=pca.rs$x[,1],pc2=pca.rs$x[,2],batch.id=batch.id) 
#     #ggplot(df) + geom_point(aes(x=pc1,y=pc2,col=batch.id %>% factor))
#     col.vec <- brewer.pal(9,'Set1')
#     pairs(pca.rs$x[,pc.idx],col=col.vec[color.label %>% factor %>% as.integer],pch=19)
#     pca.rs
# }
# 
# ########
# 
# pca.rs <- PCA.explore(expr.matrix = MET500.log2.fpkm.matrix[,breast.cancer.sample],color.label = 'red',pc.idx = 1:2,feature.number = 1000,rank.based = TRUE)
# 
# 
# 
# 
# polyA.sample <-  rownames(MET500.sample.meta)[MET500.sample.meta$LibrarySelection == 'POLYA'] %>% as.character
# pca.rs <- PCA.explore(expr.matrix = MET500.log2.fpkm.matrix[,polyA.sample],color.label = 'red',pc.idx = 1:2,feature.number = 1000,rank.based = TRUE)
# 
# pca.rs <- prcomp(MET500.log2.fpkm.matrix[,polyA.sample] %>% t)
# pairs(pca.rs$x[,1:4])
# 
# #################
# load('server-side/RData/Skin Cutaneous Melanoma.RData')
# annotation.df          <- GDC.sample.meta[colnames(log2.read.count.matrix) %>% as.character,]
# flag                   <- annotation.df$sample.type == 'Primary Solid Tumor'
# annotation.df          <- annotation.df[flag,]
# log2.fpkm.matrix       <- log2.fpkm.matrix[,flag]
# log2.read.count.matrix <- log2.read.count.matrix[,flag]
# 
# pca.rs <- PCA.explore(expr.matrix = log2.fpkm.matrix,
#                       color.label = annotation.df$sample.type,
#                       pc.idx = 1:2,
#                       feature.number = 500,
#                       rank.based = TRUE)
# 
# tumor.purity           <- read.csv("~/Project/Cancer2CellLine/tumor.purity.csv", stringsAsFactors=FALSE)
# rownames(tumor.purity) <- tumor.purity$Sample.ID
# hcc.purity             <- tumor.purity[rownames(pca.rs$x) %>% as.character,]
# plot(x=pca.rs$x[,1],y=hcc.purity$CPE)
# pure.samples           <-rownames(hcc.purity) [hcc.purity$CPE >=0.9]
# hist(pca.rs$x[pure.samples %>% as.character,1])
# plot(pca.rs$x[pure.samples %>% as.character,1] %>% sort)
# 
# pc1 <- pca.rs$x[pure.samples,'PC1'] %>% sort
# 
# tail.sample <- tail(names(pc1),9)
# top.sample  <- head(names(pc1),9)
# 
# tmp  <- 2^log2.read.count.matrix -1
# tmp  <- tmp[,as.character(c(top.sample,tail.sample))] %>% as.integer
# flag <- apply(tmp,1,function(x) sum(x>=1) >= 6)
# tmp <- tmp[flag,]
# require(DESeq2)
# 
# df <- data.frame(#purity = hcc.purity[c(tail.sample,top.sample) %>% as.character,'CPE'],
#                  condition  = c(rep(x = 'tail',9),rep(x='top',9))
#                  )
# dds <- DESeqDataSetFromMatrix(countData = round(tmp),
#                               colData = df,
#                               design= ~ condition )
# dds <- DESeq(dds)
# resultsNames(dds) # lists the coefficients
# res <- results(dds, contrast=c("condition","tail","top"))
# res <- res[complete.cases(res),]
# res <- res[order(res$padj),]
# 
# de.gene <- rownames(res)[res$padj < 0.001 & abs(res$log2FoldChange) > 1]
# write.csv(x=de.gene,file='de.gene.csv',quote=F)
# ######
# load('server-side/RData/metastasis_rsem_FPKM.RData')
# 
# 
# rownames(exprs) <- sapply(rownames(exprs),remove.ensemble.version.id)
# 
# common.gene        <- intersect(marker.gene,rownames(exprs))
# cor.matrix         <- cor(exprs[common.gene,],CCLE.log2.read.count.matrix[common.gene,],method='spearman')
# CCLE.cell.line.idx <- apply(cor.matrix,1,function(x) which(x == max(x))[1])
# CCLE.max.cor       <- apply(cor.matrix,1,function(x) max(x))
# rs.df              <- cbind(sample.name = rownames(cor.matrix),CCLE.max.cor=CCLE.max.cor,CCLE.sample.meta[CCLE.cell.line.idx,])
# 
# rownames(rs.df) <- sapply(rownames(rs.df),function(x) strsplit(x,split = '_')[[1]][2])
# 
# rs.df$cancer.type <- SraRunTable[rownames(rs.df) %>% as.character,'histological_type']
# 
# #############
# 
# cancer.RData.file <- 'server-side/RData/Skin Cutaneous Melanoma.RData'
# load(cancer.RData.file)
# 
# annotation.df      <- GDC.sample.meta[colnames(log2.read.count.matrix) %>% as.character,]
# flag               <- annotation.df$sample.type == 'Primary Solid Tumor'
# annotation.df      <- annotation.df[flag,]
# log2.fpkm.matrix    <- log2.fpkm.matrix[,flag]
# log2.read.count.matrix <- log2.read.count.matrix[,flag]
# 
# 
# common.gene        <- intersect(marker.gene,rownames(log2.read.count.matrix))
# cor.matrix         <- cor(log2.read.count.matrix[common.gene,],CCLE.log2.read.count.matrix[common.gene,],method='spearman')
# CCLE.cell.line.idx <- apply(cor.matrix,1,function(x) which(x == max(x))[1])
# CCLE.max.cor       <- apply(cor.matrix,1,function(x) max(x))
# rs.df              <- cbind(sample.name = rownames(cor.matrix),CCLE.max.cor=CCLE.max.cor,CCLE.sample.meta[CCLE.cell.line.idx,])
# 
# cor.median         <- apply(cor.matrix,2,median ) 
# cell.line.idx      <- sort(cor.median) %>% names
# 
# ova.score           <- read.csv("~/Project/Cancer2CellLine/ova.score.csv", stringsAsFactors=FALSE)
# rownames(ova.score) <- paste(ova.score$cell.line,"_OVARY",sep="")
# 
# it <- intersect(rownames(ova.score),names(cor.median))
# plot(x=cor.median[it] %>% rank,
#      y=ova.score[it %>% as.character,'Score'] %>% rank,
#      xlab='rank of cor median',
#      ylab='rank of Suitability',
#      cex=1,pch=19,cex.lab=1.5)
# lines(c(1,40),c(1,40),lwd=3,col='red')
# 
# 
# 
# sample.batch <- sapply(annotation.df$sample.name,function(x) strsplit(x = x %>% as.character,split ='\\-')[[1]][2] )
# sample.gender <- ifelse(log2.read.count.matrix['ENSG00000012817',] > 1,'Male','Female')
# sample.antigen <- ifelse(log2.read.count.matrix['ENSG00000176566',] > 1,'Male','Female')
# sample.type <- annotation.df$sample.type
# 
# pca.rs <- PCA.explore(expr.matrix = log2.read.count.matrix,color.label = sample.type,pc.idx = 1:2)
# plot(x=pca.rs$x[,1],y=pca.rs$x[,2],xlab='PC1',ylab='PC2',col=c('red','black','purple')[sample.type %>% factor %>% as.integer],pch=19)
# 
# pca.rs <- explore.batch.effects(expr.matrix = log2.read.count.matrix,batch.id = annotation.df$sample.type)
# pca.rs <- explore.batch.effects(expr.matrix = log2.read.count.matrix,batch.id = sample.gender)
# pca.rs <- explore.batch.effects(expr.matrix = log2.read.count.matrix,batch.id = sample.type)
# 
# 
# 
# 
# plot(cor.median[cell.line.idx], col = ifelse(grepl(x=(cell.line.idx),pattern='OVARY'),'red','black'),pch=19,xlab='rank',ylab='median.of.correlation',cex.lab=2,cex.axis=2,cex=0.5)
# plot(cor.median[cell.line.idx], col = ifelse(grepl(x=(cell.line.idx),pattern='OVARY'),'red','white'),pch=19,xlab='rank',ylab='median.of.correlation',cex.lab=2,cex.axis=2,cex=0.5)
# 
# 
# 
# plot(cor.median[cell.line.idx], col = ifelse(grepl(x=(cell.line.idx),pattern='OVARY'),'red','grey'),pch=19,xlab='rank',ylab='median.of.correlation',cex.lab=2,cex.axis=2,cex=0.5)
# points(x=1:length(cor.median),y=cor.25.quantile[cell.line.idx],col = ifelse(grepl(x=(cell.line.idx),pattern='LIVER'),'red','black'),pch=19,xlab='rank',ylab='median.of.correlation')
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# boxplot(cor.matrix[,cell.line.idx],
#         col=ifelse( grepl(x=cell.line.idx,pattern='LIVER'),'purple','black')
# )
# boxplot(cor.matrix[,cell.line.idx],
#         col=ifelse( grepl(x=cell.line.idx,pattern='LIVER'),'red','white')
# )
# 
# 
# 
# m <- foreach(i = 1:nrow(cor.matrix),.combine='rbind') %do% {
#     row        <- cor.matrix[i,]
#     names(row) <- colnames(cor.matrix)
#     row        <- sort(row,decreasing = TRUE)
#     v          <- ifelse(grepl(x = names(row),pattern = 'CENTRAL_NERVOUS_SYSTEM
# '),1,0)
#     #v[grepl(x = names(row),pattern = 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE')] <- 1
#     v
# }
# rownames(m) <- rownames(cor.matrix)
# 
# 
# #pheatmap(m,cluster_rows = T,cluster_cols = F,show_rownames = F)
# annotation.row.df = data.frame(annotation.df[rownames(m) %>% as.character,'sample.type'])
# rownames(annotation.row.df) <- rownames(m)
# pheatmap(m,cluster_rows = T,cluster_cols = F,annotation_row = annotation.row.df,show_rownames = F)
# 
# 
# 
# 
# 
# pca.rs <- prcomp(cor.matrix)
# pairs(pca.rs$x[,1:8],col=c('red','purple','green','yellow','orange')[ sample.batch%>% factor %>% as.integer],pch=19)
# 
# 
# pca.rs <- prcomp(log2.fpkm.matrix %>% t)
# ggplot(data.frame(pc1=pca.rs$x[,1],pc2=pca.rs$x[,2],color=sample.batch))+geom_point(aes(x=pc1,y=pc2,color=color))
# pairs(pca.rs$x[,1:10],col=c('red','purple','green','yellow','orange','grey')[ sample.batch%>% factor %>% as.integer],pch=19)
# pairs(pca.rs$x[,1:8],col=c('red','purple','green','yellow','orange')[ annotation.df$sample.type %>% factor %>% as.integer],pch=19)
# 
# 
# 
# 
# 
# # 
# # r2 <- rs$rotation[,'PC2']
# # r2 <- sort(r2)
# # 
# # pc1 <- sort(rs$x[,1])
# # pheatmap(m[names(pc1),],cluster_rows = F,cluster_cols = F)
# # 
# # weired.sample.name <- tail(names(pc1),24)
# # weired.sample.name <- weired.sample.name[17:24]
# 
# ###require(pheatmap)
# cor.matrix <- cor((CCLE.log2.read.count.matrix[marker.gene,184:350]),method='spearman')
# dist.matrix <- 1 - cor.matrix
# 
# hclust(as.dist(dist.matrix)) %>% plot
# pheatmap(cor.matrix,cluster_rows = T,cluster_cols = T)
# 
# rs.1 <- prcomp(CCLE.log2.rpkm.matrix[marker.gene,137:183] %>% t)
# plot(rs.1$x[,1:2])
# 
# 
# rs.2 <- prcomp(CCLE.log2.rpkm.matrix[! (rownames(CCLE.log2.rpkm.matrix) %in% marker.gene),137:183] %>% t)
# plot(rs.2$x[,1:2])
# 
# 
# plot(rs.1$x[,1],rs.2$x[,1])
# 
# plot(rs.1$x[,2],rs.2$x[,2])
# 
# 
# 
# pheatmap(cor.matrix[names(sort(rs$x[,1])), names(sort(rs$x[,1]))],cluster_rows = F,cluster_cols = F)
# 
# 
# #pheatmap(CCLE.log2.rpkm.matrix[marker.gene,1:100],cluster_rows = T ,cluster_cols = T,clustering_distance_cols = 'correlation')
# 
# # require(RUVSeq)
# # 
# # rs  <- RUVg(x = 2^log2.read.count.matrix,k=1)
# # nd  <- rs$normalizedCounts
# # pca.rs <- prcomp(log2(nd+1) %>% t)
# # plot(pca.rs$x[,1:2],col=c('red','purple','blue')[as.integer(annotation.df$sample.type %>% factor)],pch=19)
# # 
# # 
# #pca.rs <- prcomp(log2.read.count.matrix %>% t)
# #pairs(pca.rs$x[,c(1:6)],col=c('red','purple')[( rownames(pca.rs$x) %in% weired.sample.name ) %>% factor %>% as.integer],pch=19)
# #plot(pca.rs$x[,1:2],col=c('red','purple','blue')[as.integer(annotation.df$sample.type %>% factor)],pch=19)
# # pca.rs <- prcomp(log2.read.count.matrix %>% t)
# # plot(pca.rs$x[,1:2],col=c('red','purple','blue')[as.integer(annotation.df$sample.type %>% factor)],pch=19)
# # 
# # 
# # pca.rs <- prcomp(log2.read.count.matrix[1:1000,] %>% t)
# # plot(pca.rs$x[,1:2],col=c('red','purple','blue')[as.integer(annotation.df$sample.type %>% factor)],pch=19)
# # 
# # pca.rs <- prcomp(log2.read.count.matrix[1001:2000,] %>% t)
# # plot(pca.rs$x[,1:2],col=c('red','purple','blue')[as.integer(annotation.df$sample.type %>% factor)],pch=19)
# # 
# # pca.rs <- prcomp(log2.read.count.matrix[2001:3000,] %>% t)
# # plot(pca.rs$x[,1:2],col=c('red','purple','blue')[as.integer(annotation.df$sample.type %>% factor)],pch=19)
# 
# f <- CCLE.sample.meta$cell.line.tumor.site == 'SKIN' 
# explore.batch.effects(CCLE.log2.read.count.matrix[,f],CCLE.sample.meta$cell.line.tumor.site[f])

