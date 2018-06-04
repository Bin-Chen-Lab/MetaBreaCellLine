load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
require(ggplot2)
require(pheatmap)
require(quantreg)
require(plyr)
require(dplyr)
require(PairedData)
require(foreach)



#############Some functions used in the analysis ############
pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
    marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))  
    marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene)) 
    correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
    cell.line.median.cor  <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
    best.cell.line        <- names(cell.line.median.cor)[1]
    p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
        v       <- correlation.matrix[,cell.line]
        p.value <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
      
    }
    names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
    fdr.vec            <- p.adjust(p.value.vec,method='fdr')
    list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}






#####   Pick out Breast Invasive Ductal Carcinoma samples ######
# There are two RNASeq protocals used in the MET500 project:polyA and hybrid#
# I will mainly focus on polyA, but the results of hybrid is similar #

flag                           <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') & MET500.sample.meta$LibrarySelection == 'PolyA'
polyA.breast.cancer.sample     <- rownames(MET500.sample.meta)[flag] %>% as.character
flag                           <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') & MET500.sample.meta$LibrarySelection == 'Hybrid Selection'
hybrid.breast.cancer.sample    <- rownames(MET500.sample.meta)[flag] %>% as.character


#####   Pick out cell line which shows largest correlation with samples ######
polyA.rs  <- pick.out.cell.line(MET500.log2.fpkm.matrix[,polyA.breast.cancer.sample], CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
hybrid.rs <- pick.out.cell.line(MET500.log2.fpkm.matrix[,hybrid.breast.cancer.sample],CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)

#### Show the results are highly correlated. Therefore, I only care PolyA protocal ####
plot(x=polyA.rs$cell.line.median.cor,
     y=hybrid.rs$cell.line.median.cor[names(polyA.rs$cell.line.median.cor)],
     xlab='PolyA',ylab='Hybrid',pch=19,cex=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(0,0.45),ylim=c(0,0.45)
     )
lines(c(0,0.45),c(0,0.45),lwd=3,lty=2)

###### Draw a plot to show that top-ranked cell lines are all breast cancer cell lines########
v   <- polyA.rs$cell.line.median.cor %>% sort
col <- ifelse(grepl(names(v),pattern='BREAST'),'red','black')
plot(v,col=col,cex=0.5,pch=19,cex.lab=1.5,cex.axis=1.5,xlab='rank',ylab='median.of.correlation')
points(x=(1:length(v))[col=='red'],v[col=='red'],col='red',cex=0.9,pch=19)
breast.cancer.cell.line.median.cor <-v[col == 'red'] %>% sort(decreasing = TRUE)
View(breast.cancer.cell.line.median.cor) # show the breast cancer cell line


###### Run paired wilcoxcon rank test to show MCF7 cell line is not the best ############
MCF7     <- polyA.rs$correlation.matrix[,'MCF7_BREAST']
MDMAB415 <- polyA.rs$correlation.matrix[,'MDAMB415_BREAST']
EFM192A  <- polyA.rs$correlation.matrix[,'EFM192A_BREAST']
BT483    <- polyA.rs$correlation.matrix[,'BT483_BREAST']
plot(paired(MCF7,MDMAB415),type='profile')+  theme(axis.text=element_text(size=30,face='bold'))
wilcox.test(x=MCF7,y=MDMAB415,paired = TRUE)
plot(paired(MCF7,EFM192A),type='profile')+  theme(axis.text=element_text(size=30,face='bold'))
wilcox.test(x=MCF7,y=EFM192A,paired = TRUE)
plot(paired(MCF7,BT483),type='profile')+  theme(axis.text=element_text(size=30,face='bold'))
wilcox.test(x=MCF7,y=BT483,paired = TRUE)


####  In the above analysis, all samples from different biopsysites were pooled together          ######
#### Let us do the analysis for each biopsy site, check whether we should use different cell line ###### 
### for different biopsy sites.                                                                    ######

biopsy.site                                  <-  MET500.sample.meta[polyA.breast.cancer.sample,'biopsy.site']
polyA.breast.cancer.liver.biopsy.sample      <-  polyA.breast.cancer.sample[biopsy.site == 'LIVER']
polyA.breast.cancer.lymph.node.biopsy.sample <-  polyA.breast.cancer.sample[biopsy.site == 'LYMPH_NODE']

polyA.liver.rs      <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,polyA.breast.cancer.liver.biopsy.sample],     expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)
polyA.lymph.node.rs <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,polyA.breast.cancer.lymph.node.biopsy.sample],expr.of.cell.lines = CCLE.log2.rpkm.matrix,CCLE.rna.seq.marker.gene.1000)

v   <- polyA.liver.rs$cell.line.median.cor %>% sort
col <- ifelse(grepl(names(v),pattern='BREAST'),'red','black')
plot(v,col=col,cex=0.5,pch=19,cex.lab=2,cex.axis=2,xlab='rank',ylab='median.of.correlation')
points(x=(1:length(v))[col=='red'],v[col=='red'],col='red',cex=0.9,pch=19)

v   <- polyA.lymph.node.rs$cell.line.median.cor %>% sort
col <- ifelse(grepl(names(v),pattern='BREAST'),'red','black')
plot(v,col=col,cex=0.5,pch=19,cex.lab=2,cex.axis=2,xlab='rank',ylab='median.of.correlation')
points(x=(1:length(v))[col=='red'],v[col=='red'],col='red',cex=0.9,pch=19)


col <- ifelse(grepl(x=names(polyA.liver.rs$cell.line.median.cor),pattern = 'BREAST'),'red','black')
plot(x=polyA.liver.rs$cell.line.median.cor,
     y=polyA.lymph.node.rs$cell.line.median.cor[names(polyA.liver.rs$cell.line.median.cor)],
     xlab='Liver',ylab='Lymph.node',cex=1,cex.lab=1.5,cex.axis=1.5,pch=19,col=col,xlim=c(0,0.45),ylim=c(0,0.45))
lines(c(0,0.45),c(0,0.45),lty=2,lwd=2)

intersect(polyA.lymph.node.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>% head(10),
          polyA.liver.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>% head(10)
          )
intersect(polyA.lymph.node.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>% head(20),
          polyA.liver.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>% head(20)
) # well, the results are quite similar


######### According to the above analysis, MDAMB415_BREAST seems to be the best#######
### Take it as example, check whether samples obtained from different biopsy sites correlate differently with it ####
liver.df      <- data.frame(biopsy.site='LIVER',     cor.value=polyA.liver.rs$correlation.matrix[,'MDAMB415_BREAST'])
lymph.node.df <- data.frame(biopsy.site='LYMPH_NODE',cor.value=polyA.lymph.node.rs$correlation.matrix[,'MDAMB415_BREAST'])

#boxplot of correlation value with MDAMB415_BREAST and rank test, no significant difference
draw.df       <- rbind(liver.df,lymph.node.df)
ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE'),]) + geom_boxplot(aes(y=cor.value,x=biopsy.site,color=biopsy.site),size=3) + theme(axis.text=element_text(size=21),axis.title=element_text(size=21),legend.text=element_text(size=21))
p.value <- wilcox.test(liver.df$cor.value,lymph.node.df$cor.value,paired=FALSE)$p.value



######## load tumor purity data ########
MET500.tumor.content <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/MET500.tumor.content.csv", stringsAsFactors=FALSE)
flag                 <- MET500.tumor.content$sequecing.tumor.content == "<30%"
MET500.tumor.content[flag,'sequecing.tumor.content'] <- '-1%'
MET500.tumor.content$sequecing.tumor.content <- sub(x = MET500.tumor.content$sequecing.tumor.content,replacement = '',pattern = '%') %>% as.numeric
rownames(MET500.tumor.content)               <- MET500.tumor.content$MET.500.id

polyA.breast.cancer.sample.MET500.id      <- MET500.sample.meta[polyA.breast.cancer.sample %>% as.character,'MET500.id']
polyA.breast.cancer.sample.tumor.content  <- MET500.tumor.content[polyA.breast.cancer.sample.MET500.id,'sequecing.tumor.content'] / 100
cor.value                                 <- polyA.rs$correlation.matrix[polyA.breast.cancer.sample,'MDAMB415_BREAST']

flag <- polyA.breast.cancer.sample.tumor.content > 0
plot(polyA.breast.cancer.sample.tumor.content[flag],cor.value[flag],xlim=c(0,1),ylim=c(0,0.6),xlab='tumor.purity',ylab='correlation.with.MDAMB415',cex.lab=1.5,cex.axis=1.5,pch=19)


# fit.data       <- data.frame(tumor.content=polyA.breast.cancer.sample.tumor.content,cor.value = cor.value)
# fit.data       <- fit.data[fit.data$tumor.content >0,]
# plot(fit.data$tumor.content,fit.data$cor.value,xlim=c(0,1),ylim=c(0,0.6))
# 
# 
# line.fit <- lm(formula=cor.value ~ tumor.content,data=fit.data)
# plot(x=fit.data$tumor.content,y=fit.data$cor.value,xlim=c(0,1),ylim=c(0,0.6),cex.lab=2.5,xlab='tumor.content',ylab='cor.value',cex.axis=1.5,pch=19)
# lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)
# summary(line.fit)# p-value of regression
# 
# line.fit <- rq(formula=cor.value ~ tumor.content,data=fit.data,tau=0.5)
# plot(x=fit.data$tumor.content,y=fit.data$cor.value,xlim=c(0,1),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='cor.value',cex.axis=1.5,pch=19)
# lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)
# summary(line.fit, se = "boot")# p-value of regression


polyA.breast.cancer.sample.rank.vec <- foreach(s = polyA.breast.cancer.sample,.combine='c') %do% {
    v <- polyA.rs$correlation.matrix[s,]      
    rank(v)['MDAMB415_BREAST']
}
names(polyA.breast.cancer.sample.rank.vec) <- polyA.breast.cancer.sample
#rank.vec <- sort(rank.vec)

plot(x=polyA.breast.cancer.sample.tumor.content[flag],y=polyA.breast.cancer.sample.rank.vec[flag],cex.lab=1.5,cex.axis=1.5,ylab='rank.of.MDAMB415',xlab='tumor.content',pch=19)
lowess(x=polyA.breast.cancer.sample.tumor.content[flag],y=polyA.breast.cancer.sample.rank.vec[flag]) %>% lines(lwd=5,col='red')


z.score <-( polyA.breast.cancer.sample.rank.vec - median(polyA.breast.cancer.sample.rank.vec) )/ mad(polyA.breast.cancer.sample.rank.vec)

specious.sample <- names(polyA.breast.cancer.sample.rank.vec)[z.score < -3]
abline(h=polyA.breast.cancer.sample.rank.vec[specious.sample] %>% max)
# col <- ifelse(rownames(fit.data) %in% specious.sample,'red','black')
# plot(fit.data$tumor.content,fit.data$cor.value,col=col,pch=19)
# 
# ## well , let us arbitraily remove some points , the regression line is much better!##
# ## there might be heterogenerity in the data, we need to further confirm###
# fit.data <- fit.data[ !(fit.data$tumor.content > 80 & fit.data$cor.value<0.4),]
# line.fit <- lm(formula=cor.value ~ tumor.content,data=fit.data)
# plot(x=fit.data$tumor.content,y=fit.data$cor.value,xlim=c(0,100),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='cor.value',cex.axis=1.5,pch=19)
# lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5,pch=19)
# summary(line.fit)





################# Analysis with microarray data ######################
#                                                                 #
#                                                                 #
####################################################################
load('server-side/RData/microarray.RData')

tmp                         <- CCLE.microarray.expr.matrix
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.microarray.marker.gene <- names(sort(rank.sd,decreasing =TRUE))[1:1000]

microarray.rs <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix,expr.of.cell.lines = CCLE.microarray.expr.matrix,marker.gene = CCLE.microarray.marker.gene)


col.name <- colnames(microarray.rs$correlation.matrix)
rank.vec <- foreach(i=1:nrow(microarray.rs$correlation.matrix),.combine='c') %do% {
    v <-  microarray.rs$correlation.matrix[i,]
    rank(v)['ZR7530_BREAST']
    
}
names(rank.vec) <- rownames(microarray.rs$correlation.matrix)

v   <- microarray.rs$cell.line.median.cor %>% sort
col <- ifelse(grepl(names(v),pattern='BREAST'),'red','black')
plot(v,col=col,cex=0.5,pch=19,cex.lab=2,cex.axis=2,xlab='rank',ylab='median.of.correlation')
points(x=(1:length(v))[col=='red'],v[col=='red'],col='red',cex=0.9,pch=19)



#### differnt biopsy site, different cell line?#######
biopsy.site       <-  metastatic.breast.cancer.microarray.meta.df$biopsy.site
liver.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LIVER']  %>% as.character
lung.samples      <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LUNG']   %>% as.character
bone.samples      <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'BONE']   %>% as.character
brain.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'BRAIN']  %>% as.character


microarray.liver.rs  <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,liver.samples], expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene)
microarray.lung.rs   <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,lung.samples],  expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene)
microarray.bone.rs   <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,bone.samples],  expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene)
microarray.brain.rs  <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,brain.samples], expr.of.cell.lines = CCLE.microarray.expr.matrix,CCLE.microarray.marker.gene)

liver.cell.line <- microarray.liver.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>% head(10)
lung.cell.line <- microarray.lung.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>%head(10)
bone.cell.line <- microarray.bone.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>%head(10)
brain.cell.line <- microarray.brain.rs$cell.line.median.cor %>% sort(decreasing = TRUE) %>% names %>% head(10)
require(gplots)
venn(list(liver=liver.cell.line,lung=lung.cell.line,bone=bone.cell.line,brain.cell.line))


require(data.table)
sum.df.score   <- fread('./client-side/output/estimate/sum.df.score.gct',skip=2) %>% as.data.frame
estimate.score <- sum.df.score[4,names(rank.vec)] %>% unlist

col.vec <- ifelse(names(rank.vec) %in% bone.samples,'red','black')
col.vec <- ifelse(names(rank.vec) %in% brain.samples,'red','black')
col.vec <- ifelse(names(rank.vec) %in% liver.samples,'red','black')
col.vec <- ifelse(names(rank.vec) %in% lung.samples,'red','black')


plot(x=estimate.score,y=rank.vec,xlab='estimate.score',ylab='rank of ZR7530',cex=1.5,pch=19,cex.axis=1.5,cex.lab=1.5)
#lowess(x=estimate.score,y=rank.vec) %>% lines(lwd=5,col='red')
abline(h=median(rank.vec) - 3 * mad(rank.vec))

z.score          <- (rank.vec - median(rank.vec)) / mad(rank.vec)
specious.samples <- names(z.score)[z.score < -3]

specious.samples.rank.vec <- rank.vec[specious.samples] %>% sort



############### Load TCGA data, show that many samples match best to HS822T_BONE
load('./server-side/RData/TCGA.breast.cancer.RData')
TCGA.rs <- pick.out.cell.line(expr.of.samples = TCGA.breast.cancer.log2.fpkm.matrix,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = marker.gene.1000)
cor.rank <- apply(TCGA.rs$correlation.matrix,1,rank) %>% t
View(sort(cor.rank[,'HS822T_BONE'],decreasing = TRUE))

TCGA.breast.cancer.log2.fpkm.matrix.rank <- apply(TCGA.breast.cancer.log2.fpkm.matrix,2,rank)
var.value <- apply(TCGA.breast.cancer.log2.fpkm.matrix.rank,1,mad) %>% sort(decreasing = TRUE)
top.gene <- names(var.value)[1:100]
pca.rs <- prcomp(TCGA.breast.cancer.log2.fpkm.matrix.rank[top.gene,] %>% t)
plot(pca.rs$x[,1:2])


##########################################  Trash  ###############################################
######## PCA plot of CCLE cell line, the PC1 obtained by all genes or PAM50 genes are consistent #####
f <- grepl(x=colnames(CCLE.log2.rpkm.matrix),pattern='BREAST')
breast.cancer.cell.line <- colnames(CCLE.log2.rpkm.matrix)[f]
pca.rs <- prcomp(CCLE.log2.rpkm.matrix[,breast.cancer.cell.line] %>% t)
plot(pca.rs$x[,1:2])
breast.cancer.cell.line.df <- data.frame(cell.line.name=breast.cancer.cell.line,
                                         PC1 = pca.rs$x[breast.cancer.cell.line,1],
                                         PC2 = pca.rs$x[breast.cancer.cell.line,2],
                                         median.cor = polyA.cell.line.median.cor[breast.cancer.cell.line]
                                         )
ggplot(breast.cancer.cell.line.df) + geom_text(aes(x=PC1,y=PC2,color=median.cor,label=cell.line.name),size=4) + scale_color_gradient(low='black',high='yellow')





pam50.gene <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')$ensemble.gene.id %>% as.character
pca.rs <- prcomp(CCLE.log2.rpkm.matrix[pam50.gene,breast.cancer.cell.line] %>% t)
df <- data.frame(
           pam50.PC1 = pca.rs$x[breast.cancer.cell.line,1],
           pam50.PC2 = pca.rs$x[breast.cancer.cell.line,2]
   )
breast.cancer.cell.line.df <- cbind(breast.cancer.cell.line.df,df)
ggplot(breast.cancer.cell.line.df) + geom_point(aes(x=pam50.PC1,y=pam50.PC2,color=median.cor),size=4) + scale_color_gradient(low='black',high='yellow')
r1 <- pca.rs$rotation[,'PC1'] %>% sort
CCLE.log2.rpkm.matrix['ENSG00000129514',breast.cancer.cell.line] %>% sort %>% plot


ggplot(breast.cancer.cell.line.df) + geom_point(aes(y=pam50.PC1,x=PC1,color=median.cor),size=4) + scale_color_gradient(low='black',high='yellow')

####### Do similar analysis on TCGA breast cancer data########
load('server-side/RData/breast.cancer.RData')

common.gene              <- intersect(marker.gene,rownames(TCGA.breast.cancer.log2.fpkm.matrix))
cor.matrix               <- cor(TCGA.breast.cancer.log2.fpkm.matrix[common.gene,],CCLE.log2.rpkm.matrix[common.gene,],method='spearman')
cell.line.cor.median.with.TCGA.breast.cancer.sample <- apply(cor.matrix,2,median) %>% sort(decreasing = TRUE)

########

# TCGA.tumor.purity <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/TCGA.tumor.purity.csv", stringsAsFactors=FALSE)
# TCGA.tumor.purity$Sample.ID <- gsub(pattern = 'A$',replacement = '',x=TCGA.tumor.purity$Sample.ID)
# rownames(TCGA.tumor.purity) <- TCGA.tumor.purity$Sample.ID
# 
# purity <- TCGA.tumor.purity[colnames(TCGA.breast.cancer.log2.fpkm.matrix),'ESTIMATE'] * 100
# cor.vec <- cor.matrix[colnames(TCGA.breast.cancer.log2.fpkm.matrix),'MDAMB415_BREAST']
# 
# fit.data       <- data.frame(tumor.content=purity ,cor.value = cor.vec)
# fit.data       <- fit.data[complete.cases(fit.data),]
# 
# line.fit <- lm(formula=cor.value ~ tumor.content,data=fit.data)
# plot(x=fit.data$tumor.content ,y=fit.data$cor.value,xlim=c(0,100),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='cor.value',cex.axis=1.5,pch=19)
# lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)
# summary(line.fit)# p-value of regression
# 
# 
# line.fit <- rq(formula=cor.value ~ tumor.content,data=fit.data,tau=0.5)
# plot(x=fit.data$tumor.content ,y=fit.data$cor.value,xlim=c(0,100),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='cor.value',cex.axis=1.5,pch=19)
# lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)
# summary(line.fit)# p-value of regression
# 
# 
# 
# rank.vec <- foreach(s = colnames(TCGA.breast.cancer.log2.fpkm.matrix),.combine='c') %do% {
#   v <- cor.matrix[s,]      
#   rank(v)['MDAMB415_BREAST']
# }
# names(rank.vec) <- colnames(TCGA.breast.cancer.log2.fpkm.matrix)
# 
# z <- (rank.vec - median(rank.vec) )/ mad(rank.vec)
# 
# 
# hehe.matrix<- foreach(s = rownames(cor.matrix),.combine='rbind') %do% {
#   v <- cor.matrix[s,]   %>% sort    
#   ifelse(grepl(x=names(v),pattern = 'BREAST'),1,0)
# }
# rownames(hehe.matrix) <- rownames(cor.matrix)
# 
# 
# 
# 
# # good.breast.cancer.cell.line <- names(cell.line.median.cor)[cell.line.median.cor > median(cell.line.median.cor) + 1.5 * IQR(cell.line.median.cor)]    
# # good.breast.cancer.cell.line <- good.breast.cancer.cell.line[grepl(x=good.breast.cancer.cell.line,pattern = 'BREAST')]
# # 
# # good.breast.cancer.cell.line <- names(cell.line.median.cor)[cell.line.median.cor > 0]    
# # good.breast.cancer.cell.line <- good.breast.cancer.cell.line[grepl(x=good.breast.cancer.cell.line,pattern = 'BREAST')]
# 
# 
# good.breast.cancer.cell.line <- tail(names(sort(cell.line.median.cor)),1)
# 
# data                  <- cor.matrix[polyA.breast.cancer.sample,good.breast.cancer.cell.line]
# pheatmap(data,show_rownames = F,show_colnames = T,cluster_cols = T)
# 
# #data          <- breast.cor.matrix[,breast.cancer.cell.line]
# #pheatmap(data,show_rownames = F,show_colnames = T,cluster_cols = T,clustering_distance_rows = 'correlation')
# 
# #pca.rs <- prcomp(data)
# #plot(pca.rs$x[,1:2])
# 
# 
# biopsy.site <- MET500.sample.meta[ rownames(data),'biopsy.site']
# #pheatmap(data[biopsy.site == 'LIVER',],show_rownames = F,show_colnames = T,cluster_cols = T)
# #pheatmap(data[biopsy.site == 'BONE',],show_rownames = F,show_colnames = T,cluster_cols = T)
# #pheatmap(data[biopsy.site == 'LYMPH_NODE',],show_rownames = F,show_colnames = T,cluster_cols = T)
# 
# 
# # pca.rs         <- prcomp(data)
# # plot(pca.rs$x[,1:2])
# # pc1 <- pca.rs$x[,1] 
# # pc2 <- pca.rs$x[,2] 
# draw.df <- data.frame(biopsy.site = biopsy.site %>% as.character,mean.cor=apply(data,1,mean))
# 
# breast.cancer.sample.MET500.id   <- MET500.sample.meta[rownames(draw.df) %>% as.character,'MET500.id']
# tumor.content                    <- MET500.tumor.content[breast.cancer.sample.MET500.id,'sequecing.tumor.content']
# draw.df$tumor.content            <- tumor.content
# 
# #ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE','BONE'),]) + geom_point(aes(x=PC1,y=PC2,color=biopsy.site),size=3)
# ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE'),]) + geom_boxplot(aes(y=mean.cor,x=biopsy.site,color=biopsy.site),size=3) + theme(axis.text=element_text(size=21),axis.title=element_text(size=21),legend.text=element_text(size=21))
# ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE'),]) + geom_boxplot(aes(y=tumor.content,x=biopsy.site,color=biopsy.site),size=3)+ theme(axis.text=element_text(size=21),axis.title=element_text(size=21),legend.text=element_text(size=21))
# 
# 
# t.test(     draw.df$mean.cor[draw.df$biopsy.site == 'LIVER'],    draw.df$mean.cor[draw.df$biopsy.site == 'LYMPH_NODE'])
# wilcox.test(draw.df$mean.cor[draw.df$biopsy.site == 'LIVER'],    draw.df$mean.cor[draw.df$biopsy.site == 'LYMPH_NODE'])
# 
# t.test(     draw.df$tumor.content[draw.df$biopsy.site == 'LIVER'],draw.df$tumor.content[draw.df$biopsy.site == 'LYMPH_NODE'])
# wilcox.test(draw.df$tumor.content[draw.df$biopsy.site == 'LIVER'],draw.df$tumor.content[draw.df$biopsy.site == 'LYMPH_NODE'])
# 
# draw.df  <- draw.df[draw.df$tumor.content >0,]
# fit.data <- draw.df
# line.fit <- rq(formula=mean.cor ~ tumor.content,data=fit.data,tau=0.5)
# plot(x=fit.data$tumor.content,y=fit.data$mean.cor,xlim=c(0,100),ylim=c(0,0.6),cex.lab=2,xlab='tumor.content',ylab='mean.cor',cex.axis=1.5,pch=19)
# lines(x=fit.data$tumor.content,y=line.fit$fitted.values,col='red',lwd=5)
# summary(line.fit, se = "boot")# p-value of regression
# draw.df$residual <- line.fit$residuals
# ggplot(draw.df[draw.df$biopsy.site %in% c('LIVER','LYMPH_NODE'),]) + geom_boxplot(aes(y=residual,x=biopsy.site,color=biopsy.site),size=3) + theme(axis.text=element_text(size=21),axis.title=element_text(size=21),legend.text=element_text(size=21))
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
# draw.df <- draw.df[breast.cancer.sample,]
# draw.df$cancer.type <- MET500.sample.meta[match( rownames(draw.df) %>% as.character, MET500.sample.meta$Run %>% as.character),'cancer.type']
# 
# #### check expression of breast marker genes
# GTEX.tissue.specific.gene <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/GTEX.tissue.specific.gene.csv", stringsAsFactors=FALSE)
# breast.specific.gene    <- GTEX.tissue.specific.gene$gene_id[GTEX.tissue.specific.gene$Tissue == 'Breast' ] %>% as.character
# remove.ensemble.version.id <- function(x){
#   v <- strsplit(x = x,split = "\\.") %>% unlist  
#   v[1]
# }
# breast.specific.gene <- sapply(breast.specific.gene,remove.ensemble.version.id) %>% as.character
# 
# 
# 
# 
# data   <- CCLE.log2.rpkm.matrix[breast.specific.gene,breast.cancer.cell.line]
# data   <- CCLE.log2.rpkm.matrix[,breast.cancer.cell.line]
# 
# pca.rs <- prcomp(data %>% t)
# plot(pca.rs$x[,1:2],cex.axis=1.5,pch=19,cex.lab=1.5)
# pc1 <- pca.rs$x[,1]
# pc2 <- pca.rs$x[,2]
# r1 <- pca.rs$rotation[,'PC1'] %>% sort
# r2 <- pca.rs$rotation[,'PC2'] %>% sort
# 
# CCLE.log2.rpkm.matrix['ENSG00000240800',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000235687',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000160182',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000235123',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000132746',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000167754',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000138685',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000164692',breast.cancer.cell.line] %>% sort %>% plot
# CCLE.log2.rpkm.matrix['ENSG00000106366',breast.cancer.cell.line] %>% sort %>% plot
# 
# 
# 
# hehe <- breast.specific.gene[breast.specific.gene %in% rownames(MET500.log2.fpkm.matrix)]
# data <- MET500.log2.fpkm.matrix[hehe,polyA.breast.cancer.sample] 
# pca.rs <- prcomp(data %>% t)
# plot(pca.rs$x[,1:2])






# hybrid.lib.id <- rownames(MET500.sample.meta)[MET500.sample.meta$LibrarySelection == 'Hybrid Selection']
# hybrid.lib.id <- hybrid.lib.id[hybrid.lib.id %in% breast.cancer.sample]
# h <- MET500.sample.meta[hybrid.lib.id %>% as.character,'MET500.id']
# 
# data <- MET500.log2.fpkm.matrix[,hybrid.lib.id] 
# tt <- MET500.tumor.content[h %>% as.character,'sequecing.tumor.content']
# data <- data[,tt>=70]
# 
# pca.rs <- prcomp(data %>% t)
# plot(pca.rs$x[,1:2])
# pc1 <- pca.rs$x[,1]

# CCLE.cell.line.idx     <- apply(breast.cor.matrix,1,function(x) which(x == max(x))[1])
# CCLE.max.cor           <- apply(breast.cor.matrix,1,function(x) max(x))
# rs.df                  <- cbind(sample.name = rownames(breast.cor.matrix),CCLE.max.cor=CCLE.max.cor,CCLE.sample.meta[CCLE.cell.line.idx,])
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

