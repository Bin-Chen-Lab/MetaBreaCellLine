require(dplyr)
require(Rtsne)
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output//CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')
load('client-side/output/correlation.analysis.with.expression.microarray.R.output/correlation.analysis.with.expression.microarray.RData')
source('client-side/code/for.figure/ggplot.style.R')
require('PerformanceAnalytics')
source('client-side/code/for.figure/chart.Correlation.R') # I steal the code from PerformanceAnalytics package and modify it,remove the weired stars


######################################################################################################################
#   
#  Fig 3
#  
######################################################################################################################


######################################
### Fig 3a: Show ranking of ALL CCLE cell lines according totranscriptome correlation analysis
######################################
draw.df <- data.frame(y = MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor,
                      x = length(MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor):1,
                      is.breast.cancer.cell.line = ifelse( names(MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor) %in% CCLE.breast.cancer.cell.line,'Y','N'),
                      size =ifelse(names(MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor) %in% CCLE.breast.cancer.cell.line,4,0.5)
                      )
ggplot(draw.df,aes(x=x,y=y,size=size,color=factor(is.breast.cancer.cell.line))) + geom_point(show.legend=F) + 
xlab('rank') + ylab('transcriptome similarity') + 
scale_color_manual(values=c('Y' = 'red','N' = 'grey')) + ggplot.style + scale_size_area(max_size = 4)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/MET500.vs.CCLE.polyA.expression.correlation.result.pdf',width=20,height=15)




######################################
### Fig 3b Show biopsy-site-specific transcriptome correaltion analysis resutls are correlated between liver and lymph node
######################################
draw.df <- data.frame(y = MET500.vs.CCLE.polyA.liver.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                      x = MET500.vs.CCLE.polyA.lymph.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(draw.df,aes(x=x,y=y)) + geom_point(show.legend=F,size=6) + 
xlab('Lymph Node') + ylab('Liver') + geom_abline(slope = 1,intercept = 0) + 
xlim(0,0.5) +ylim(0,0.5) +ggplot.style 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/MET500.result.correlation.between.biopsy.site.pdf',width=20,height=10)
sp.correlation <- cor(draw.df$y,draw.df$x,method='spearman')




########################################
### Fig 3c: tsne plot to show that PAM50 subtype information is still maintained, DOES NOT depend on tissues (MET500)
########################################
load('server-side/RData/MET500.RData')
require(Rtsne)
pam50.gene.df  <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
pam50.gene     <- pam50.gene.df$ensemble.gene.id %>% as.character
dist.obj       <- as.dist(1- cor(MET500.log2.fpkm.matrix[pam50.gene,MET500.breast.cancer.polyA.sample],method='spearman'))
set.seed(8) # I want to reproduce the tsne results, 12 is just a arbitrary numnber, it DOES NOT change the conclusion
tsne.rs        <- Rtsne(dist.obj,perplexity = 15)
MET500.pam50.subtype.vec <- c( rep(x='LuminalA',  times= MET500.breast.cancer.polyA.LumA.sample %>% length),
                               rep(x='LuminalB',  times= MET500.breast.cancer.polyA.LumB.sample %>% length),
                               rep(x='Her2-enriched',  times= MET500.breast.cancer.polyA.Her2.sample  %>% length),
                               rep(x='Basal-like', times= MET500.breast.cancer.polyA.Basal.sample %>% length),
                               rep(x='Normal-like',times= MET500.breast.cancer.polyA.Normal.sample %>% length)
)
names(MET500.pam50.subtype.vec) <- c(MET500.breast.cancer.polyA.LumA.sample,
                                     MET500.breast.cancer.polyA.LumB.sample,
                                     MET500.breast.cancer.polyA.Her2.sample,
                                     MET500.breast.cancer.polyA.Basal.sample,
                                     MET500.breast.cancer.polyA.Normal.sample
)
draw.df <- data.frame(dim1=tsne.rs$Y[,1],
                      dim2=tsne.rs$Y[,2],
                      subtype = MET500.pam50.subtype.vec[MET500.breast.cancer.polyA.sample],
                      biopsy.site   = MET500.sample.meta[MET500.breast.cancer.polyA.sample,'biopsy.site']
)
ggplot(draw.df,aes(x=dim1,y=dim2,shape=subtype,color=biopsy.site)) + geom_point(size=6) +  scale_shape_manual(values=c(8,15:18))+
  ggplot.style+ theme(legend.text=element_text(size=20,face='bold'),legend.title=element_text(size=20,face = 'bold'),legend.key.size = unit(1.5, 'lines'))+
  ylim(-16,16)  
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/tsne.pam50.subtype.pdf',width=20,height=10)




###############################################
######Fig 3d.Show the correlation of expression correlation analysis with samples of different subtype
###############################################


mat <- cbind( MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line], 
              MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
colnames(mat) <- c('Basal-like','LuminalB','LuminalA','Her2\nenriched')

pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/MET500.result.correlation.between.subtype.pdf',width=20,height=20)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.cor.tuning = 0.7,cor.range = c(-1,1))
dev.off()
my_palette <- colorRampPalette(c("blue", "red"))(n = 101)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)
  
    dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }
}
pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/color.bar.pdf',width=20,height=20)
color.bar(my_palette,min = 0.2,max=1) # draw the color bar
dev.off()


######################################################################################################################
#   
#  Fig S2
#  
######################################################################################################################



######################################
### Fig S2a  Boxplot showing that samples from liver and lymph node DO NOT show significant differentce co-expression with  cell line MDAMB415
######################################

biopsy.site      <-c( rep(times=MET500.vs.CCLE.polyA.liver.expression.correlation.result$correlation.matrix %>% nrow,x='Liver'),
                      rep(times=MET500.vs.CCLE.polyA.lymph.expression.correlation.result$correlation.matrix %>% nrow,x='Lymph node')
                    )
MDAMB415.cor.vec <- c(MET500.vs.CCLE.polyA.liver.expression.correlation.result$correlation.matrix[,'MDAMB415_BREAST'],
                      MET500.vs.CCLE.polyA.lymph.expression.correlation.result$correlation.matrix[,'MDAMB415_BREAST']
                      )

draw.df <- data.frame(MDAMB415.cor.vec = MDAMB415.cor.vec,
                      biopsy.site = biopsy.site
                      )
ggplot(draw.df,aes(x=factor(biopsy.site),y=MDAMB415.cor.vec)) + geom_boxplot(outlier.shape = NA,lwd=1.5) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=biopsy.site))+ ggplot.style + xlab('') +ylab('')
wilcox.test(MET500.vs.CCLE.polyA.liver.expression.correlation.result$correlation.matrix[,'MDAMB415_BREAST'],MET500.vs.CCLE.polyA.lymph.expression.correlation.result$correlation.matrix[,'MDAMB415_BREAST'])
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/boxplot.liver.vs.lymph.with.MDAMB415.cor.pdf',width=20,height=10)




########################################
### Fig S2b: tsne plot to show that PAM50 subtype information is still maintained, DOES NOT depend on tissues (TCGA and MET500)
########################################
load('server-side/RData/TCGA.breast.cancer.RData')
load('server-side/RData/MET500.RData')
require(genefu)
combined.data               <- cbind(MET500.log2.fpkm.matrix[pam50.gene,MET500.breast.cancer.polyA.sample],TCGA.breast.cancer.log2.fpkm.matrix[pam50.gene,])
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- combined.data %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )

dist.obj       <- as.dist(1- cor(combined.data,method='spearman'))
set.seed(8) # I want to reproduce the tsne results, 12 is just a arbitrary numnber, it DOES NOT change the conclusion
tsne.rs        <- Rtsne(dist.obj)
plot(tsne.rs$Y)
mapping.vec <- c('LumA'='LuminalA','LumB' = 'LuminalB','Normal'= 'Normal-like','Basal'='Basal-like','Her2'='Her2-enriched')
draw.df     <- data.frame(dim1=tsne.rs$Y[,1],dim2=tsne.rs$Y[,2],subtype=pam50.subtype.rs$subtype)
draw.df$data.source <- ifelse(grepl(x=rownames(draw.df),pattern='SRR'),'MET500','TCGA')
draw.df$subtype     <- mapping.vec[draw.df$subtype %>% as.character]

ggplot(draw.df,aes(x=dim1,y=dim2,shape=subtype,color=data.source)) + geom_point(size=4) +  scale_shape_manual(values=c(8,15:18))+scale_color_manual(values=c("MET500" =  "#800000", "TCGA" = "#008080"))+
ggplot.style+ theme(legend.text=element_text(size=20,face='bold'),legend.title=element_text(size=20,face = 'bold'),legend.key.size = unit(1.5, 'lines'))  
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/tsne.pam50.subtype.MET500.and.TCGA.pdf',width=20,height=10)


######################################################################################################################
#   
#  Fig S3
#  
######################################################################################################################


##############################################
#### Fig S3a: consistency of biopsy-site specific correlation analysis between microarray and rna-seq
##############################################

#show the consistency of liver-specific transcripteom  correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.liver.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.liver.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/Liver.consistency.pdf',width=20,height=10)

#show the consistency of lymph-specific transcripteom  correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.lymph.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.lymph.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/Lymph.consistency.pdf',width=20,height=10)




##############################################
#### Fig S3b: consistency of sub-type specific correlation analysis between microarray and rna-seq
##############################################
#show the consistency of LumB expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/LumB.consistency.pdf',width=20,height=10)


#show the consistency of LumA expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/LumA.consistency.pdf',width=20,height=10)

#show the consistency of Her2 expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/Her2.consistency.pdf',width=20,height=10)

#show the consistency of Basal expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/Basal.consistency.pdf',width=20,height=10)



######################################################################################################################
#   
#  Fig S4
#  
######################################################################################################################

################################################################
### Fig S4, show gene expression correlation analysis results among differnt biopsy sites (microarray data)
################################################################

mat <- cbind( microarray.vs.CCLE.bone.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line], 
              microarray.vs.CCLE.liver.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              microarray.vs.CCLE.lymph.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              microarray.vs.CCLE.brain.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              microarray.vs.CCLE.lung.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
              
)
colnames(mat) <- c('Bone','Liver','Lymph\nNode','Brain','Lung')
pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/microarray.result.correlation.between.biopsy.site.pdf',width=20,height=20)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.cor.tuning = 0.7,cor.range = c(-1,1))
dev.off()




######################################################################################################################
#   
#  Fig S5
#  
######################################################################################################################



################################################################
### Fig S5, show gene expression correlation analysis results among differnt subtypes (microarray data)
###############################################################

mat <- cbind( microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line], 
              microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
              microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
              )
colnames(mat) <- c('Basal-like','LuminalB','LuminalA','Her2\nenriched')
pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/microarray.result.correlation.between.subtype.pdf',width=20,height=20)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.cor.tuning = 0.7,cor.range = c(-1,1))
dev.off()


######################################################################################################################
#   
#  Fig S6  just paste the microarray.data.tumor.purity.pdf into FigS6.ai
#  
######################################################################################################################



################################################
#### Show tumor purity of microarray data ##############
################################################
load('client-side/output/run.estimate.R.output/microarray.data.tumor.purity.RData')
microarray.data.tumor.purity.df$biopsy.site <- factor(microarray.data.tumor.purity.df$biopsy.site,levels = c('BRAIN','LYMPH_NODE','LIVER','LUNG','BONE'))
ggplot(microarray.data.tumor.purity.df,aes(x=biopsy.site,y=estimate.score),show.legned=FALSE) + geom_boxplot(lwd=1.2,outlier.shape = NA) + geom_point(size=4,aes(fill=biopsy.site),position=position_jitterdodge(jitter.width = 0.2)) + ggplot.style + xlab('') + theme(axis.text.x = element_text(face="bold",size=30)) + ylab('tumor purity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section2/microarray.data.tumor.purity.pdf',width=20,height=10)

wilcox.test(microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'BONE'],microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'BRAIN'])
wilcox.test(microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'BONE'],microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'LUNG'])
wilcox.test(microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'BONE'],microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'LIVER'])
wilcox.test(microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'BONE'],microarray.data.tumor.purity.df$estimate.score[microarray.data.tumor.purity.df$biopsy.site == 'LYMPH_NODE'])

