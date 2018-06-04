load('client-side/output/pick.out.cell.line.by.subtype.R.output//pick.out.cell.line.by.subtype.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')
load('client-side/output/correlation.analysis.with.expression.microarray.R.output/correlation.analysis.with.expression.microarray.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData')
load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.cnv.RData')
load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')

source('client-side/code/for.figure/ggplot.style.R')
require(dplyr)
require(reshape2)


MET500.sra.run                      <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/MET500.SraRunTable.txt", stringsAsFactors=FALSE)
MET500.sra.run                      <- MET500.sra.run[,c('Assay_Type','LibrarySelection','LibrarySource','Run','is_tumor','body_site','submitted_subject_id')]
MET500.sra.run$dbGap.subject.id     <- MET500.sra.run$submitted_subject_id
MET500.sra.run$submitted_subject_id <- NULL

MET500.to.dbGap                     <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/MET500_dbGaP_SubjectID.csv")

MET500.sample.phenotype                    <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/MET500.sample.phenotype.csv", stringsAsFactors=FALSE)
MET500.sample.phenotype$primary.tumor.site <- MET500.sample.phenotype$primary.site
MET500.sample.phenotype$primary.tumor.site <- toupper(MET500.sample.phenotype$primary.tumor.site)
MET500.sample.phenotype$primary.site       <- NULL
MET500.sample.phenotype                    <- MET500.sample.phenotype[,c('MET500.id','cancer.type','primary.tumor.site')]

MET500.meta             <- merge(x=MET500.sra.run,y=MET500.to.dbGap,by = 'dbGap.subject.id')
MET500.meta$biopsy.site <- toupper(MET500.meta$body_site)
MET500.meta$body_site   <- NULL
MET500.meta             <- merge(MET500.meta,MET500.sample.phenotype,by='MET500.id')
MET500.RNASeq.meta      <- MET500.meta[MET500.meta$Assay_Type =='RNA-Seq',]
id                      <- MET500.RNASeq.meta$MET500.id[match(MET500.breast.cancer.polyA.Basal.sample,MET500.RNASeq.meta$Run)] %>% as.character %>% unique

############################################################ 
#### Fig 3b: right  panel, boxplot of CNV correlation with Basal-like MET500 breast cancer samples, MDAMB231 and HCC70
###########################################################

CNV.cor.matrix <- cor(CCLE.gene.cnv[,CCLE.breast.cancer.cell.line],MET500.breast.cancer.gene.cnv[,id],method='spearman')
m              <- apply(CNV.cor.matrix,1,median) %>% sort
names(m)       <- gsub(x=names(m),replacement = '',pattern = '_BREAST')

df                     <- as.data.frame(CNV.cor.matrix)
df$cell.line.name      <- gsub(x=rownames(df),replacement = '',pattern = '_BREAST')
long.df                <- melt(df,id.vars = 'cell.line.name')
long.df$cell.line.name <- factor(long.df$cell.line.name,levels = names(m))
long.df$color          <- 'black'
long.df[long.df$cell.line.name == 'MDAMB231','color'] <- 'red'
long.df <- long.df[long.df$cell.line.name %in% c('MDAMB231','HCC70'),]

ggplot(long.df) + geom_boxplot(aes(x=cell.line.name,y=value,color=factor(color)),size=1.5,show.legend=F) + ggplot.style + ylab('correlation')+ xlab('') + scale_color_manual(values=c('red'='red','black'='black'))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/cnv.correlation.MDAMB231.vs.HCC70.boxplot.pdf',width=20,height=10)


############################################################ 
#### Fig 3b: left  panel, ranking cell lines according to median CNV correlation with Basal-like MET500 breast cancer samples
###########################################################
CNV.cor.matrix <- cor(CCLE.gene.cnv,MET500.breast.cancer.gene.cnv[,id],method='spearman')
m              <- apply(CNV.cor.matrix,1,median) %>% sort
names(m)       <- gsub(x=names(m),replacement = '',pattern = '_BREAST')
df             <- data.frame(rank=1:length(m),median.correlation=m)
df$color       <- 'black'
df[rownames(df) == 'MDAMB231','color'] <- 'red'
df$size       <- 1
df[rownames(df) == 'MDAMB231','size'] <- 4
ggplot(df) + geom_point(aes(x=rank,y=median.correlation,color=factor(color),size=size),show.legend=F) + ggplot.style + ylab('median.correlation')+ xlab('rank') + scale_color_manual(values=c('red'='red','black'='black')) + scale_size_area(max_size = 4)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/MET500.basal.like.sample.vs.CCLE.cnv.correlation.rank.plot.pdf',width=20,height=10)



############################################################ 
#### Fig 3a: right  panel, boxplot of expression correlation with Basal-like MET500 breast cancer samples, MDAMB231 and HCC70
###########################################################
expr.cor.matrix <- MET500.vs.CCLE.polyA.Basal.expression.correlation.result$correlation.matrix %>% t
m               <- apply(expr.cor.matrix,1,median) %>% sort
names(m)        <- gsub(x=names(m),replacement = '',pattern = '_BREAST')

df                     <- as.data.frame(expr.cor.matrix)
df$cell.line.name      <- gsub(x=rownames(df),replacement = '',pattern = '_BREAST')
long.df                <- melt(df,id.vars = 'cell.line.name')
long.df$cell.line.name <- factor(long.df$cell.line.name,levels = names(m))
long.df$color          <- 'black'
long.df[long.df$cell.line.name == 'MDAMB231','color'] <- 'red'
long.df <- long.df[long.df$cell.line.name %in% c('MDAMB231','HCC70'),]

ggplot(long.df) + geom_boxplot(aes(x=cell.line.name,y=value,color=factor(color)),size=1.5,show.legend=F) + ggplot.style + ylab('correlation') + xlab('') + scale_color_manual(values=c('red'='red','black'='black'))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/expr.correlation.MDAMB231.vs.HCC70.boxplot.pdf',width=20,height=10)


############################################################ 
#### Fig 3a: left  panel, ranking cell lines according to median expression correlation with Basal-like MET500 breast cancer samples
###########################################################
m              <- apply(expr.cor.matrix,1,median) %>% sort
names(m)       <- gsub(x=names(m),replacement = '',pattern = '_BREAST')
df             <- data.frame(rank=1:length(m),median.correlation=m)
df$color       <- 'black'
df[rownames(df) == 'MDAMB231','color'] <- 'red'
df$size       <- 1
df[rownames(df) == 'MDAMB231','size'] <- 4
ggplot(df) + geom_point(aes(x=rank,y=median.correlation,color=factor(color),size=size),show.legend=F) + ggplot.style + ylab('median.correlation')+ xlab('rank') + scale_color_manual(values=c('red'='red','black'='black')) + scale_size_area(max_size = 4)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/MET500.basal.like.sample.vs.CCLE.expr.correlation.rank.plot.pdf',width=20,height=10)




#####################################################################
#### Fig 3a: the sub-panel of the left panel, expression scatter-plot between MDAMB231 and with the "median" Basal-like samples
####################################################################
load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]

load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')
basal.sample <- MET500.vs.CCLE.polyA.Basal.expression.correlation.result$correlation.matrix %>% rownames
common.gene <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(MET500.log2.fpkm.matrix))
median.sample <- apply(MET500.log2.fpkm.matrix[common.gene,basal.sample],1,median)

df <- data.frame(median.sample = median.sample, 
                 MDAMB231 = CCLE.log2.rpkm.matrix[common.gene,'MDAMB231_BREAST'],
                 HCC70 = CCLE.log2.rpkm.matrix[common.gene,'HCC70_BREAST']
)

ggplot(df) + geom_point(aes(x=MDAMB231,y=median.sample),size=4) + geom_smooth(aes(x=MDAMB231,y=median.sample),method='loess')  + ggplot.style + xlab('MDAMB231') + ylab('MET500 Basal-like ') 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/MDAMB231.vs.basal.expr.scatter.plot.pdf',width=20,height=10)

### Well, I keep this plot, although it is not used in the manuscript
ggplot(df) + geom_point(aes(x=HCC70,   y=median.sample),size=4) + geom_smooth(aes(x=HCC70,y=median.sample),method='loess')     + ggplot.style + xlab('HCC70')    + ylab('MET500 Basal-like ')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/HCC70.vs.basal.expr.scatter.plot.pdf',width=20,height=10)






############################################################ 
#### Fig 3b: visulize mutation profiles of the 25 highly mutated genes in Basal-like MET500 samples
###########################################################
x <- cbind(MET500.breast.cancer.gene.mutation.profile[,id],CCLE.breast.cancer.gene.mutation.profile[,'MDAMB231_BREAST'])
s <- apply(MET500.breast.cancer.gene.mutation.profile[,id],1,sum) %>% sort(TRUE)
s <- s/length(id)
gg <- names(s)[s >= 0.1]
x  <- cbind(MET500.breast.cancer.gene.mutation.profile[gg,id],CCLE.breast.cancer.gene.mutation.profile[gg,'MDAMB231_BREAST'])

MET500.dist              <- dist(x ,method='binary')
MET500.hclust.rs         <- hclust(MET500.dist)
gg                       <- MET500.hclust.rs$labels[MET500.hclust.rs$order]
x                        <- x[gg,]

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h, gp = gpar(fill = "#CCCCCC", col = NA))
  },
  
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.8, gp = gpar(fill = "black", col = NA))
  }
)
col <- c("MUT" = "black",'background' = '#CCCCCC')

source.vec <- c( rep(times= id %>% length, x='MET500'),
                 rep(times= 1, x='MDAMB231_BREAST')
                 )
names(source.vec) <- c(id,'MDAMB231_BREAST')

library(circlize)
library(ComplexHeatmap)
pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/oncoprint.MDAMB231.pdf",width = 30,height=15 )
col.annotation  <-  HeatmapAnnotation(type=source.vec[c(id,'MDAMB231_BREAST')],which='column',col = list(type = c("MET500" =  "#800000", "MDAMB231_BREAST" = "#e6beff")) ,show_legend = FALSE)
oncoPrint(x, get_type = function(x) {ifelse(x == 1,'MUT','background')},row_names_side = 'left',
          show_heatmap_legend = FALSE,
          alter_fun = alter_fun, col = col, row_order = NULL, column_order=NULL, remove_empty_columns = T,show_row_names = TRUE,
          show_column_names = FALSE,show_row_barplot = F,top_annotation = col.annotation,show_pct = FALSE,row_names_gp = gpar(fontsize = 25,fontface=2)
) 
dev.off()



############################################################ 
#### Fig 3d Boxplot of expression correlation between cell lines (including lung-derived MDAMB231) and lung-metastais derived samples  #############
###########################################################

load('client-side/output/double.check.MDAMB231.R.output/double.check.MDAMB231.RData')

expr.cor.matrix <- microarray.lung.vs.CCLE.plus.lung.metastasis.MDAMB231$correlation.matrix %>% t
m               <- apply(expr.cor.matrix,1,median) %>% sort
names(m)        <- gsub(x=names(m),replacement = '',pattern = '_BREAST')

df                     <- as.data.frame(expr.cor.matrix)
df$cell.line.name      <- gsub(x=rownames(df),replacement = '',pattern = '_BREAST')
long.df                <- melt(df,id.vars = 'cell.line.name')
long.df$cell.line.name <- factor(long.df$cell.line.name,levels = names(m))
long.df$color          <- 'black'
long.df[long.df$cell.line.name %in% lung.derived.MDAMB231,'color'] <- 'red'

ggplot(long.df) + geom_boxplot(aes(x=cell.line.name,y=value,color=factor(color)),size=1.5,show.legend=F) + ggplot.style + ylab('correlation')+theme(axis.text.x=element_text(angle=90, hjust=1,size=15)) + xlab('') + scale_color_manual(values=c('red'='red','black'='black'))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/CCLE.and.MDAMB231.expr.cor.with.lung.samples.boxplot.pdf',width=20,height=10)



##############################################
#### FigS6 show the qqplot to demnostrated that we could use normal distribution to compute p-value
##############################################

load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')

rs <- qqnorm(MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor)
df <- data.frame(x=rs$x,y=rs$y)
ggplot(df) + geom_point(aes(x=x,y=y),size=4) + ggplot.style + xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +ylim(0,0.6)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/LumA.qqplot.pdf',width=20,height=10)


rs <- qqnorm(MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor)
df <- data.frame(x=rs$x,y=rs$y)
ggplot(df) + geom_point(aes(x=x,y=y),size=4) + ggplot.style + xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +ylim(0,0.6)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/LumB.qqplot.pdf',width=20,height=10)


rs <- qqnorm(MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor,main='LumB',cex.axis=2,font.axis=2)
df <- data.frame(x=rs$x,y=rs$y)
ggplot(df) + geom_point(aes(x=x,y=y),size=4) + ggplot.style + xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +ylim(0,0.6)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/Her2.qqplot.pdf',width=20,height=10)


rs <- qqnorm(MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor,main='LumB',cex.axis=2,font.axis=2)
df <- data.frame(x=rs$x,y=rs$y)
ggplot(df) + geom_point(aes(x=x,y=y),size=4) + ggplot.style + xlab('Theoretical Quantiles') + ylab('Sample Quantiles') +ylim(0,0.6)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/Basal.qqplot.pdf',width=20,height=10)




#############################################
### Fig S7a  BT483 shows significantly higer correlation  with MET500 LumB samples than MCF7
#############################################

cor.vec <- c(
            MET500.vs.CCLE.polyA.LumB.expression.correlation.result$correlation.matrix[,'MCF7_BREAST'],
            MET500.vs.CCLE.polyA.LumB.expression.correlation.result$correlation.matrix[,'BT483_BREAST']
           )
cell.line <- c(
               rep(x='MCF7',MET500.vs.CCLE.polyA.LumB.expression.correlation.result$correlation.matrix %>% nrow),
               rep(x='BT483',MET500.vs.CCLE.polyA.LumB.expression.correlation.result$correlation.matrix %>% nrow)
              )
df <- data.frame(cor = cor.vec,cell.line=cell.line)
ggplot(df) + geom_boxplot(aes(x=cell.line,y=cor),size=1.5) + ggplot.style + xlab('') + ylab('correlation') + ylim(c(0,0.6))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/MCF7.vs.BT483.LumB.boxplot.pdf',width=20,height=10)

rs <- wilcox.test(x=MET500.vs.CCLE.polyA.LumB.expression.correlation.result$correlation.matrix[,'MCF7_BREAST'],
            y=MET500.vs.CCLE.polyA.LumB.expression.correlation.result$correlation.matrix[,'BT483_BREAST'],
            paired=TRUE
            )
rs$p.value



###################################################
# Fig S7b  MDAMB415 does NOT show significantly higer correlation with MET500 LumA samples than T47D
###################################################
cor.vec <- c(
    MET500.vs.CCLE.polyA.LumA.expression.correlation.result$correlation.matrix[,'T47D_BREAST'],
    MET500.vs.CCLE.polyA.LumA.expression.correlation.result$correlation.matrix[,'MDAMB415_BREAST']
    )
cell.line <- c(
  rep(x='T47D',    MET500.vs.CCLE.polyA.LumA.expression.correlation.result$correlation.matrix %>% nrow),
  rep(x='MDAMB415',MET500.vs.CCLE.polyA.LumA.expression.correlation.result$correlation.matrix %>% nrow)
)
df <- data.frame(cor = cor.vec,cell.line=factor(x=cell.line,levels=c('MDAMB415','T47D')))
ggplot(df) + geom_boxplot(aes(x=cell.line,y=cor),size=1.5) + ggplot.style + xlab('') + ylab('correlation') + ylim(c(0,0.6))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/T47D.vs.MDAMB415.LumA.boxplot.pdf',width=20,height=10)

rs <- wilcox.test(x=MET500.vs.CCLE.polyA.LumA.expression.correlation.result$correlation.matrix[,'T47D_BREAST'],
                  y=MET500.vs.CCLE.polyA.LumA.expression.correlation.result$correlation.matrix[,'MDAMB415_BREAST'],
                  paired=TRUE
)
rs$p.value



###################################################
#### Fig S7c  EFM192A  show significantly higer correlation with MET500 Her2 samples than T47D
###################################################

cor.vec <- c(
    MET500.vs.CCLE.polyA.Her2.expression.correlation.result$correlation.matrix[,'T47D_BREAST'],
    MET500.vs.CCLE.polyA.Her2.expression.correlation.result$correlation.matrix[,'EFM192A_BREAST']
)
cell.line <- c(
    rep(x='T47D',   MET500.vs.CCLE.polyA.Her2.expression.correlation.result$correlation.matrix %>% nrow),
    rep(x='EFM192A',MET500.vs.CCLE.polyA.Her2.expression.correlation.result$correlation.matrix %>% nrow)
)
df <- data.frame(cor = cor.vec,cell.line=factor(x=cell.line,levels=c('EFM192A','T47D')))
ggplot(df) + geom_boxplot(aes(x=cell.line,y=cor),size=1.5) + ggplot.style + xlab('') + ylab('correlation') + ylim(c(0,0.6))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/T47D.vs.EFM192A.Her2.boxplot.pdf',width=20,height=10)

rs <- wilcox.test(x=MET500.vs.CCLE.polyA.Her2.expression.correlation.result$correlation.matrix[,'T47D_BREAST'],
                  y=MET500.vs.CCLE.polyA.Her2.expression.correlation.result$correlation.matrix[,'EFM192A_BREAST'],
                  paired=TRUE
)
rs$p.value


##############################################
#### Fig S8: consistency of sub-type specific correlation analysis between microarray and rna-seq
##############################################

#show the consistency of LumB expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
                 )
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/LumB.consistency.pdf',width=20,height=10)


#show the consistency of LumA expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/LumA.consistency.pdf',width=20,height=10)

#show the consistency of Her2 expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/Her2.consistency.pdf',width=20,height=10)

#show the consistency of Basal expression correlation analysis between microarray and rna-seq
df <- data.frame(MET500 = MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
                 GEO    = microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)
ggplot(df) + geom_point(aes(x=MET500,y=GEO),size=4) + ggplot.style + xlab('MET500') + ylab('microarray') 
cor(df$MET500,df$GEO,method='spearman')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section3/Basal.consistency.pdf',width=20,height=10)















