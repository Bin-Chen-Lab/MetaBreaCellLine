load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output/organize.breast.cancer.organoid.data.R.output/organoid.RData')
source('client-side/code/util.R')
source('client-side/code/for.figure/ggplot.style.R')
require(dplyr)


# Pick out most varied genes across CCLE cell lines
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]


# Correlating cell lines and organoids to MET500 samples (different subtypes) 
s                   <- intersect(rownames(CCLE.log2.rpkm.matrix),rownames(organoid.log2.rpkm.matrix))
merged.model.matrix <- cbind(CCLE.log2.rpkm.matrix[s,],organoid.log2.rpkm.matrix[s,])

MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumB.sample],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Basal.sample], expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.Her2.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.Her2.sample],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,MET500.breast.cancer.polyA.LumA.sample],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)
MET500.vs.CCLE.and.organoid.polyA.non.Basal.expression.correlation.result  <- pick.out.cell.line(expr.of.samples = MET500.log2.fpkm.matrix[,c(MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.Her2.sample)],  expr.of.cell.lines = merged.model.matrix,CCLE.rna.seq.marker.gene.1000)


organoid.line <- colnames(organoid.log2.rpkm.matrix)



######################################################################################################################
#   
#  Fig 5
#  
######################################################################################################################

### Fig 5a, visualize of organoid ranking results between different subtypes 
require('PerformanceAnalytics')
source('client-side/code/for.figure/chart.Correlation.R') # I steal the code from PerformanceAnalytics package and modify it,remove the weired stars

mat <- cbind( MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[organoid.line], 
              MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result$cell.line.median.cor[organoid.line],
              MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result$cell.line.median.cor[organoid.line],
              MET500.vs.CCLE.and.organoid.polyA.Her2.expression.correlation.result$cell.line.median.cor[organoid.line]
)
colnames(mat) <- c('Basal-like','LuminalB','LuminalA','Her2\nenriched')

pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/organoid.result.correlation.between.subtype.pdf',width=30,height=30)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.labels=10,cor.range=c(-1,1),cex.cor.tuning = 0.7)
dev.off()


### Fig 5b, right panel: boxplot to show that statisticaly organoid population shows higher correlation with MET500 samples
#Basal-like organoid vs Basal-like cell line
df1 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[Basal.organoid],
                  type = rep(x='organoid',times = Basal.organoid %>% length)
                  )

df2 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.Basal.cell.line],
                  type = rep(x='cell.line',times = CCLE.breast.cancer.Basal.cell.line %>% length)
)
ggplot(rbind(df1,df2),aes(x=type,y=median.cor)) + geom_boxplot(outlier.shape=NA,lwd=1.2) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=type)) + ggplot.style + xlab('') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/Basal.organoid.vs.cell.line.population.pdf',width=20,height=15)
wilcox.test(df1$median.cor,df2$median.cor,paired = FALSE)


### Fig 5b, left panel: boxplot to show that statisticaly organoid population shows higher correlation with MET500 samples
#non-Basal-like organoid vs non-Basal-like cell line
non.Basal.organoid                     <- c(LumA.organoid,LumB.organoid,Her2.organoid)
CCLE.breast.cancer.non.Basal.cell.line <- c(CCLE.breast.cancer.LumA.cell.line,CCLE.breast.cancer.LumB.cell.line,CCLE.breast.cancer.Her2.cell.line)

df1 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.non.Basal.expression.correlation.result$cell.line.median.cor[non.Basal.organoid],
                  type = rep(x='organoid',times = non.Basal.organoid %>% length)
)

df2 <- data.frame(median.cor = MET500.vs.CCLE.and.organoid.polyA.non.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.non.Basal.cell.line],
                  type = rep(x='cell.line',times = CCLE.breast.cancer.non.Basal.cell.line %>% length)
)
ggplot(rbind(df1,df2),aes(x=type,y=median.cor)) + geom_boxplot(outlier.shape=NA,lwd=1.2) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=type)) + ggplot.style + xlab('') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/non.Basal.organoid.vs.cell.line.population.pdf',width=20,height=15)
wilcox.test(df1$median.cor,df2$median.cor,paired = FALSE)




# Fig 5c, for each of the subtypes show that most correlated organoid is better than most correlated cell line
#LumA
cell.line <- 'MDAMB415_BREAST'
organoid  <- 'MMC01031_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.LumA.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='MMC01031',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='MDAMB415',times = nrow(rs$correlation.matrix))
)
ggplot(rbind(df1,df2),aes(x=type,y=correlation)) + geom_boxplot(outlier.shape=NA,lwd=1.2) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=type)) + ggplot.style + xlab('') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/LumA.pdf',width=20,height=15)

wilcox.test(df1$correlation,df2$correlation,paired = TRUE)


#LumB
cell.line <- 'BT483_BREAST'
organoid  <- 'MMC01031_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.LumB.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='MMC01031',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='BT483',times = nrow(rs$correlation.matrix))
)
ggplot(rbind(df1,df2),aes(x=type,y=correlation)) + geom_boxplot(outlier.shape=NA,lwd=1.2) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=type)) + ggplot.style + xlab('') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/LumB.pdf',width=20,height=15)

wilcox.test(df1$correlation,df2$correlation,paired = TRUE)


#Her2
cell.line <- 'EFM192A_BREAST'
organoid  <- 'MMC01031_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.Her2.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='MMC01031',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='EFM192A',times = nrow(rs$correlation.matrix))
)
ggplot(rbind(df1,df2),aes(x=type,y=correlation)) + geom_boxplot(outlier.shape=NA,lwd=1.2) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=type)) + ggplot.style + xlab('') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/Her2.pdf',width=20,height=15)

wilcox.test(df1$correlation,df2$correlation,paired = TRUE)


#Basal
cell.line <- 'HCC70_BREAST'
organoid  <- 'W1009_BREAST_ORGANOID'
rs        <- MET500.vs.CCLE.and.organoid.polyA.Basal.expression.correlation.result

df1 <- data.frame(correlation = rs$correlation.matrix[,organoid],
                  type = rep(x='W1009',times = nrow(rs$correlation.matrix))
)

df2 <- data.frame(correlation = rs$correlation.matrix[,cell.line],
                  type = rep(x='HCC70',times = nrow(rs$correlation.matrix))
)
ggplot(rbind(df1,df2),aes(x=type,y=correlation)) + geom_boxplot(outlier.shape=NA,lwd=1.2) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=type)) + ggplot.style + xlab('') + ylab('transcriptome similarity')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/Basal.pdf',width=20,height=15)
wilcox.test(df1$correlation,df2$correlation,paired = TRUE)






