load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output/pick.out.cell.line.by.subtype.R.output//pick.out.cell.line.by.subtype.RData')
load('client-side/output/ssgsea.analysis.R.output//ssgsea.analysis.RData')
require(ComplexHeatmap)
require(foreach)
require(plyr)
require(dplyr)
require(RColorBrewer)
require(circlize)
source('client-side/code/for.figure/ggplot.style.R')



######################################################################################################################
#   
#  Fig 6
#  
######################################################################################################################



### Fig 6a : heatmap to visulaize ssGSEA results of non-basal MET500, organoid, cell.line
MET500.number    <- sum(grepl(x=colnames(non.basal.gsea.results),pattern='SRR'))
organoid.number  <- sum(grepl(x=colnames(non.basal.gsea.results),pattern='ORGANOID'))
cell.line.number <- ncol(non.basal.gsea.results) - MET500.number - organoid.number


col.annotation  <-  HeatmapAnnotation(type=c(rep(times=MET500.number,x='MET500'),
                                             rep(times=organoid.number,x='organoid'),
                                             rep(times=cell.line.number,x='CCLE')
                                             ),
                                      which='column',
                                      col = list(type = c("MET500" =  "#800000", "organoid" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE
)


tmp           <- apply(non.basal.gsea.results,1,scale) %>% t
rownames(tmp) <- gsub(rownames(tmp),pattern = 'HALLMARK_',replacement = '')

pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/non.Basal.hallmark.ssGSEA.heatmap.pdf",width = 30,height=15 )
Heatmap(tmp,show_heatmap_legend = FALSE,  
        top_annotation = col.annotation,cluster_rows = TRUE,cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 15,fontface=2),
        width = unit(20, "in"),
        col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev),
        clustering_distance_rows = 'pearson',
        clustering_distance_columns = 'pearson'
        
)
dev.off()

# color bar of heatmap
my_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev
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
pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/color.bar.pdf',width=20,height=20)
color.bar(my_palette,min = -5,max=5,nticks = 30) # draw the color bar
dev.off()


### Fig 6b,left panel: scatter plot to show DA anlysis results (non-basal-like subtype)
fdr.df <- as.data.frame(non.basal.aov.fdr)
ggplot(fdr.df) + geom_point(aes(x=organoid,y=cell.line),size=4) + ggplot.style + xlim(0,10) + ylim(0,10) +geom_vline(xintercept=3,linetype=2) + geom_hline(yintercept=3,linetype=2) + xlab('MET500 vs ORGANOIDS') +ylab('MET500 vs CCLE') + geom_abline(intercept = 0, slope = 1,color='red')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/non.basal.fdr.scatter.plot.pdf',width=20,height=20)


### Fig 6b,right panel: scatter plot to show DA anlysis results (basal-like subtype)
fdr.df <- as.data.frame(basal.aov.fdr)
ggplot(fdr.df) + geom_point(aes(x=organoid,y=cell.line),size=4) + ggplot.style + xlim(0,10) + ylim(0,10) +geom_vline(xintercept=3,linetype=2) + geom_hline(yintercept=3,linetype=2) + xlab('MET500 vs ORGANOIDS') +ylab('MET500 vs CCLE') + geom_abline(intercept = 0, slope = 1,color='red')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/basal.fdr.scatter.plot.pdf',width=20,height=20)




###############
#
# Fig S10
#
#### ssGSEA heatmap of basal MET500, organoid, cell.line
MET500.number    <- sum(grepl(x=colnames(basal.gsea.results),pattern='SRR'))
organoid.number  <- sum(grepl(x=colnames(basal.gsea.results),pattern='ORGANOID'))
cell.line.number <- ncol(basal.gsea.results) - MET500.number - organoid.number


col.annotation  <-  HeatmapAnnotation(type=c(rep(times=MET500.number,x='MET500'),
                                             rep(times=organoid.number,x='organoid'),
                                             rep(times=cell.line.number,x='CCLE')
),
which='column',
col = list(type = c("MET500" =  "#800000", "organoid" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE
)


tmp           <- apply(basal.gsea.results,1,scale) %>% t
rownames(tmp) <- gsub(rownames(tmp),pattern = 'HALLMARK_',replacement = '')

pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/Basal.hallmark.ssGSEA.heatmap.pdf",width = 30,height=15 )
Heatmap(tmp,show_heatmap_legend = FALSE,  
        top_annotation = col.annotation,cluster_rows = TRUE,cluster_columns = TRUE,
        row_names_gp = gpar(fontsize = 15,fontface=2),
        width = unit(20, "in"),
        col=colorRamp2(seq(from=-5,to=5,by=0.5),colorRampPalette(brewer.pal(11, "RdBu"))(21) %>% rev)
)
dev.off()







MET500.breast.cancer.polyA.non.Basal.sample           <- c(MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.Her2.sample)
MET500.breast.cancer.polyA.Basal.sample               <- c(MET500.breast.cancer.polyA.Basal.sample)
MET500.non.Basal.good.cell.line                       <- c(MET500.LumA.good.cell.line,MET500.LumB.good.cell.line,MET500.Her2.good.cell.line) %>% unique


#### determine gene set class for each subtype#######
# class I:    only significant in organoids vs MET500 comparision
# class II:   only significant in CCLE vs MET500 comparision
# class III:  significant in both  comparisions
# class IV:   not significant in either comparision

a           <- as.data.frame(non.basal.aov.fdr)
b           <- as.data.frame(basal.aov.fdr)
colnames(a) <- c('non.basal.organoids','non.basal.cell.line')
colnames(b) <- c('basal.organoids','basal.cell.line')

d <- cbind(a,b)
d$non.basal.organoids <- ifelse(d$non.basal.organoids > 3,1,0)
d$non.basal.cell.line <- ifelse(d$non.basal.cell.line > 3,1,0)
d$basal.organoids     <- ifelse(d$basal.organoids > 3,1,0)
d$basal.cell.line     <- ifelse(d$basal.cell.line > 3,1,0)
m.d                   <- as.matrix(d)

non.basal.m.d <- m.d[,1:2]
row.sum       <- apply(non.basal.m.d,1,sum) 
f4 <- which(row.sum == 0) # class VI ,   NOT significant in either comparision
f3 <- which(row.sum == 2) # class III ,  significant in both
f1 <- which(non.basal.m.d[,1] == 1 & non.basal.m.d[,2] == 0 ) # class I , only significant in organoids
f2 <- which(non.basal.m.d[,1] == 0 & non.basal.m.d[,2] == 1 ) # class II , only significant in cell lines
require(pheatmap)
pheatmap(non.basal.m.d[c(f1,f2,f3,f4),],cluster_rows = F,cluster_cols = F,cellwidth = 150,border_color = 'black')
non.basal.c1 <- rownames(non.basal.m.d)[f1]
non.basal.c2 <- rownames(non.basal.m.d)[f2]
non.basal.c3 <- rownames(non.basal.m.d)[f3]
non.basal.c4 <- rownames(non.basal.m.d)[f4]

basal.m.d <- m.d[,3:4]
row.sum       <- apply(basal.m.d,1,sum) 
f4 <- which(row.sum == 0)
f3 <- which(row.sum == 2)
f1 <- which(basal.m.d[,1] == 1 & basal.m.d[,2] == 0 )
f2 <- which(basal.m.d[,1] == 0 & basal.m.d[,2] == 1 )
pheatmap(basal.m.d[c(f1,f2,f3,f4),],cluster_rows = F,cluster_cols = F,cellwidth = 150,border_color = 'black')
basal.c1 <- rownames(basal.m.d)[f1]
basal.c2 <- rownames(basal.m.d)[f2]
basal.c3 <- rownames(basal.m.d)[f3]
basal.c4 <- rownames(basal.m.d)[f4]

require( gplots)
venn(list(non.basal.c1=non.basal.c1,basal.c1=basal.c1))
venn(list(non.basal.c2=non.basal.c2,basal.c2=basal.c2))
venn(list(non.basal.c3=non.basal.c3,basal.c3=basal.c3))
venn(list(non.basal.c4=non.basal.c4,basal.c4=basal.c4))



# Fig 6c
gs1 <- 'HALLMARK_ANDROGEN_RESPONSE'
gs2 <- 'HALLMARK_E2F_TARGETS'
gs3 <- 'HALLMARK_COMPLEMENT'
gs4 <- 'HALLMARK_FATTY_ACID_METABOLISM'


MET500.number    <- sum(grepl(x=colnames(non.basal.gsea.results),pattern='SRR'))
organoid.number  <- sum(grepl(x=colnames(non.basal.gsea.results),pattern='ORGANOID'))
cell.line.number <- ncol(non.basal.gsea.results) - MET500.number - organoid.number

data.source <- c(rep(times=MET500.number,x='MET500'),
                 rep(times=organoid.number,x='ORGANOIDS'),
                 rep(times=cell.line.number,x='CCLE')
)
df1 <- data.frame(gsea.score = non.basal.gsea.results[gs1,],
                  data.source = data.source,
                  gene.set    = gs1
)
df2 <- data.frame(gsea.score = non.basal.gsea.results[gs2,],
                  data.source = data.source,
                  gene.set    = gs2
)
df3 <- data.frame(gsea.score = non.basal.gsea.results[gs3,],
                  data.source = data.source,
                  gene.set    = gs3
)
df4 <- data.frame(gsea.score = non.basal.gsea.results[gs4,],
                  data.source = data.source,
                  gene.set    = gs4
)
df             <- rbind(df1,df2,df3,df4)
df$data.source <- factor(df$data.source,levels = c('MET500','ORGANOIDS','CCLE'))
df$subtype     <- 'non.Basal'
non.basal.df   <- df


MET500.number    <- sum(grepl(x=colnames(basal.gsea.results),pattern='SRR'))
organoid.number  <- sum(grepl(x=colnames(basal.gsea.results),pattern='ORGANOID'))
cell.line.number <- ncol(basal.gsea.results) - MET500.number - organoid.number

data.source <- c(rep(times=MET500.number,x='MET500'),
                 rep(times=organoid.number,x='ORGANOIDS'),
                 rep(times=cell.line.number,x='CCLE')
)
df1 <- data.frame(gsea.score = basal.gsea.results[gs1,],
                  data.source = data.source,
                  gene.set    = gs1
)
df2 <- data.frame(gsea.score = basal.gsea.results[gs2,],
                  data.source = data.source,
                  gene.set    = gs2
)
df3 <- data.frame(gsea.score = basal.gsea.results[gs3,],
                  data.source = data.source,
                  gene.set    = gs3
)
df4 <- data.frame(gsea.score = basal.gsea.results[gs4,],
                  data.source = data.source,
                  gene.set    = gs4
)
df             <- rbind(df1,df2,df3,df4)
df$data.source <- factor(df$data.source,levels = c('MET500','ORGANOIDS','CCLE'))
df$subtype     <- 'Basal'
basal.df       <- df


draw.df             <- rbind(basal.df,non.basal.df)
draw.df$data.source <- factor(draw.df$data.source,levels = c('ORGANOIDS','MET500','CCLE'))
draw.df$subtype     <- factor(draw.df$subtype,levels = c('non.Basal','Basal'))

ggplot(data = draw.df[draw.df$gene.set == gs1,],aes(x=subtype,color=data.source,y=gsea.score)) + geom_boxplot(outlier.shape = NA) +  geom_jitter(position=position_jitterdodge(jitter.width = 0.2),size=4) + ggplot.style + scale_color_manual(values=c("MET500" =  "#800000", "ORGANOIDS" = "#008080","CCLE" = "#e6beff")) + xlab('') + ylab('ssGSEA.score') + ggtitle(gs1)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/gs1.box.plot.pdf',width=20,height=20)

ggplot(data = draw.df[draw.df$gene.set == gs2,],aes(x=subtype,color=data.source,y=gsea.score)) + geom_boxplot(outlier.shape = NA) +  geom_jitter(position=position_jitterdodge(jitter.width = 0.2),size=4) + ggplot.style + scale_color_manual(values=c("MET500" =  "#800000", "ORGANOIDS" = "#008080","CCLE" = "#e6beff")) + xlab('') + ylab('ssGSEA.score') + ggtitle(gs2)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/gs2.box.plot.pdf',width=20,height=20)


ggplot(data = draw.df[draw.df$gene.set == gs3,],aes(x=subtype,color=data.source,y=gsea.score)) + geom_boxplot(outlier.shape = NA) +  geom_jitter(position=position_jitterdodge(jitter.width = 0.2),size=4) + ggplot.style + scale_color_manual(values=c("MET500" =  "#800000", "ORGANOIDS" = "#008080","CCLE" = "#e6beff")) + xlab('') + ylab('ssGSEA.score') + ggtitle(gs3)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/gs3.box.plot.pdf',width=20,height=20)

ggplot(data = draw.df[draw.df$gene.set == gs4,],aes(x=subtype,color=data.source,y=gsea.score)) + geom_boxplot(outlier.shape = NA) +  geom_jitter(position=position_jitterdodge(jitter.width = 0.2),size=4) + ggplot.style + scale_color_manual(values=c("MET500" =  "#800000", "ORGANOIDS" = "#008080","CCLE" = "#e6beff")) + xlab('') + ylab('ssGSEA.score') + ggtitle(gs4)
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section5/gs4.box.plot.pdf',width=20,height=20)
