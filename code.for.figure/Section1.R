require(plyr)
require(dplyr)
require(pheatmap)
require(ggplot2)
library(ComplexHeatmap)
require(foreach)
source('client-side/code/for.figure/ggplot.style.R')

load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')
load('~/Project/Cancer2CellLine/client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.cnv.RData')
load('~/Project/Cancer2CellLine/client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('~/Project/Cancer2CellLine/client-side/output/for.revision.nature.communications.round1.R.output//for.revision.nature.communications.round1.RData')

######################################################################################################################
#   
#  Fig 2
#  
######################################################################################################################




############################################################################################################
###### Fig 2a, Use R Oncoprint package to show the mutation porifles of hotspot mutated and DE mutated genes  
############################################################################################################

#### Compute differentially-mutated gene between MET500 and TCGA breast cancer samples
MET500.sample.number <- ncol(MET500.breast.cancer.gene.mutation.profile)
common.gene.list     <- intersect(names(TCGA.breast.cancer.gene.mutation.freq),names(MET500.breast.cancer.gene.mutation.freq)) # Well, this operation only removes gene ZNF668, whose mutation data is not returned by cBioportal
dm.df <- foreach(g= common.gene.list,.combine='rbind') %do% {
    q       <- MET500.breast.cancer.gene.mutation.profile[g,] %>% sum
    prob    <- TCGA.breast.cancer.gene.mutation.freq[g]    
    if(prob ==0) {
        prob <- min(TCGA.breast.cancer.gene.mutation.freq[TCGA.breast.cancer.gene.mutation.freq > 0])  
    }
    p.value <- pbinom(q=q,size=MET500.sample.number,prob = prob,lower.tail = F)#right side p-value,metastatic cancer has more mutation burdens
    data.frame(TCGA = prob,MET500=q/MET500.sample.number,CCLE=CCLE.breast.cancer.gene.mutation.freq[g],p.value=p.value,ratio=q/(MET500.sample.number * prob))
}
rownames(dm.df)  <- common.gene.list
dm.df$fdr        <- p.adjust(dm.df$p.value,method='fdr') 
dm.df            <- dm.df[order(dm.df$fdr),]
dm.gene          <- rownames(dm.df)[dm.df$fdr  < 0.001 ]


#### Check the overlap between highly mutated and differntially mutated genes, and pool them together
library(gplots)
venn(list(highly.mutated.gene=MET500.breast.cancer.top.mutated.gene,differentialy.mutated.gene=dm.gene))
pooled.gene              <- c(dm.gene,setdiff(MET500.breast.cancer.top.mutated.gene,dm.gene)) %>% unique

TCGA.mutation.cnt        <- apply(TCGA.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) 
TCGA.mutation.frequency  <- apply(TCGA.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) / ncol(TCGA.breast.cancer.gene.mutation.profile)

MET500.mutation.cnt        <- apply(MET500.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) 
MET500.mutation.frequency  <- apply(MET500.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) / ncol(MET500.breast.cancer.gene.mutation.profile)


CCLE.mutation.cnt        <- apply(CCLE.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) 
CCLE.mutation.frequency  <- apply(CCLE.breast.cancer.gene.mutation.profile[pooled.gene,],1,sum) / ncol(CCLE.breast.cancer.gene.mutation.profile)




### Rearrange MET500 and TCGA samples according to clustering results of mutation profile,rearrange CCLE samples according to sites (metastaic or primary) 
MET500.dist              <- dist(MET500.breast.cancer.gene.mutation.profile[pooled.gene,] %>% t,method='binary')
MET500.hclust.rs         <- hclust(MET500.dist)
MET500.rearranged.sample <- MET500.hclust.rs$labels[MET500.hclust.rs$order]

TCGA.dist                <- dist(TCGA.breast.cancer.gene.mutation.profile[pooled.gene,] %>% t,method='binary')
TCGA.hclust.rs           <- hclust(TCGA.dist)
TCGA.rearranged.sample   <- TCGA.hclust.rs$labels[TCGA.hclust.rs$order]

CCLE.rearranged.sample   <- c(CCLE.breast.cancer.metastatic.cell.line,CCLE.breast.cancer.non.metastatic.cell.line)


### Rearrange pooled gene according to clustering results
pooled.gene.dist         <- dist(cbind(MET500.breast.cancer.gene.mutation.freq[pooled.gene],CCLE.breast.cancer.gene.mutation.freq[pooled.gene]),method='euclidean')
pooled.gene.rs           <- hclust(pooled.gene.dist)
pooled.gene              <- pooled.gene.rs$labels[pooled.gene.rs$order]


#### Assemble the mutation profile matrix for the pooled genes
pooled.mutation.profile      <- cbind(MET500.breast.cancer.gene.mutation.profile[pooled.gene,MET500.rearranged.sample],
                                      CCLE.breast.cancer.gene.mutation.profile[pooled.gene,  CCLE.rearranged.sample],
                                      TCGA.breast.cancer.gene.mutation.profile[pooled.gene,  TCGA.rearranged.sample]
)

#### Show the oncoprint and save it
source.vec <- c( rep(times= MET500.rearranged.sample %>% length, x='MET500'),
                 rep(times= CCLE.rearranged.sample   %>% length, x='CCLE'),
                 rep(times= TCGA.rearranged.sample   %>% length, x='TCGA')
                 
)
names(source.vec) <- c(MET500.rearranged.sample,CCLE.rearranged.sample,TCGA.rearranged.sample)

alter_fun <- list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h, gp = gpar(fill = "#CCCCCC", col = NA))
    },
  
    MUT = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h*0.8, gp = gpar(fill = "black", col = NA))
    }
)
col <- c("MUT" = "black",'background' = '#CCCCCC')

library(circlize)
pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/oncoprint.MET500.CCLE.pdf",width = 30,height=15 )
col.annotation  <-  HeatmapAnnotation(type=source.vec[c(MET500.rearranged.sample,CCLE.rearranged.sample)],which='column',col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE)
oncoPrint(pooled.mutation.profile[,c(MET500.rearranged.sample,CCLE.rearranged.sample)], get_type = function(x) {ifelse(x == 1,'MUT','background')},row_names_side = 'left',
          show_heatmap_legend = FALSE,
          alter_fun = alter_fun, col = col, row_order = NULL, column_order=NULL, remove_empty_columns = T,
          show_column_names = FALSE,show_row_barplot = F,top_annotation = col.annotation,show_pct = FALSE,row_names_gp = gpar(fontsize = 13,fontface=2)
) + 
  Heatmap(cbind(MET500.breast.cancer.gene.mutation.freq[pooled.gene],CCLE.breast.cancer.gene.mutation.freq[pooled.gene]),show_column_dend = FALSE, 
          col=colorRamp2(c(0,1),c('blue','red')),
          row_names_gp = gpar(fontsize = 13,fontface=2),
          width = unit(4, "cm"),show_heatmap_legend = FALSE,
          top_annotation=HeatmapAnnotation( show_legend = FALSE,type=c('MET500','CCLE'),which='column',col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) )
  )
dev.off()


#### Draw the color bar
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
color.f    <- colorRamp2(c(0,1),c('blue','red'))
my_palette <- sapply(seq(from=0,to=1,by=0.05),color.f)
pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/color.bar.blue.to.red.from.0.to.1.pdf',width=20,height=20)
color.bar(my_palette,min = 0.0,max=1.0) 
dev.off()



############################################################################################################
###### Fig 2b, show analysis results of CNV profiles  
############################################################################################################


#### Fig 2b,first panel: the scatter plot of median.cnv between MET500 and CCLE cell lines 
gene.cnv.df <- data.frame(MET500= MET500.breast.cancer.gene.cnv.median,
                          TCGA  = TCGA.breast.cancer.gene.cnv.median ,
                          CCLE  = CCLE.breast.cancer.gene.cnv.median ,
                          CCLE.metastatic     = CCLE.breast.cancer.gene.cnv.median.metastatic ,
                          CCLE.non.metastatic = CCLE.breast.cancer.gene.cnv.median.non.metastatic 
                          )
sp.correaltion <-  cor(gene.cnv.df$MET500,gene.cnv.df$CCLE,method='spearman')
annot.text     <-  sprintf('Spearman.rank.cor=%s',round(sp.correaltion,3) %>% as.character)
ggplot(gene.cnv.df,aes(x=CCLE,y=MET500)) + geom_point(size=2.0,show.legend=F) +
  xlim(-0.8,1.0) + ylim(-0.25,0.8)+ xlab('CCLE') + ylab('MET500') +
  geom_abline(slope = 1,intercept = 0)  + ggplot.style 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/MET500.vs.CCLE.all.cnv.pdf',width=20,height=10)


#### Fig 2b, second and third panels: show cnv plot, metastatic.cell.line vs MET500, non.metastatic.cell.line vs MET500
flag <- gene.cnv.df$MET500 >= 0.4
gene.cnv.df$color         <- 'low.cnv'
gene.cnv.df[flag,'color'] <- 'high.cnv'


#### second panel
sp.correaltion <-  cor(gene.cnv.df$MET500,gene.cnv.df$CCLE.metastatic,method='spearman')
annot.text     <- sprintf('Spearman.rank.cor=%s',round(sp.correaltion,3) %>% as.character)
ggplot(gene.cnv.df,aes(x=CCLE.non.metastatic,y=MET500,col=color)) + geom_point(size=4.0,show.legend=F) +
  xlim(-1,1.5) + ylim(-1,1.5)+ xlab('CCLE (primary.site)') + ylab('MET500') +
  geom_abline(slope = 1,intercept = 0) +  ggplot.style +scale_color_manual(values=c('low.cnv' = 'black','high.cnv' = 'red')) 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/MET500.vs.CCLE.non.metastatic.cnv.pdf',width=20,height=10)


#### third panel
sp.correaltion <-  cor(gene.cnv.df$MET500,gene.cnv.df$CCLE.non.metastatic,method='spearman')
annot.text     <- sprintf('Spearman.rank.cor=%s',round(sp.correaltion,3) %>% as.character)
ggplot(gene.cnv.df,aes(x=CCLE.metastatic,y=MET500,col=color)) + geom_point(size=3.0,show.legend=F) + 
  xlim(-1,1.5) + ylim(-1,1.5)+xlab('CCLE (metastatic.site)') + ylab('MET500') +
  geom_abline(slope = 1,intercept = 0)  + ggplot.style +scale_color_manual(values=c('low.cnv' = 'black','high.cnv' = 'red')) 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/MET500.vs.CCLE.metastatic.cnv.pdf',width=20,height=10)


######################################################################################################################
#   
#  Fig S1
#  
######################################################################################################################


###########################################################
###### Fig S1a, Long-tailed mutation spectrum ###### 
###########################################################
qplot(MET500.breast.cancer.gene.mutation.freq, geom="histogram",binwidth=0.005) + ggplot.style + xlab('mutation frequency')
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/Long.tailed.mutation.spectrum.pdf',width=20,height=10) 




###########################################################
######Fig S1b,draw the valcano plot of differential mutation analysis between MET500 and TCGA
###########################################################
dm.df$is.de.gene <- ifelse(rownames(dm.df) %in% dm.gene,'Y','N')
draw.df          <- dm.df
ggplot(draw.df,aes(x=log2(ratio+1),y=-1 * log10(fdr),col=is.de.gene)) + geom_point(size=6.0,show.legend=F) + 
xlab('log2((MET500/TCGA)+1)') + ylab('-log10(fdr)') + xlim(0,6.5) + geom_hline(yintercept = 3,linetype="dashed") + 
ggplot.style +scale_color_manual(values=c('Y' = 'red','N' = 'grey')) 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/MET500.vs.TCGA.differential.mutation.frequency.volcano.plot.pdf',width=20,height=10)




###########################################################
####Compute number of genes whose mutational status could be resembled in at least 50% of CCLE cell lines
###########################################################
s         <- CCLE.breast.cancer.gene.mutation.freq[pooled.gene]
s[s>=0.5] # genes mutated in at least 50% of cell lines
s[s == 0] # genes NOT mutated in cell lines 
mut.cnt   <- apply(CCLE.gene.mutation.profile[pooled.gene,CCLE.breast.cancer.cell.line],2,sum) %>% sort #  compute the number of mutations carried by each cell line




###########################################################
###### Fig S1c, show mutation frequency of CCLE-cell-line-specific hypermutated genes in TCGA,MET500 and CCLE.  
###########################################################
CCLE.highly.mutated.gene           <-  names(CCLE.breast.cancer.gene.mutation.freq)[CCLE.breast.cancer.gene.mutation.freq >= 0.5]
CCLE.specific.highly.mutated.gene  <-  setdiff(CCLE.highly.mutated.gene,pooled.gene)
col.annotation                     <-  HeatmapAnnotation(type=c('MET500','TCGA','CCLE'),which='column',col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE)


require(circlize)
pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/CCLE.specific.highly.mutated.gene.mutation.frequency.heatmap.pdf",width = 30,height=15 )
Heatmap(log10(0.01+cbind(MET500.breast.cancer.gene.mutation.freq[CCLE.specific.highly.mutated.gene],TCGA.breast.cancer.gene.mutation.freq[CCLE.specific.highly.mutated.gene],CCLE.breast.cancer.gene.mutation.freq[CCLE.specific.highly.mutated.gene]) ),show_heatmap_legend = FALSE, 
        top_annotation = col.annotation,cluster_rows = FALSE,cluster_columns = FALSE,
        col=colorRamp2(c(-2,-0.05),c('blue','red')),row_names_gp = gpar(fontsize = 35,fontface=2),
        width = unit(20, "in")
)
dev.off()




###########################################################
###### Fig S1d: copy number of genes showing loss-of-copy number in CCLE cell lines were not altered in MET500
###########################################################
flag1 <- gene.cnv.df$CCLE > 0
flag2 <- gene.cnv.df$CCLE < 0


gain.and.loss.df <- rbind(   data.frame(cnv =gene.cnv.df$CCLE[flag1],  data.source ='CCLE',type='gain' ),
                             data.frame(cnv =gene.cnv.df$MET500[flag1],data.source ='MET500',type='gain'),
                             data.frame(cnv =gene.cnv.df$CCLE[flag2], data.source ='CCLE',type='loss' ),
                             data.frame(cnv =gene.cnv.df$MET500[flag2],data.source ='MET500',type='loss')
                             )
ggplot(gain.and.loss.df,aes(x=type,y=cnv,color=data.source)) +geom_boxplot(outlier.shape = NA,lwd=1.5) +  geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=data.source)) + ggplot.style + xlab('') + ylab('median cnv') + scale_color_manual(values=c('CCLE'="#e6beff",'MET500'='#800000'))
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/MET500.vs.CCLE.cnv.gain.and.loss.gene.boxplot.pdf',width = 20,height=10)
wilcox.test(gene.cnv.df$CCLE[flag1],gene.cnv.df$MET500[flag1],paired=TRUE)
wilcox.test(gene.cnv.df$CCLE[flag2],gene.cnv.df$MET500[flag2],paired=TRUE)




###########################################################
##### Fig S1e: show  cell lines derived from metastatic site more closely resemble genes with high CNV (in MET500 dataset)
###########################################################
flag <- gene.cnv.df$MET500 >= 0.4
a    <- (gene.cnv.df$MET500 - gene.cnv.df$CCLE.non.metastatic)[flag] %>% abs
b    <- (gene.cnv.df$MET500 - gene.cnv.df$CCLE.metastatic)[flag]     %>% abs
df   <- rbind(data.frame(abs.cnv.diff=a,  site ='primary site'),
              data.frame(abs.cnv.diff=b,  site ='metastatic site')
             )
ggplot(data=df,mapping = aes(x=site,y=abs.cnv.diff)) + geom_boxplot(outlier.shape = NA,lwd=1.5) + geom_point(size=4,position=position_jitterdodge(jitter.width = 0.2),aes(fill=site)) + ggplot.style + xlab('') + ylab('cnv.difference') 
ggsave(filename = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section1/abs.cnv.diff.metastatic.vs.non.metastatic.boxplot.pdf',width = 20,height=10)
wilcox.test(a,b,paired = TRUE)




##############################################
##### TableS1: mutation frequency of highly mutated and DE mutated genes
##### TableS2: just a copy of CCLE.cell.line.characteristic.csv under  /Project/Cancer2CellLine/client-side/meta.data/CCLE.cell.line.characteristic.csv
##############################################

table.s1.df <- data.frame( TCGA=TCGA.breast.cancer.gene.mutation.freq[pooled.gene],
                           MET500=MET500.breast.cancer.gene.mutation.freq[pooled.gene],
                           CCLE=CCLE.breast.cancer.gene.mutation.freq[pooled.gene]
                        )

table.s1.df <- data.frame( TCGA.mutation.cnt = TCGA.mutation.cnt,
                           TCGA.mutation.frequency = TCGA.mutation.frequency,
                           MET500.mutation.cnt = MET500.mutation.cnt,
                           MET500.mutation.frequency = MET500.mutation.frequency,
                           CCLE.mutation.cnt = CCLE.mutation.cnt,
                           CCLE.mutation.frequency = CCLE.mutation.frequency
                           
)

rownames(table.s1.df) <- pooled.gene
write.csv(file='/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Table/Section1/TableS1.csv',quote=F,
          x=table.s1.df
)


