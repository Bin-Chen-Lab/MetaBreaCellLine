require(ggplot2)
require(pheatmap)


load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData')

## Draw the mutation heatmap plot. Within each group,samples were clusterd according to the mutation profile 
pooled.top.mutated.gene      <- c(MET500.breast.cancer.top.mutated.gene,TCGA.breast.cancer.top.mutated.gene) %>% unique
pooled.top.mutated.gene.freq <- MET500.breast.cancer.gene.mutation.freq[pooled.top.mutated.gene] %>% sort(decreasing = TRUE)
pooled.top.mutated.gene      <- names(pooled.top.mutated.gene.freq)

MET500.dist         <- dist(MET500.breast.cancer.gene.mutation.profile[,pooled.top.mutated.gene],method='binary')
MET500.hclust.rs    <- hclust(MET500.dist)
MET500.rearranged.sample <- MET500.hclust.rs$labels[MET500.hclust.rs$order]
TCGA.dist           <- dist(TCGA.breast.cancer.gene.mutation.profile[,pooled.top.mutated.gene],method='binary')
TCGA.hclust.rs      <- hclust(TCGA.dist)
TCGA.rearranged.sample <- TCGA.hclust.rs$labels[TCGA.hclust.rs$order]
CCLE.dist           <- dist(CCLE.breast.cancer.gene.mutation.profile[,pooled.top.mutated.gene],method='binary')
CCLE.hclust.rs      <- hclust(CCLE.dist)
CCLE.rearranged.sample   <- CCLE.hclust.rs$labels[CCLE.hclust.rs$order]

pooled.mutation.profile      <- rbind(MET500.breast.cancer.gene.mutation.profile[MET500.rearranged.sample,pooled.top.mutated.gene],
                                      TCGA.breast.cancer.gene.mutation.profile[TCGA.rearranged.sample,pooled.top.mutated.gene],
                                      CCLE.breast.cancer.gene.mutation.profile[CCLE.rearranged.sample,pooled.top.mutated.gene]
                                      )

source.vec <- c( rep(times=c(MET500.breast.cancer.gene.mutation.profile %>% nrow),x='MET500'),
                 rep(times=c(TCGA.breast.cancer.gene.mutation.profile %>% nrow),  x='TCGA'),
                 rep(times=c(CCLE.breast.cancer.gene.mutation.profile %>% nrow),  x='CCLE')
               )
name.vec <- c( c(MET500.breast.cancer.gene.mutation.profile %>% rownames),
               c(TCGA.breast.cancer.gene.mutation.profile %>% rownames),
               c(CCLE.breast.cancer.gene.mutation.profile %>% rownames)
              )
col.annotation           <- data.frame(sample.source=as.factor(source.vec))
rownames(col.annotation) <- name.vec
pheatmap(t(pooled.mutation.profile),cluster_rows = F,cluster_cols = F,show_rownames = T,show_colnames = F,legend = F,annotation_col = col.annotation)


####Draw the paired scatter plot of mutation frequency
gene.mutation.freq.df <- data.frame(MET500=MET500.breast.cancer.gene.mutation.freq[pooled.top.mutated.gene],
                                    TCGA  =TCGA.breast.cancer.gene.mutation.freq[pooled.top.mutated.gene],
                                    CCLE  =CCLE.breast.cancer.gene.mutation.freq[pooled.top.mutated.gene]
                                   )
ggplot(gene.mutation.freq.df,aes(x=CCLE,y=MET500)) + geom_point(size=4) + xlim(0,1) + ylim(0,1)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('MET500 vs CCLE')
ggplot(gene.mutation.freq.df,aes(x=CCLE,y=TCGA))   + geom_point(size=4) + xlim(0,1) + ylim(0,1)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('TCGA vs CCLE')
ggplot(gene.mutation.freq.df,aes(x=TCGA,y=MET500)) + geom_point(size=4) + xlim(0,1) + ylim(0,1)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('MET500 vs TCGA')


#####Draw the paired scatter plot of cnv#############
### Well, CCLE have outliergenes whose cnv is < -3 or >1, so I generate two versions (with or without the outlier)
load('client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.cnv.RData')
gene.cnv.df <- data.frame(MET500=MET500.breast.cancer.gene.cnv.median,
                          TCGA  =TCGA.breast.cancer.gene.cnv.median,
                          CCLE  =CCLE.breast.cancer.gene.cnv.median
                          )
ggplot(gene.cnv.df,aes(x=CCLE,y=MET500)) + geom_point(size=4) + xlim(-0.8,0.8) + ylim(-0.8,0.8)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('MET500 vs CCLE')
ggplot(gene.cnv.df,aes(x=CCLE,y=MET500)) + geom_point(size=4) + xlim(-4,4) + ylim(-1,1)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('MET500 vs CCLE')

ggplot(gene.cnv.df,aes(x=CCLE,y=TCGA)) + geom_point(size=4) + xlim(-0.8,0.8) + ylim(-0.8,0.8)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('TCGA vs CCLE')
ggplot(gene.cnv.df,aes(x=CCLE,y=TCGA)) + geom_point(size=4) + xlim(-4,4) + ylim(-1,1)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('TCGA vs CCLE')

ggplot(gene.cnv.df,aes(x=TCGA,y=MET500)) + geom_point(size=4) + xlim(-0.8,0.8) + ylim(-0.8,0.8)+geom_abline(slope = 1,intercept = 0) + theme_grey(base_size = 45)+theme(axis.text=element_text(face="bold"),axis.title=element_text(face="bold"))+ggtitle('MET500 vs TCGA')







