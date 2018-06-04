load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output/pick.out.cell.line.by.subtype.R.output//pick.out.cell.line.by.subtype.RData')
load('client-side/output/ssgsea.analysis.R.output//ssgsea.analysis.RData')
require(ComplexHeatmap)
require(foreach)

MET500.breast.cancer.polyA.non.Basal.sample           <- c(MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.Her2.sample)
MET500.breast.cancer.polyA.Basal.sample               <- c(MET500.breast.cancer.polyA.Basal.sample)
MET500.non.Basal.good.cell.line                       <- c(MET500.LumA.good.cell.line,MET500.LumB.good.cell.line,MET500.Her2.good.cell.line) %>% unique

#1. Wilcox rank test, p-value between MET500 samples and non.Basal good cell line, hallmark gene set
p.value.vec <- foreach(i = 1:nrow(hallmark.geneset.gsea.results),.combine='c') %do% {
    rs <- wilcox.test(hallmark.geneset.gsea.results[i,MET500.breast.cancer.polyA.non.Basal.sample],hallmark.geneset.gsea.results[i,MET500.non.Basal.good.cell.line])  
    rs$p.value
}
names(p.value.vec) <- rownames(hallmark.geneset.gsea.results)
fdr.vec            <- p.adjust(p.value.vec,method='fdr')
diff.activity.vec <- foreach(i  = 1:nrow(hallmark.geneset.gsea.results),.combine='c') %do% {
  median(hallmark.geneset.gsea.results[i,MET500.breast.cancer.polyA.non.Basal.sample]) - median(hallmark.geneset.gsea.results[i,MET500.non.Basal.good.cell.line]) 
}
names(diff.activity.vec) <- rownames(hallmark.geneset.gsea.results)
hallmark.df <- data.frame(gene.set = rownames(hallmark.geneset.gsea.results),
                          p.value  = p.value.vec,
                          fdr      = fdr.vec,
                          effect.size = diff.activity.vec
)
hallmark.df <- hallmark.df[order(hallmark.df$fdr),]
write.csv(file='/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Table/Section4/hallmark.test.result.csv',x=hallmark.df,quote=F,row.names=F)

##################################################
### Fig4: visualization of ssGSEA scores of hallmark genesets across TCGA,CCLE and MET500
#################################################

rownames(hallmark.geneset.gsea.results) <- gsub(rownames(hallmark.geneset.gsea.results),pattern = 'HALLMARK_',replacement = '')
TCGA.random.sample <- colnames(hallmark.geneset.gsea.results)[grepl(x=colnames(hallmark.geneset.gsea.results),pattern='TCGA')][1:37]
col.annotation  <-  HeatmapAnnotation(type=c(rep(times=MET500.breast.cancer.polyA.non.Basal.sample %>% length,x='MET500'),rep(times=TCGA.random.sample %>% length,x='TCGA'),rep(times=MET500.non.Basal.good.cell.line %>% length,x='CCLE')),
                                      which='column',
                                      col = list(type = c("MET500" =  "#800000", "TCGA" = "#008080","CCLE" = "#e6beff")) ,show_legend = FALSE
)

tmp                                     <- apply(hallmark.geneset.gsea.results,1,scale) %>% t
colnames(tmp)                           <- colnames(hallmark.geneset.gsea.results)
hallmark.geneset.gsea.results           <- tmp
names(fdr.vec)                          <- gsub(names(fdr.vec),pattern = 'HALLMARK_',replacement = '')

MET500.dist              <- dist(hallmark.geneset.gsea.results[,MET500.breast.cancer.polyA.non.Basal.sample] %>% t,method='euclidean')
MET500.hclust.rs         <- hclust(MET500.dist,method = 'complete')
MET500.breast.cancer.polyA.non.Basal.sample <- MET500.hclust.rs$labels[MET500.hclust.rs$order]

pdf("/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/hallmark.ssGSEA.heatmap.pdf",width = 30,height=15 )
Heatmap(hallmark.geneset.gsea.results[names(fdr.vec %>% sort),c(MET500.breast.cancer.polyA.non.Basal.sample,TCGA.random.sample,MET500.non.Basal.good.cell.line)],show_heatmap_legend = FALSE,  
        top_annotation = col.annotation,cluster_rows = FALSE,cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 15,fontface=2),
        width = unit(20, "in"),
        #col=colorRamp2(c(hallmark.geneset.gsea.results %>% min,hallmark.geneset.gsea.results %>% max),c('green','red')),
        #col=colorRamp2(c(-10,10),c('green','red')),
        col=colorRamp2(seq(from=-11,to=11,by=2),colorRampPalette(brewer.pal(11, "RdBu"))(12) %>% rev),
)
dev.off()

my_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(12) %>% rev
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
pdf(file  = '/Users/liuke/OneDrive/OneDrive - Michigan State University/Project/Cancer2CellLine/Manuscript/section.Figure/Section4/color.bar.pdf',width=20,height=20)
color.bar(my_palette,min = -11,max=11,nticks = 30) # draw the color bar
dev.off()














