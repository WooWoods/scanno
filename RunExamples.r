
if (!require('scAnnoTest')) {
  install.packages('scAnnoTest_1.0.0.tar.gz',repos = NULL, type = 'source')
}

suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(SCP))

# score visualization
ScorePlot <- function(AnnoRes,source = 'SingleR') {
  if (source == 'SingleR') {
    score <- AnnoRes$score
    score <- score[,c(1:(ncol(score)-5))]
    colnames(score) <- gsub('scores.','',colnames(score))
    score$cluster <- rownames(score)
    score <- reshape2::melt(score)
    score$cluster <- factor(score$cluster, levels = unique(score$cluster))
    colnames(score) <- c('Cluster','Label','Score')
  } else if (source == 'scType') {
    score <- AnnoRes$score[,c(1:3)]
    colnames(score) <- c('Cluster','Label','Score')
  }
  score_max <- score %>% dplyr::group_by(Cluster) %>% dplyr::top_n(n = 1, wt = Score)
  mypal <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,'Spectral'))(101))
  
  ggplot(data = score, aes(x = Cluster, y = Label, size = Score, fill = Score)) +
    geom_point(data = score, aes(x = Cluster, y = Label, size = Score, fill = Score),shape = 21) +
    scale_size_area(name = 'Score', max_size = 6, n.breaks = 4) +
    geom_point(data = score_max, aes(x = Cluster, y = Label, size = Score),
               shape = 8, color = 'red', stroke = 1.5) +
    guides(size = guide_legend(override.aes = list(fill = 'grey30',shape = 21),order = 1)) +
    scale_fill_gradientn(
      name = 'Score', limits = c(0,max(score$Score)),n.breaks = 3,
      colors = mypal,
      na.value = 'grey80',
      guide = guide_colorbar(frame.colour = 'black',ticks.colour = 'black',barheight = 4,barwidth = 1,order = 2)
    ) +
    scale_color_manual(values = NA, na.value = 'black') + 
    theme(
      panel.border = element_rect(size = 1,color = 'black',fill = NA),
      legend.position = 'right',
      legend.direction = "vertical",
      panel.grid.major = element_line(colour = "grey80", linetype = 2),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(lineheight = 0.5, hjust = 1,size = 12))
}

dat <- scAnnoTest::example_data
res_1 <- scAnnoTest::scAnnotation(
  seurat_obj = dat, annot.method = 'dataset',datasetREF = 'BlueprintEncode'
)
res_2 <- scAnnoTest::scAnnotation(
  seurat_obj = dat, annot.method = 'marker',
  Tissue = 'Lung', db.type = 'scType'
)
res_3 <- scAnnoTest::scAnnotation(
  seurat_obj = dat, annot.method = 'marker',Tissue = 'Lung',
  db.type = 'CellMarker',Species = 'Human',Cancer = T
)
res_4 <- scAnnoTest::scAnnotation(
  seurat_obj = dat, annot.method = 'marker',Tissue = 'Lung',
  db.type = 'CellTaxonomy',Species = 'Human'
)
res_5 <- scAnnoTest::scAnnotation(
  seurat_obj = dat, annot.method = 'marker',Tissue = 'Lung',
  db.type = 'scMayoMap',Species = 'Human'
)
res_6 <- scAnnoTest::scAnnotation(
  seurat_obj = dat, annot.method = 'marker',Tissue = 'Lung',
  db.type = 'VITA CytBase'
)

p1 <- ScorePlot(AnnoRes = res_1, source = 'SingleR')
p2 <- ScorePlot(AnnoRes = res_2, source = 'scType')
p3 <- ScorePlot(AnnoRes = res_3, source = 'scType')
p4 <- ScorePlot(AnnoRes = res_4, source = 'scType')
p5 <- ScorePlot(AnnoRes = res_5, source = 'scType')
p6 <- ScorePlot(AnnoRes = res_6, source = 'scType')

ggsave(filename = 'test/SingleR_score_bubble.pdf',p1,width = 7,height = 12)
ggsave(filename = 'test/scType_score_bubble.pdf',p2,width = 7,height = 6)
ggsave(filename = 'test/CellMarker_score_bubble.pdf',p3,width = 7,height = 7)
ggsave(filename = 'test/CellTaxonomy_score_bubble.pdf',p4,width = 7,height = 8)
ggsave(filename = 'test/scMayoMap_score_bubble.pdf',p5,width = 7,height = 8)
ggsave(filename = 'test/VITACytBase_score_bubble.pdf',p6,width = 7,height = 8)

dat$SingleR <- plyr::mapvalues(dat$seurat_clusters,from = result_1$annot.result$ClusterID,to = result_1$annot.result$SingleRAnnot)
annot <- cbind(result_2$annot.result['scTypeAnnot'],result_3$annot.result['CellMarkerAnnot'],
               result_4$annot.result['CellTaxonomyAnnot'],
               result_5$annot.result['scMayoMapAnnot'], result_6$annot.result['VITACytBaseAnnot'])
colnames(annot) <- gsub('Annot','',colnames(annot))
dat <- AddMetaData(dat,metadata = annot)

p1 <- CellDimPlot(dat, group.by = 'seurat_clusters', label = T, label.size = 5)
p2 <- CellDimPlot(dat, group.by = 'SingleR', label = T, label.size = 5)
p3 <- CellDimPlot(dat, group.by = 'scType', label = T, label.size = 5)
p4 <- CellDimPlot(dat, group.by = 'CellMarker', label = T, label.size = 5)
p5 <- CellDimPlot(dat, group.by = 'CellTaxonomy', label = T, label.size = 5)
p6 <- CellDimPlot(dat, group.by = 'scMayoMap', label = T, label.size = 5)
p7 <- CellDimPlot(dat, group.by = 'VITACytBase', label = T, label.size = 5)
ggsave(filename = 'test/seurat_clusters.pdf',p1,width = 8,height = 8)
ggsave(filename = 'test/SingleR.pdf',p2,width = 8,height = 8)
ggsave(filename = 'test/scType.pdf',p3,width = 8,height = 8)
ggsave(filename = 'test/CellMarker.pdf',p4,width = 8,height = 8)
ggsave(filename = 'test/CellTaxonomy.pdf',p5,width = 8,height = 8)
ggsave(filename = 'test/scMayo.pdf',p6,width = 8,height = 8)
ggsave(filename = 'test/VITACytBase.pdf',p7,width = 8,height = 8)