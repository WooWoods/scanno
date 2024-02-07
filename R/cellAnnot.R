#' @title cellAnnot
#' @param seurat_obj A Seurat object
#' @param db.type scType, CellMarker, CellTaxonomy, scMayoMap or VITA CytBase.
#' @param Species Human or Mouse.
#' @param Tissue 
#' @param Cancer TRUE or FALSE
#' 
#' @import Seurat
#' @importFrom dplyr `%>%` group_by top_n
#' @importFrom plyr mapvalues
#' 
#' @export
#' 
cellAnnot <- function(seurat_obj, 
                      db.type = 'CellMarker',
                      Tissue,
                      Species = 'Human', 
                      Cancer = TRUE) {
  
  if (db.type == 'scType') {
    gs <- prepareGeneset(db.type = db.type,Tissue = Tissue)
  } else if (db.type == 'CellMarker') {
    gs <- prepareGeneset(
      db.type = db.type, Species = Species,
      Tissue = Tissue, Cancer = Cancer
    )
  } else if (db.type == 'CellTaxonomy') {
    gs <- prepareGeneset(db.type = db.type, Species = Species,Tissue = Tissue)
  } else if (db.type == 'VITA CytBase') {
    gs <- prepareGeneset(db.type = db.type)
  } else if (db.type == 'scMayoMap') {
    gs <- prepareGeneset(db.type = db.type, Tissue = Tissue)
  }
  
  DefaultAssay(seurat_obj) <- 'RNA'
  seurat_obj <- ScaleData(seurat_obj)
  es <- Scoring(expr = seurat_obj@assays$RNA@scale.data, gs = gs)
  result <- do.call('rbind',lapply(unique(seurat_obj$seurat_clusters),function(x){
    es.max.cl <- sort(rowSums(es[,rownames(seurat_obj@meta.data[seurat_obj$seurat_clusters == x,])]),decreasing = TRUE)
    head(data.frame(cluster = x, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj$seurat_clusters==x)),10)
  }))
  score <- result %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 3, wt = scores)
  top1 <- result %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 1, wt = scores)
  
  if (db.type == 'scType') {
    seurat_obj$scTypeAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = top1$cluster, to = top1$type)
    annot.result <- seurat_obj@meta.data[,c('seurat_clusters','scTypeAnnot')]
  } else if (db.type == 'CellMarker') {
    seurat_obj$CellMarkerAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = top1$cluster, to = top1$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','CellMarkerAnnot')]
  } else if (db.type == 'CellTaxonomy') {
    seurat_obj$CellTaxonomyAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = top1$cluster, to = top1$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','CellTaxonomyAnnot')]
  } else if (db.type == 'VITA CytBase') {
    seurat_obj$VITACytBaseAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = top1$cluster, to = top1$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','VITACytBaseAnnot')]
  } else if (db.type == 'scMayoMap') {
    seurat_obj$scMayoMapAnnot <- plyr::mapvalues(seurat_obj$seurat_clusters, from = top1$cluster, to = top1$type)
    annot.result <- seurat_obj@meta.data[, c('seurat_clusters','scMayoMapAnnot')]
  }
  
  return(list(annot.result = annot.result, score = score))
}
