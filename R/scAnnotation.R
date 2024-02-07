#' @title scAnnotation
#' @param seurat_obj
#' @param annot.method dataset or marker
#' @param dataset_ref If annot.method = 'dataset', built-in
#' datasets will be used. Options are: `BlueprintEncode`,
#' `DICE`, `HPCA`, `ImmGen`, `MonacoImmune`, `MouseRNAseq`,
#' `Hematopoiesis`.
#' @param db.type If annot.method = 'marker',Options are:
#' `scType`, `CellMarker`, `CellTaxonomy`,`scMayoMap`, `VITA CytBase`.
#' Only `scType` not support Mouse.
#' @param Tissue
#' @param Species
#' @param Cancer
#' 
#' @importFrom SingleR SingleR
#' @import Seurat
#' 
#' @examples 
#' \dontrun{
#' library(scAnnotation)
#' result <- scAnnotation(example_data,annot.method='dataset')
#' result_2 <- scAnnotation(example_data,annot.method='marker')
#' }
#' 
#' @export
#' 
scAnnotation <- function(seurat_obj,
                         annot.method = 'marker',
                         datasetREF = 'BlueprintEncode',
                         db.type = 'CellMarker',
                         Tissue,
                         Species = 'Human',
                         Cancer = TRUE) {
  
  if (annot.method == 'dataset') {
    refData <- switch(datasetREF,
                      BlueprintEncode = scAnnoTest::BlueprintEncode,
                      DICE = scAnnoTest::DICE,
                      HPCA = scAnnoTest::HPCA,
                      ImmGen = scAnnoTest::ImmGen,
                      MonacoImmune = scAnnoTest::MonacoImmune,
                      MouseRNAseq = scAnnoTest::MouseRNAseq,
                      Hematopoiesis = scAnnoTest::Hematopoiesis
    )
    
    pred <- SingleR::SingleR(
      test = seurat_obj@assays$RNA@data,
      ref = refData,
      labels = refData$label.fine,
      clusters = seurat_obj$seurat_clusters
    )
    pred_df <- data.frame(pred)
    cell.annot <- data.frame(ClusterID = rownames(pred_df), SingleRAnnot = pred_df$labels)
    return(list(annot.result = cell.annot, score = pred_df))
    
  } else if (annot.method == 'marker') {
    cell.annot <- cellAnnot(
      seurat_obj = seurat_obj, db.type = db.type, Tissue = Tissue, 
      Species = Species, Cancer = Cancer
    )
  }
  return(cell.annot)
}
