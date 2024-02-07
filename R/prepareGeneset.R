#' @title prepareGeneset
#' @param db.type `scType`, `CellMarker`, `CellTaxonomy`,`scMayoMap`, or `VITA CytBase`.
#' @param Species Human or Mouse
#' @param Tissue
#' @param Cancer TRUE or FALSE
#' 
#' @export
#' 
prepareGeneset <- function(db.type,Tissue,Species='Human',Cancer=TRUE) {
  
  if (db.type == 'scType') {
    message('[[', Sys.time(), ']]: --- scTypeDB currently support Species: Human')
    message('[[', Sys.time(), ']]: --- scTypeDB currently support tissue: ')
    ScTypeDB_full <- scAnnoTest::ScTypeDB_full
    print(unique(ScTypeDB_full$tissueType))
    print(paste0('Select Tissue: ', Tissue))
    cellmarker <- ScTypeDB_full[grep(Tissue, ScTypeDB_full$tissueType, ignore.case = TRUE),]
    cellmarker$geneSymbolmore1 <- gsub(' ','',cellmarker$geneSymbolmore1)
    cellmarker$geneSymbolmore1 <- sapply(1:nrow(cellmarker), function(i){
      markers_all <- gsub(' ', '', unlist(strsplit(cellmarker$geneSymbolmore1[i],',')))
      markers_all <- markers_all[markers_all != 'NA' & markers_all != '']
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    cellmarker$geneSymbolmore1 <- gsub('///',',',cellmarker$geneSymbolmore1)
    cellmarker$geneSymbolmore1 <- gsub(' ', '', cellmarker$geneSymbolmore1)
    gs <- lapply(1:nrow(cellmarker), function(j){
      gsub(' ','',unlist(strsplit(toString(cellmarker$geneSymbolmore1[j]),',')))
    })
    names(gs) <- cellmarker$cellName
    
  } else if (db.type == 'CellMarker') {
    message('[[', Sys.time(), ']]: --- CellMarker currently support Species: ')
    Cell_marker_All <- scAnnoTest::Cell_marker_All
    print(unique(Cell_marker_All$species))
    print(paste0('Select Species: ', Species))
    cellmarker <- subset(Cell_marker_All, species == Species)
    message('[[', Sys.time(), ']]: --- CellMarker currently support tissue: ')
    print(unique(cellmarker$tissue_type))
    print(paste0('Select Tissue: ', Tissue))
    cellmarker <- cellmarker[grep(Tissue,cellmarker$tissue_type, ignore.case = T),]
    message('[[', Sys.time(), ']]: --- CellMarker currently support Cells: ')
    print(unique(cellmarker$cell_type))
    if (Cancer == TRUE) {
      print(paste0('Select Cells: Cancer cell'))
      cellmarker <- subset(cellmarker, cell_type == 'Cancer cell')
    } else {
      print(paste0('Select Cells: Normal cell'))
      cellmarker <- subset(cellmarker, cell_type == 'Normal cell')
    }
    celltype <- unique(cellmarker$cell_name)
    gs <- lapply(1:length(celltype), function(i){
      markers_all <- cellmarker[cellmarker$cell_name == celltype[i],]$Symbol
      markers_all <- as.character(na.omit(markers_all))
      markers_all <- sort(markers_all)
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    gs <- lapply(1:length(gs), function(j){
      gsub(' ', '', unlist(strsplit(toString(gs[j]),',')))
    })
    names(gs) <- celltype
    
  } else if (db.type == 'CellTaxonomy') {
    message('[[', Sys.time(), ']]: --- CellTaxonomy currently support Species: ')
    Cell_Taxonomy <- scAnnoTest::Cell_Taxonomy
    print(unique(Cell_Taxonomy$species))
    print(paste0('Select Species: ', Species))
    cellmarkers <- subset(Cell_Taxonomy, species == Species)
    message('[[', Sys.time(), ']]: --- CellTaxonomy currently support tissue: ')
    print(unique(cellmarkers$Tissue_standard))
    print(paste0('Select Tissue: ', Tissue))
    cellmarkers <- cellmarkers[grep(Tissue, cellmarkers$Tissue_standard, ignore.case = T),]
    
    celltype <- unique(cellmarkers$Cell_standard)
    gs <- lapply(1:length(celltype), function(i){
      markers_all <- cellmarkers[cellmarkers$Cell_standard == celltype[i],]$Cell_Marker
      markers_all <- as.character(na.omit(markers_all))
      markers_all <- sort(markers_all)
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    gs <- lapply(1:length(gs), function(j){
      gsub(' ', '', unlist(strsplit(toString(gs[j]),',')))
    })
    names(gs) <- celltype
    
  } else if (db.type == 'VITA CytBase') {
    VITACytBase <- scAnnoTest::VITACytBase
    if (Species=='Human') {
      VITACytBase <- VITACytBase[,c('Markers','Celltype')]
    } else if (Species=='Mouse') {
      VITACytBase <- VITACytBase[,c('MouseOrthologs','Celltype')]
      colnames(VITACytBase) <- c('Markers','Celltype')
    }
    celltype <- unique(VITACytBase$Celltype)
    gs <- lapply(1:length(celltype), function(i){
      markers_all <- VITACytBase[VITACytBase$Celltype == celltype[i],]$Markers
      markers_all <- as.character(na.omit(markers_all))
      markers_all <- sort(markers_all)
      markers_all <- unique(na.omit(markers_all))
      paste0(markers_all, collapse = ',')
    })
    gs <- lapply(1:length(gs), function(j){
      gsub(' ', '', unlist(strsplit(toString(gs[j]),',')))
    })
    names(gs) <- celltype
  } else if (db.type == 'scMayoMap') {
    message('[[', Sys.time(), ']]: --- scMayoMap currently support Species: Human')
    message('[[', Sys.time(), ']]: --- scMayoMap currently support tissue: ')
    
    scMayoMap <- scAnnoTest::scMayoMap
    print(paste0('Select Tissue: ', Tissue))
    cellmarker <- scMayoMap[grep(Tissue, scMayoMap$tissue, ignore.case = T),]
    cellmarker$markers <- gsub(' ','',cellmarker$markers)
    cellmarker$marker <- sapply(1:nrow(cellmarker),function(i){
      markers_all <- gsub(' ','',unlist(strsplit(cellmarker$marker[i],',')))
      markers_all <- markers_all[markers_all != 'NA' & markers_all != '']
      if (length(markers_all) > 0) {
        markers_all <- unique(na.omit(markers_all))
        paste0(markers_all, collapse = ',')
      } else {
        ''
      }
    })
    cellmarker$marker <- gsub('///',',',cellmarker$marker)
    cellmarker$marker <- gsub(' ','',cellmarker$marker)
    gs <- lapply(1:nrow(cellmarker), function(j){
      gsub(' ','',unlist(strsplit(toString(cellmarker$marker[j]),',')))
    })
    names(gs) <- cellmarker$celltype
  }
  
  return(gs)
}
