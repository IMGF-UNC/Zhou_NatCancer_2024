runSeuratSignatures <- function(obj, 
                                ref.file = c('/datastore/nextgenout2/share/labs/imgf/ref/imgf_custom_references/imgf_mm_single_cell/imgf_mm_single_cell.tsv'),
                                colnames.reserved = c('NAME', 'SHORTNAME', 'DESCRIPTION', 'CITATION'),
                                name.grep = NULL,
                                print.graphics = FALSE,
                                heading.level = NA,
                                gene.search = FALSE
                                ) {
  
  # Setup our arguments
  #ref.file = match.arg(ref.file)
  
  ### Load reference file
  # We assign column names ourselves here, as read.table-flavor functions do not support files with a mix of named and un-named columns
  sigs.raw <- read.table(file = ref.file,
                         header = FALSE,
                         strip.white = T,
                         sep = "\t",
                         row.names=NULL)
  
  # Assign column names based on the header on 1st line
  colnames(sigs.raw)[match(colnames.reserved, sigs.raw[1,])] <- colnames.reserved
  sigs.raw <- sigs.raw[-c(1),] # Remove header line
  row.names(sigs.raw) <- 1:nrow(sigs.raw) # Fix our rownames

  # Give a warning if the reference signature name column is not unique
  if(! all(sigs.raw$NAME == unique(sigs.raw$NAME)) )
    warning(paste0("NAME column is not unique, but should be (reference file ", ref.file, ") "))
  
  # Ensure our names are valid 'name' strings in R, and force uniqueness
  sigs.raw$NAME <- make.names(sigs.raw$NAME, unique = TRUE)

  # Optionally, filter our signatures
  if( !is.null(name.grep) && !isFALSE(name.grep) && !is.na(name.grep)){
    if( is.character(name.grep) ) {
      sigs.raw <- sigs.raw[grepl(name.grep, sigs.raw$NAME) | grepl(name.grep, sigs.raw$SHORTNAME), ]
    } else if( is.vector(name.grep) || is.list(name.grep) ) {
      sigs.raw.new <- data.frame()
      for(x.name.grep in name.grep) {
        sigs.raw.new <- rbind(sigs.raw.new, sigs.raw[grepl(x.name.grep, sigs.raw$NAME) | grepl(x.name.grep, sigs.raw$SHORTNAME), ])
      }
      sigs.raw <- sigs.raw.new
    } else {
      stop("I don't know how to handle name.grep argument data type")
    }
    if( !nrow(sigs.raw) > 1 )
      stop("Error: name.grep argument returned zero rows")
  } # End option: name.grep
  
  ### Drop signatures with insufficient expression
  # Seurat will run signatures where genes are missing from expression data, but if all genes are missing then it will fail with an error.
  # We need to avoid this and run the remaining (working) signatures, so we pre-filter to remove signatures that we expect will fail when passed
  # to Seurat's AddModuleScore() function.
  for( i in 1:nrow(sigs.raw) ) {
    sigs.raw$GENESFOUND[i] <- length(which(sigs.raw[i,!colnames(sigs.raw) %in% colnames.reserved]
                       %in% rownames(obj)))
    if( sigs.raw$GENESFOUND[i] == 0 )
      message(paste0('Notice: Signature ', sigs.raw$NAME[i], ' has no genes identified in the dataset and cannot be calculated'))
  }
  # Drop signatures with zero genes found
  sigs.raw <- sigs.raw[sigs.raw$GENESFOUND > 0,]
  
  ### Convert gene signatures to a list (required by Seurat::AddModuleScore)
  sigs.features <- as.list(as.data.frame(t(
    sigs.raw[,colnames(sigs.raw)[!colnames(sigs.raw) %in% colnames.reserved]]
  )))
  sigs.name <- sigs.raw$NAME
  sigs.shortname <- sigs.raw$SHORTNAME
  sigs.description <- sigs.raw$DESCRIPTION
  sigs.citation <- sigs.raw$CITATION
  sigs.feature_name <- paste0('IMGF', 1:length(sigs.name) )
  sigs.genesfound <- sigs.raw$GENESFOUND


  
  
  ### Run the signature(s)  
  obj <- Seurat::AddModuleScore(obj, features = sigs.features, name = 'IMGF', assay = 'RNA', search = gene.search)
  
  if( print.graphics ) {
    # FeaturePlot(integrated, features = paste0(x.name, '1')) +
    #   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + theme_minimal() + NoLegend() + NoAxes() + labs(title = x.title)
    # VlnPlot(integrated, features = paste0(x.name, '1'), pt.size = 0.0 ) + theme_minimal() + NoLegend() + NoGrid() + labs(title = x.title)
    
    for(i in 1:length(sigs.name)) {
      sig.feature_name <- sigs.feature_name[i]
      sig.name <- sigs.name[i]
      sig.shortname <- sigs.shortname[i]
      sig.description <- sigs.description[i]
      sig.citation <- sigs.citation[i]
      sig.genesfound <- sigs.genesfound[i]
      
      if( is.numeric(heading.level) )
        cat(strrep('#', heading.level), sig.name, "\n\n")
      
      print(
        Seurat::FeaturePlot(obj, features = sig.feature_name, order = TRUE) + 
        ggplot2::scale_colour_gradientn(colours = c(RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))) + ggplot2::theme_minimal() + Seurat::NoLegend() + Seurat::NoAxes() +
          ggplot2::labs(title = sig.shortname, subtitle = paste0(sig.description, " (", sig.genesfound, " genes)"), caption = sig.citation),
      )
      
      print(
        Seurat::VlnPlot(obj, features = sig.feature_name, pt.size = 0.0) + ggplot2::theme_minimal() + Seurat::NoLegend() + Seurat::NoGrid() + 
          ggplot2::labs(title = sig.shortname, subtitle = sig.description, caption = sig.citation) + ggplot2::geom_boxplot(outlier.shape = 19)
      )

      if( is.numeric(heading.level) )
        cat("\n\n")
      
    }
  }
  return(obj)
}


imgfReintegrateObject <- function(obj, 
                                  opt.split.by = 'orig.ident',
                                  opt.reduction = 'rpca',
                                  opt.filter.cellnum = 150,
                                  opt.FindClusters.resolution = 0.8) {
  
  ### Filter samples by number of cells present
  obj.list <- SplitObject(obj, split.by = 'orig.ident')
  if( is.numeric(opt.filter.cellnum) ) {
    for(x.name in names(obj.list) ) {
      if(ncol(obj.list[[x.name]] ) < opt.filter.cellnum) {
        message(paste0("Removing sample ", x.name, " with less than ", opt.filter.cellnum, " cells from downstream analysis."))
        obj.list[[x.name]] <- NULL
      }
    }
  }
  
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = 'vst', nfeatures = 2000)
  })
  
  features <- SelectIntegrationFeatures(object.list = obj.list)
  obj.list <- lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  integration.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, reduction = opt.reduction)
  obj.rpca <- IntegrateData(anchorset = integration.anchors)
  
  DefaultAssay(obj.rpca) <- 'integrated'
  obj.rpca <- ScaleData(obj.rpca)
  obj.rpca <- RunPCA(obj.rpca, npcs = 50)
  obj.rpca <- RunUMAP(obj.rpca, reduction = 'pca', dims = 1:50)
  obj.rpca <- FindNeighbors(obj.rpca, reduction = 'pca', dims = 1:50)
  if( is.numeric)
  obj.rpca <- FindClusters(obj.rpca, resolution = 0.8)
  DimPlot(obj.rpca, label.box = T, label = T) + NoLegend()
}
